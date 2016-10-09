#!/usr/bin/env python2
"""
Script that generates data.py
"""

import csv
import hashlib
import json
import math
import os
import re
import sys
import urllib2
from collections import OrderedDict

if __name__ == '__main__' and __package__ is None:
    sys.path.insert(1,
          os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import dimple.cell

WIKI_URL = ('https://raw.githubusercontent.com/wiki/'
            'ccp4/dimple/Crystallization-Contaminants.md')

CLUSTER_CUTOFF = 0.03

CACHE_DIR = 'cached'

def cached_urlopen(url, cache_name=None):
    if not cache_name:
        cache_name = hashlib.sha1('abc').hexdigest()[:12]
    elif cache_name == -1:
        cache_name = url.split('/')[-1]
    path = os.path.join(CACHE_DIR, cache_name)
    if not os.path.exists(path):
        if not os.path.isdir(CACHE_DIR):
            os.mkdir(CACHE_DIR)
        print '--> %s' % path
        response = urllib2.urlopen(url)
        with open(path, 'wb') as f:
            f.write(response.read())
    return open(path)

def extract_pdb_ids(page):
    pdb_ids = []
    for line in page:
        if ' - ' in line:
            _, ids = line.rsplit(' - ', 1)
            for pdb_id in ids.split(','):
                pdb_id = pdb_id.strip().upper()
                assert len(pdb_id) == 4, pdb_id
                assert pdb_id[0].isdigit(), pdb_id
                pdb_ids.append(pdb_id)
    return pdb_ids

def extract_uniprot_names(page):
    pat = re.compile(r'`[A-Z0-9]+_[A-Z0-9]+`')
    names = []
    for line in page:
        obj = pat.search(line)
        if obj:
            names.append(obj.group(0).strip('`'))
        elif '`' in line:
            sys.stderr.write('No UniProt name? ' + line)
    return names

def fetch_pdb_info_from_ebi(pdb_id):
    pdb_id = pdb_id.lower()
    # We tried to remove heteromers when quering Uniprot, but we still have
    # things like 4TTD (Hetero 3-mer with a single link to UniProt AC).
    # An assembly that is a monomer can also be classified as "hetero", e.g.
    # http://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/2war
    # From what I understand such forms are unlikely as contaminants?
    # The same about "homo" monomeric assemblies that contain 2 distinct
    # polypeptide molecules (e.g. 5AQ9) or polypeptide and DNA or RNA.
    entry_url = "http://www.ebi.ac.uk/pdbe/api/pdb/entry/"
    response = cached_urlopen(entry_url + 'summary/' + pdb_id,
                              'ebi-summary-%s.json' % pdb_id)
    summary = json.load(response)[pdb_id]
    assert isinstance(summary, list)
    assert len(summary) == 1
    summary = summary[0]
    exper_method = summary['experimental_method']
    if exper_method != ['X-ray diffraction']:
        if exper_method != ['X-ray powder diffraction']:
            print 'WARNING: got %s in %s' % (exper_method, pdb_id)
        return None
    forms = set(a['form'] for a in summary['assemblies'])
    assert forms.issubset({'homo', 'hetero'}) # some are both (2EKS)
    if forms != {'homo'}:
        return None
    entities = summary['number_of_entities']
    big_molecules = ['polypeptide', 'dna', 'dna/rna', 'rna']
    if sum(entities[x] for x in big_molecules) > 1:
        return None
    # now http://www.ebi.ac.uk/pdbe/api/pdb/entry/experiment/:pdbid
    response = cached_urlopen(entry_url + 'experiment/' + pdb_id,
                              'ebi-exper-%s.json' % pdb_id)
    summary = json.load(response)[pdb_id]
    assert isinstance(summary, list)
    assert len(summary) == 1
    summary = summary[0]
    parameters = [summary['cell'][x] for x in 'a b c alpha beta gamma'.split()]
    sg = summary['spacegroup']
    # TODO: standarize unit cell settings (e.g. 3wai, 1y6e)
    cell = dimple.cell.Cell(tuple(parameters), symmetry=sg)
    cell.pdb_id = pdb_id
    # Estimate the model quality (cell.quality). Higher is better.
    # Use simplistic criterium, similar to the one from
    # http://www.rcsb.org/pdb/statistics/clusterStatistics.do
    try:
        cell.quality = 1.0 / summary['resolution'] - summary['r_free']
    except (ValueError, TypeError):
        cell.quality = -10
    if summary['experiment_data_available'] != 'Y':
        cell.quality = -20
    # TODO: unit-cell axes in non-standard order can be confusing
    if False:
        cell.quality = -30
    return cell

def dump_uniprot_tab(pdb_ids):
    query_url = ('http://www.uniprot.org/uniprot/'
                 '?query=database:(type:pdb%%20%s)&sort=score&format=tab')
    pdb_ids = [x.strip(',') for x in sys.stdin.read().split()]
    assert all(len(x) == 4 for x in pdb_ids)
    print >>sys.stderr, 'Read %d IDs' % len(pdb_ids)
    for pdb_id in pdb_ids:
        response = cached_urlopen(query_url % pdb_id, 'up-pdb-%s.tab' % pdb_id)
        text = response.read()
        print pdb_id
        print text
        print

def uniprot_names_to_acs(names):
    query_url = ('http://www.uniprot.org/uniprot/'
                 '?query=%s&columns=id,entry%%20name&format=tab')
    acs = OrderedDict()
    for name in names:
        response = cached_urlopen(query_url % name, 'up-%s-ie.tab' % name)
        lines = response.readlines()
        assert len(lines) >= 2
        ac, entry_name = lines[1].strip().split('\t')
        assert entry_name == name
        acs[ac] = name
    return acs

def fetch_uniref_clusters(acs):
    query_url = ('http://www.uniprot.org/uniref/'
                 '?query=%s&fil=identity:1.0&format=tab')
    clusters = OrderedDict()
    for ac in acs:
        response = cached_urlopen(query_url % ac, 'ur100-%s.tab' % ac)
        reader = csv.DictReader(response, delimiter='\t')
        for d in reader:
            key = d['Cluster ID']
            name = d['Cluster name']
            if name.startswith('Cluster: '):
                name = name[9:]
            members = [a.strip() for a in d['Cluster members'].split(';')]
            if ac not in members:
                continue
            print '%s -> %s (%s) - %d members, %d residues' % (
                    ac, key, name, len(members), int(d['Length']))
            if key in clusters:
                sys.exit('Duplicated cluster: ' + key)
            clusters[key] = members
    return clusters

def read_pdbtosp():
    f = open('pdbtosp.txt')
    f = cached_urlopen('http://www.uniprot.org/docs/pdbtosp.txt', -1)
    def get_acs_from_line(line):
        ac1 = line[40:52]
        assert ac1[0] == '(', line
        assert ac1[-1] in ' )', line
        acs = [ac1.strip(' ()')]
        if len(line) > 66:
            ac2 = line[66:78]
            assert ac2[0] == '(', line
            assert ac2[-1] in ' )', line
            acs.append(ac2.strip(' ()'))
        return acs
    ret = {}
    for line in f:
        if len(line) > 25 and line[0].isdigit():
            pdbid = line[0:4]
            assert pdbid not in ret
            method = line[6:14].strip()
            acs = get_acs_from_line(line)
            if line[78:79] == ',':
                acs += get_acs_from_line(f.next())
            ret[pdbid] = (method, acs)
    return ret


def group_homomeric_crystals_by_ac(pdbtosp):
    hmap = {}
    for key, (method, acs) in pdbtosp.iteritems():
        if method == 'X-ray' and len(acs) == 1:
            hmap.setdefault(acs[0], []).append(key)
    return hmap

def get_pdb_clusters(cells):
    # simplistic clustering
    # TODO: use scipy instead
    # scipy.spatial.distance.pdist(lambda a, b: cell_distance(a, b))
    # scipy.cluster.hierarchy.linkage
    cc = [x for x in cells[::-1] if x is not None]
    cc.sort(key=lambda x: (x.symmetry, x.a))
    while cc:
        clu = [cc.pop()]
        for n in range(len(cc)-1, -1, -1):
            dist = cell_distance(cc[n], clu[0])
            #print ':::', cc[n], clu[0], dist
            if dist < CLUSTER_CUTOFF:
                clu.append(cc.pop(n))
        yield clu

def cell_distance(a, b):
    if a.symmetry != b.symmetry:
        return float('inf')
    return a.max_shift_in_mapping(b)

# Representative entry is picked as a one with at least median "quality"
# metric that is least dissimilar to other cells.
def get_representative_unit_cell(cells):
    if len(cells) < 2:
        return cells[0]
    cells.sort(key=lambda x: -x.quality)
    cutoff = int(math.ceil(0.5 * len(cells)))
    def med_metric(c):
        return max(cell_distance(c, other) for other in cells)
    best = min(cells[:cutoff], key=med_metric)
    #print '===>', med_metric(best)
    return best

def write_data_py(representants):
    with open('data.py', 'w') as data_py:
        data_py.write('\nDATA = [\n')
        for cell in representants:
            data = ['"%s"' % cell.pdb_id,
                    '"%s"' % cell.symmetry
                   ] + ['%6.2f' % x for x in cell.cell]
            data_py.write('(' + ', '.join(data) + '),\n')
        data_py.write(']')

def main():
    page = cached_urlopen(WIKI_URL, -1).readlines()
    pdb_ids_from_page = extract_pdb_ids(page)
    pdbtosp = read_pdbtosp()
    missing = [p for p in pdb_ids_from_page if p not in pdbtosp]
    print 'missing:', missing

    uniprot_names = extract_uniprot_names(page)
    acs = uniprot_names_to_acs(uniprot_names)
    clusters = fetch_uniref_clusters(acs)
    #clusters = {'UniRef100_P00698': clusters['UniRef100_P00698']}
    homomers = group_homomeric_crystals_by_ac(pdbtosp)
    representants = []
    for name, clust in clusters.items():
        pdb_set = sum([homomers[c] for c in clust if c in homomers], [])
        up_name = ' '.join(acs[a] for a in clust if a in acs)
        print '%s -> %s (%d entries) -> %d PDBs' % (
                up_name, name, len(clust), len(pdb_set))
        cells = [fetch_pdb_info_from_ebi(pdb_id) for pdb_id in pdb_set]
        for pdb_cluster in get_pdb_clusters(cells):
            r = get_representative_unit_cell(pdb_cluster)
            print '%s %-10s %s  %.2f' % (r.pdb_id, r.symmetry, r, r.quality)
            representants.append(r)
    representants.sort(key=lambda x: x.a)
    write_data_py(representants)



if __name__ == '__main__':
    main()
