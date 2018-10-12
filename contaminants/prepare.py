#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Script that generates data.py - list of unit cells of contaminants
"""

import csv
import hashlib
import itertools
import json
import math
import os
import re
import sys
import urllib2
from collections import OrderedDict
import xml.etree.ElementTree as ET
import numpy

SCIPY_CLUSTERING = True

if SCIPY_CLUSTERING:
    import scipy.cluster
    import scipy.spatial

if __name__ == '__main__' and __package__ is None:
    assert os.path.basename(os.getcwd()) == 'contaminants'
    sys.path.insert(1,
          os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import dimple.cell

WIKI_URL = ('https://raw.githubusercontent.com/wiki/'
            'ccp4/dimple/Crystallization-Contaminants.md')

CONTAMINER_JSON = 'https://strube.cbrc.kaust.edu.sa/contaminer/contabase.json'

CLUSTER_CUTOFF = 0.04

OUTPUT_PY_FILE = 'data.py'
OUTPUT_JSON_FILE = 'data.json'
CACHE_DIR = 'cached'

# ContaBase has homologs of these
NOT_IN_CONTABASE = ['POLG_HRV2', 'THRB_HUMAN', 'FA10_BOVIN']

def cached_urlopen(url, cache_name):
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

def parse_wiki_page(page):
    pat = re.compile(r'`[A-Z0-9]+_[A-Z0-9]+`')
    data = OrderedDict()
    for line in page:
        obj = pat.search(line)
        if obj:
            name = obj.group(0).strip('`')
            if name == 'LIKE_THIS':
                continue
        elif '`' in line:
            sys.exit('No UniProt name? ' + line)
        else:
            continue
        sp = line[obj.span()[1]:].rsplit(' - ', 1)
        desc = sp[0].strip(' ①② ')
        pdb_ids = []
        if len(sp) == 2:
            pdb_ids = [i.strip().upper() for i in sp[1].split(',')]
            assert all(len(i) == 4 for i in pdb_ids), pdb_ids
            assert all(i[0].isdigit() for i in pdb_ids), pdb_ids
        data[name] = (desc, pdb_ids)
    return data

def fetch_contaminer_protein_acs():
    response = cached_urlopen(CONTAMINER_JSON, -1)
    data = json.load(response)
    all_acs = {}
    for category in data['categories']:
        for c in category['contaminants']:
            all_acs[c['uniprot_id']] = c['short_name'] + ' ' + c['long_name']
    return all_acs

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
    single_crystal_methods = ['X-ray diffraction', 'Neutron Diffraction']
    #other_method = ['X-ray powder diffraction', 'Solution NMR',
    #                'Electron Microscopy', 'Electron crystallography']
    if not any(m in single_crystal_methods for m in exper_method):
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
    summaries = json.load(response)[pdb_id]
    assert isinstance(summaries, list)
    summary = summaries[0]
    par_names = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
    parameters = [summary['cell'][x] for x in par_names]
    sg = summary['spacegroup']
    for alt_sum in summaries[1:]:
        assert [alt_sum['cell'][x] for x in par_names] == parameters
        assert alt_sum['spacegroup'] == sg
    if sg == 'A 1': # pesky 1LKS
        return None
    # Correct spacegroup, but in non-standard settings.
    # We could convert it to the reference/standard settings.
    if sg in ['P 1 1 21']: # only this one for now
        return None

    # Standarize unit cell settings (e.g. 3wai, 1y6e).
    # It'd be better to use a library function (cctbx?) for this.
    # to_standard() below makes a<b<c in primitive orthorhombic cells.
    if sg in ['P 1 2 1', 'P 1 21 1'] and parameters[0] > parameters[2]:
        parameters[0], parameters[2] = parameters[2], parameters[0]

    cell = dimple.cell.Cell(tuple(parameters), symmetry=sg).to_standard()

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
        query = query_url % name
        response = cached_urlopen(query, 'up-%s-ie.tab' % name)
        lines = response.readlines()
        assert len(lines) >= 2, 'up-%s-ie.tab' % name
        for line in lines[1:]:
            ac, entry_name = line.strip().split('\t')
            if entry_name == name:
                acs[ac] = name
                break
        else:
            sys.exit(name + ' not in the results of:\n' + query)
    return acs

def fetch_uniref_clusters(acs, verbose=False):
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
            if ac not in members:  # querying P63165 also returns P63165-2
                continue
            #if key != 'UniRef100_P02931': continue
            if verbose:
                print '%s -> %s (%s) - %d members, %d residues' % (
                        ac, key, name, len(members), int(d['Length']))
            if key in clusters:
                sys.exit('Duplicated cluster: ' + key)
            clusters[key] = members
        if reader.line_num == 0:
            print 'Not found in UniRef:', ac
            clusters['_' + ac] = [ac]  # we still want to keep this AC
    return clusters

def read_current_pdb_entries():
    f = cached_urlopen('http://www.rcsb.org/pdb/rest/getCurrent',
                       'current-rcsb.xml')
    root = ET.parse(f).getroot()
    entries = set(child.attrib['structureId'] for child in root)
    print 'Current PDB entries: %d' % len(entries)
    return entries

def read_sifts_mapping():
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/csv/uniprot_pdb.csv'
    f = cached_urlopen(url, -1)
    ret = {}
    for line in f:
        if not line or line[0] == '#' or line.startswith('SP_PRIMARY,'):
            continue
        u, pp = line.strip().split(',')
        for p in pp.split(';'):
            ret.setdefault(p.upper(), []).append(u)
    return ret

def get_pdb_clusters(cells):
    ordered = sorted(cells, key=lambda x: (x.symmetry, x.a))
    for _, same_sg in itertools.groupby(ordered, key=lambda x: x.symmetry):
        same_sg = list(same_sg)
        if len(same_sg) == 1:
            yield same_sg
        elif SCIPY_CLUSTERING:
            for cluster in cluster_by_cell_size_scipy(same_sg):
                yield cluster
        else:
            for cluster in cluster_by_cell_size(same_sg):
                yield cluster

def cluster_by_cell_size(cells):
    # simplistic clustering
    while cells:
        clu = [cells.pop()]
        for n in range(len(cells)-1, -1, -1):
            if cell_distance(cells[n], clu[0]) < CLUSTER_CUTOFF:
                clu.append(cells.pop(n))
        yield clu

def cluster_by_cell_size_scipy(cells):
    n = len(cells)
    # matrix like from scipy.spatial.distance.pdist(..., metric=func)
    dist_matrix = numpy.zeros(n * (n - 1) // 2, dtype=numpy.double)
    k = 0
    for i in xrange(0, n-1):
        for j in xrange(i+1, n):
            dist_matrix[k] = cells[i].max_shift_in_mapping(cells[j])
            k += 1
    assert k == len(dist_matrix)
    linkage = scipy.cluster.hierarchy.linkage(dist_matrix, method='complete')
    r = scipy.cluster.hierarchy.fcluster(linkage, CLUSTER_CUTOFF,
                                         criterion='distance')
    for cluster_number in range(1, max(r)+1):
        clu = [c for n, c in enumerate(cells) if r[n] == cluster_number]
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

def write_output_file(representants):
    with open(OUTPUT_PY_FILE, 'w') as out:
        out.write('\n# This file was automatically generated by the '
                  'contaminants/prepare.py script.')
        out.write('\nDATA = [\n')
        for cell in representants:
            out.write('("%s", "%s", ' % (cell.pdb_id, cell.symmetry))
            out.write('%.2f, %.2f, %.2f, ' % (cell.a, cell.b, cell.c))
            for angle in (cell.alpha, cell.beta, cell.gamma):
                if angle in (90.0, 120.0):
                    out.write('%.0f., ' % angle)
                else:
                    out.write('%.2f, ' % angle)
            out.write('"%s"), # %s\n' % (cell.uniprot_name, cell.comment))
        out.write(']')

def write_json_file(uniprot_acs, extra_info, representants):
    groups = {}
    for r in representants:
        group = groups.setdefault(r.uniprot_name, {})
        ref = r.to_reference()
        leaf_name = '%s  (%s)' % (r.pdb_id.upper(), ref.parameters_as_str())
        leaf = {'name': leaf_name}
        group.setdefault(ref.symmetry, []).append(leaf)
    data = {'name': '', 'children': []}
    for ac, v in uniprot_acs.iteritems():
        item = {'name': v, 'ac': ac, 'uniref': extra_info[v][0],
                'desc': extra_info[v][1]}
        if v in groups:
            item['children'] = [{'name': sg, 'children': pdbs}
                                for sg, pdbs in groups[v].items()]
        data['children'].append(item)
    with open(OUTPUT_JSON_FILE, 'w') as out:
        json.dump(data, out, indent=1, separators=(',', ': '))

def main(verbose=False):
    page = cached_urlopen(WIKI_URL, -1).readlines()
    parsed_page = parse_wiki_page(page) # { UniProt name: [some PDB IDs] }
    uniprot_acs = uniprot_names_to_acs(parsed_page) # { AC: UniProt name }
    pdb2up = read_sifts_mapping()
    descriptions_from_wiki = {}

    # check for missing PDB IDs that were given explicitely in the wiki
    for u, (desc, pp) in parsed_page.iteritems():
        descriptions_from_wiki[u] = desc
        for p in pp:
            if p not in pdb2up:
                print 'Missing in SIFTS: %s' % p
                pdb2up[p] = [k for k, v in uniprot_acs.iteritems() if v == u]

    # check for deprecated/removed PDB entries
    current_pdb_entries = read_current_pdb_entries()
    obsolete = [k for k in pdb2up if k not in current_pdb_entries]
    print 'Obsolete PDB entries in SIFTS:', len(obsolete)
    for k in obsolete:
        del pdb2up[k]

    # reverse the PDB->UniProt mapping
    ac2pdb = {}
    for p, acs in pdb2up.iteritems():
        if len(acs) == 1: # don't care about heteromers
            ac2pdb.setdefault(acs[0], []).append(p)

    uniref_clusters = fetch_uniref_clusters(uniprot_acs, verbose)

    # mapping back uniref id to uniprot name that was used an input
    uniref_sources = {}
    for uniref_name, uclust in uniref_clusters.iteritems():
        sources = [uniprot_acs[ac] for ac in uclust if ac in uniprot_acs]
        assert len(sources) == 1, sources
        uniref_sources[uniref_name] = sources[0]

    representants = []
    empty_cnt = 0
    for uniref_name, uclust in uniref_clusters.iteritems():
        pdb_set = sum([ac2pdb[ac] for ac in uclust if ac in ac2pdb], [])
        print '%s -> %s (%d entries) -> %d PDBs' % (
                uniref_sources[uniref_name], uniref_name,
                len(uclust), len(pdb_set))
        cells = []
        for pdb_id in pdb_set:
            cell = fetch_pdb_info_from_ebi(pdb_id)
            if cell is not None:
                cells.append(cell)

        # FIXME: should we exclude some mutants? Especially if many entries
        # are deposited for the same protein.
        # Maybe remove those with "Engineered mutation". We could use
        # http://www.ebi.ac.uk/pdbe/api/pdb/entry/mutated_AA_or_NA/:pdbid
        # Or Wild Type Search from RCSB.

        prev_len = len(representants)
        for pdb_cluster in get_pdb_clusters(cells):
            r = get_representative_unit_cell(pdb_cluster)
            r.uniprot_name = uniref_sources[uniref_name]
            r.comment = uniref_name
            if verbose:
                print ' '.join(p.pdb_id for p in pdb_cluster), '=>', r.pdb_id
                print('\t%s %-9s (%6.2f %6.2f %6.2f %5.1f %5.1f %5.1f) %8.2f'
                      % ((r.pdb_id, r.symmetry) + r.cell + (r.quality,)))
            representants.append(r)
        if len(representants) == prev_len:
            empty_cnt += 1
        if not verbose:
            print '\t-> %d' % (len(representants) - prev_len)
    representants.sort(key=lambda x: x.a)
    write_output_file(representants)
    print '%d (incl. %d absent) proteins -> %d unique unit cells -> %s' % (
            len(uniprot_acs), empty_cnt, len(representants), OUTPUT_PY_FILE)
    write_json_file(uniprot_acs,
                    dict((v, (k, descriptions_from_wiki[v]))
                         for k, v in uniref_sources.items()),
                    representants)
    a, b = min(itertools.combinations(representants, 2),
               key=lambda x: cell_distance(*x))
    print 'Min. dist: %.2f%% between %s and %s' % (100*cell_distance(a, b),
                                                   a.pdb_id, b.pdb_id)
    contaminer_acs = fetch_contaminer_protein_acs()
    for ac in contaminer_acs:
        if not any(ac in uclust for uclust in uniref_clusters.values()):
            print 'Only in ContaBase: %s %s' % (ac, contaminer_acs[ac])
    for uniref_name, uclust in uniref_clusters.iteritems():
        if len(set(uclust).intersection(contaminer_acs.viewkeys())) == 0:
            src = uniref_sources[uniref_name]
            if src not in NOT_IN_CONTABASE: # these are not news
                print 'Not in ContaBase:', src

if __name__ == '__main__':
    verbose = ('-v' in sys.argv[1:])
    main(verbose)
