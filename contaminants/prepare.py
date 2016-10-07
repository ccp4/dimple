#!/usr/bin/env python2
"""
Script that generates data.py
"""

import csv
import math
import os
import re
import sys
import urllib
import urllib2
from collections import OrderedDict

try:
    import scratch
except ImportError:
    scratch = None

if __name__ == '__main__' and __package__ is None:
    sys.path.insert(1,
          os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import dimple.cell

WIKI_URL = ('https://raw.githubusercontent.com/wiki/'
            'ccp4/dimple/Crystallization-Contaminants.md')

CLUSTER_CUTOFF = 0.03

def download_page_from_dimple_wiki():
    return open(WIKI_URL.split('/')[-1])
    response = urllib2.urlopen(WIKI_URL)
    return response

def extract_pdb_ids(page):
    pdb_ids = []
    for line in page:
        if ' - ' in line:
            desc, ids = line.rsplit(' - ', 1)
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

def fetch_pdb_info(pdb_id):
    par_names = ('lengthOfUnitCellLatticeA,lengthOfUnitCellLatticeB,'
                 'lengthOfUnitCellLatticeC,unitCellAngleAlpha,'
                 'unitCellAngleBeta,unitCellAngleGamma')
    url = ('http://www.rcsb.org/pdb/rest/customReport.csv?pdbids=%s'
           '&customReportColumns=structureId,spaceGroup,' + par_names +
           ',Z_PDB,releaseDate,revisionDate,'
           'experimentalTechnique,resolution,rAll,rFree'
           '&service=wsfile&format=csv')
    response = urllib2.urlopen(url % pdb_id)
    reader = csv.DictReader(response)
    d = reader.next()
    assert pdb_id == d['structureId']
    if d['experimentalTechnique'] == 'X-RAY DIFFRACTION':
        sg = d['spaceGroup']
        parameters = [float(d[x]) for x in par_names.split(',')]
        d['cell'] = dimple.cell.Cell(tuple(parameters), symmetry=sg)
        # TODO: standarize unit cell settings
    else:
        print 'No unit cell for %s (method: %s)' % (
                    pdb_id, d['experimentalTechnique'])
    return d

def dump_uniprot_tab(pdb_ids):
    query_url = ('http://www.uniprot.org/uniprot/'
                 '?query=database:(type:pdb%%20%s)&sort=score&format=tab')
    pdb_ids = [x.strip(',') for x in sys.stdin.read().split()]
    assert all(len(x) == 4 for x in pdb_ids)
    print >>sys.stderr, 'Read %d IDs' % len(pdb_ids)
    for pdb_id in pdb_ids:
        response = urllib2.urlopen(query_url % pdb_id)
        text = response.read()
        print pdb_id
        print text
        print

def uniprot_names_to_acs(names):
    return scratch.ACS #XXX
    query_url = ('http://www.uniprot.org/uniprot/'
                 '?query=%s&columns=id,entry%%20name&format=tab')
    acs = OrderedDict()
    for name in names:
        response = urllib2.urlopen(query_url % name)
        lines = response.readlines()
        assert len(lines) == 2
        ac, entry_name = lines[1].strip().split('\t')
        assert entry_name == name
        acs[ac] = name
    return acs

def fetch_uniref_clusters(acs):
    return scratch.CLUSTERS
    query_url = ('http://www.uniprot.org/uniref/'
                 '?query=%s&fil=identity:1.0&format=tab')
    clusters = OrderedDict()
    for ac in acs:
        #print query_url % ac
        response = urllib2.urlopen(query_url % ac)
        reader = csv.DictReader(response, delimiter='\t')
        for n, d in enumerate(reader):
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
    #f = urllib2.urlopen('http://www.uniprot.org/docs/pdbtosp.txt')
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

def get_pdb_clusters(infos):
    # simplistic clustering
    # TODO: use scipy instead
    # scipy.spatial.distance.pdist(lambda a, b: cell_distance(a, b))
    # scipy.cluster.hierarchy.linkage
    rinfos = infos[::-1]
    rinfos.sort(key=lambda x: (x['cell'].symmetry, x['cell'].a))
    while rinfos:
        clu = [rinfos.pop()]
        for n in range(len(rinfos)-1, -1, -1):
            dist = cell_distance(rinfos[n], clu[0])
            #print ':::', rinfos[n]['cell'], clu[0]['cell'], dist
            if dist < CLUSTER_CUTOFF:
                clu.append(rinfos.pop(n))
        yield clu

def cell_distance(a_info, b_info):
    a = a_info['cell']
    b = b_info['cell']
    if a.symmetry != b.symmetry:
        return float('inf')
    return a.max_shift_in_mapping(b)

def estimate_map_quality(info):
    # use simplistic criterium, similar to the one from
    # http://www.rcsb.org/pdb/statistics/clusterStatistics.do
    # TODO: for our purpose we regard those with unit-cell axes in non-standard
    # order as low quality
    # return -20
    try:
        resolution = float(info['resolution'])
        rfree = float(info['rFree'])
    except ValueError:
        return -10
    return 1.0 / resolution - rfree  # higher is better

# Representative entry is picked as a one with at least median "quality"
# metric that is least dissimilar to other cells.
def get_representative_unit_cell(infos):
    if len(infos) < 2:
        return infos[0]
    infos.sort(key=lambda x: -estimate_map_quality(x))
    cutoff = int(math.ceil(len(infos) / 2))
    def med_metric(inf):
        return max(cell_distance(inf, other) for other in infos)
    best = min(infos[:cutoff], key=med_metric)
    #print '===>', med_metric(best)
    return best

def write_data_py(representants):
    with open('data.py', 'w') as data_py:
        data_py.write('\nDATA = [\n')
        for info in representants:
            cell = info['cell']
            data = ['"%s"' % info['structureId'],
                    '"%s"' % cell.symmetry
                   ] + ['%6.2f' % x for x in cell.cell]
            data_py.write('(' + ', '.join(data) + '),\n')
        data_py.write(']')

def main():
    page = download_page_from_dimple_wiki()
    pdb_ids_from_page = extract_pdb_ids(page)
    pdbtosp = read_pdbtosp()
    missing = [p for p in pdb_ids_from_page if p not in pdbtosp]
    print 'missing:', missing

    uniprot_names = extract_uniprot_names(page)
    acs = uniprot_names_to_acs(uniprot_names)
    clusters = fetch_uniref_clusters(acs)
    homomers = group_homomeric_crystals_by_ac(pdbtosp)
    representants = []
    for name, clust in clusters.items():
        pdb_set = sum([homomers[c] for c in clust if c in homomers], [])
        print name, len(pdb_set)
        infos = [fetch_pdb_info(pdb_id) for pdb_id in pdb_set]
        for pdb_cluster in get_pdb_clusters(infos):
            r = get_representative_unit_cell(pdb_cluster)
            print r['structureId'], r['cell'].symmetry, r['cell']
            representants.append(r)
    representants.sort(key=lambda x: x['cell'].a)
    write_data_py(representants)



if __name__ == '__main__':
    main()
