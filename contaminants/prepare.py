#!/usr/bin/env python2
"""
Script that generates data.py
"""

import csv
import os
import sys
import urllib2

if __name__ == '__main__' and __package__ is None:
    sys.path.insert(1,
          os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import dimple.cell

def download_page_from_dimple_wiki():
    url = ('https://raw.githubusercontent.com/wiki/'
           'ccp4/dimple/Crystallization-Contaminants.md')
    response = urllib2.urlopen(url)
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

def main():
    #page = open('Crystallization-Contaminants.md')
    page = download_page_from_dimple_wiki()
    pdb_ids = extract_pdb_ids(page)
    infos = [fetch_pdb_info(pdb_id) for pdb_id in pdb_ids]
    for info in infos:
        if 'cell' not in info:
            sys.exit('Wrong item: %s' % info)
    # TODO: standarize unit cell settings
    infos.sort(key=lambda x: x['cell'].a)
    with open('data.py', 'w') as data_py:
        data_py.write('\nDATA = [\n')
        for info in infos:
            cell = info['cell']
            data = ['"%s"' % info['structureId'],
                    '"%s"' % cell.symmetry
                   ] + ['%6.2f' % x for x in cell.cell]
            data_py.write('(' + ', '.join(data) + '),\n')
        data_py.write(']')

if __name__ == '__main__':
    main()
