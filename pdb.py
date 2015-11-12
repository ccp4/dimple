import gzip
import os
import sys
import urllib2

from dimple.cell import Cell
from dimple.utils import comment, put_error


class PdbMeta(Cell):
    def __init__(self, cryst1_line):
        assert cryst1_line.startswith("CRYST1")
        a = float(cryst1_line[6:15])
        b = float(cryst1_line[15:24])
        c = float(cryst1_line[24:33])
        alpha = float(cryst1_line[33:40])
        beta = float(cryst1_line[40:47])
        gamma = float(cryst1_line[47:54])
        symmetry = cryst1_line[55:66].strip()
        Cell.__init__(self, (a, b, c, alpha, beta, gamma), symmetry)
        try:
            self.z = int(cryst1_line[66:70])
        except ValueError:
            self.z = None
        self.has_hetatm_x = None

    def __str__(self):
        return '''\
cell: %(cell)s
symmetry: "%(symmetry)s"''' % self.__dict__


def read_metadata(pdb):
    if pdb.endswith('.gz'):
        f = gzip.open(pdb, 'rb')
    else:
        f = open(pdb)
    meta = None
    for line in f:
        if line.startswith("CRYST1"):
            meta = PdbMeta(line)
            break
    if meta is None:
        if f.tell() == 0:
            put_error("empty file: %s" % pdb)
        else:
            put_error("CRYST1 line not found in %s" % pdb)
    f.close()
    return meta


def check_hetatm_x(filename, meta):
    if meta and meta.has_hetatm_x is not None:
        return meta.has_hetatm_x
    with open(filename) as f:
        has_hetatm_x = any(line[:6] == "HETATM" and line[76:78] == " X"
                           for line in f)
    if meta:
        meta.has_hetatm_x = has_hetatm_x
    return has_hetatm_x


def remove_hetatm(filename_in, file_out, remove_all):
    """Remove HETATM and related lines.
    If remove_all is False, remove only element X which happens in many PDB
    entries but is not accepted by pointless, refmac, phaser, etc
    """
    # we could instead zero occupancy of the atoms and replace X with Y
    file_in = open(filename_in)
    removed = set()

    def is_removed(serial):
        return serial and not serial.isspace() and int(serial) in removed

    for line in file_in:
        record = line[:6]
        if record == "HETATM":
            if remove_all or line[76:78] == " X":
                atom_serial_num = int(line[6:11])
                removed.add(atom_serial_num)
                continue
        elif record in ("HET   ", "HETNAM", "HETSYN", "FORMUL"):
            continue
        elif line.startswith("ANISOU"):
            if is_removed(line[6:11]):
                continue
        elif line.startswith("CONECT"):
            if any(is_removed(line[p:p+5]) for p in (6, 11, 16, 21, 26)):
                continue
        file_out.write(line)
    return len(removed)

def is_pdb_id(a):
    return len(a) == 4 and a[0].isdigit() and a[1:].isalnum()

def download_pdb(pdb_id, output_dir):
    filename = pdb_id.upper()+'.pdb'
    path = os.path.join(output_dir, filename)
    if os.path.exists(path):
        comment('%s: using existing file %s\n' % (pdb_id, filename))
    else:
        comment('Downloading %s from RCSB...  ' % pdb_id)
        url = 'http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=' + pdb_id.upper()
        try:
            u = urllib2.urlopen(url)
        except urllib2.HTTPError as e:
            put_error(str(e))
            sys.exit(1)
        content = u.read()
        with open(path, 'wb') as f:
            f.write(content)
        comment('done.\n')
    return path


def main():
    if len(sys.argv) < 2:
        sys.stderr.write("Usage: pdb.py [get|vol|nohet] file1.pdb ...\n")
        sys.exit(1)
    if sys.argv[1] == "nohet":
        remove_hetatm(sys.argv[2], sys.stdout, remove_all=True)
    elif sys.argv[1] == "vol":
        print read_metadata(sys.argv[2]).get_volume()
    elif sys.argv[1] == "get":
        for arg in sys.argv[2:]:
            if is_pdb_id(arg):
                path = download_pdb(arg, os.getcwd())
                print '-> ' + path
            else:
                sys.stderr.write('Error: %s is not a pdb code.\n' % arg)
                sys.exit(1)
    else:
        for arg in sys.argv[1:]:
            print "File: %s" % arg
            print read_metadata(arg)

if __name__ == '__main__':
    main()
