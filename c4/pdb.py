import gzip
import os
import subprocess
import sys
import urllib2

from c4.cell import Cell
from c4.utils import comment, put_error


class PdbMeta(Cell):
    def __init__(self, cryst1_line):
        assert cryst1_line.startswith("CRYST1")
        a = float(cryst1_line[6:15])
        b = float(cryst1_line[15:24])
        c = float(cryst1_line[24:33])
        alpha = float(cryst1_line[33:40])
        beta = float(cryst1_line[40:47])
        gamma = float(cryst1_line[47:54])
        Cell.__init__(self, [a, b, c, alpha, beta, gamma])
        self.symmetry = cryst1_line[55:66].strip()
        try:
            self.z = int(cryst1_line[66:70])
        except ValueError:
            self.z = None

    def __str__(self):
        return '''\
cell: %(cell)s
symmetry: "%(symmetry)s"''' % self.__dict__


def read_metadata(pdb):
    if pdb.endswith('.gz'):
        f = gzip.open(pdb, 'rb')
    else:
        f = open(pdb)
    for line in f:
        if line.startswith("CRYST1"):
            f.close()
            return PdbMeta(line)
    f.close()
    sys.stderr.write("\nCRYST1 line not found in %s\n" % pdb)

def get_protein_mw(pdb):
    p = subprocess.Popen(["rwcontents", "XYZIN", pdb],
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = p.communicate(input="END\n")
    retcode = p.poll()
    if retcode:
        raise RuntimeError("rwcontents of %s failed." % pdb)
    key = 'Molecular Weight of protein:'
    start = stdoutdata.find(key)
    if start != -1:
        start += len(key)
        return float(stdoutdata[start:start+50].split()[0])

def remove_hetatm(filename_in, file_out):
    "remove HETATM and related lines"
    file_in = open(filename_in)
    removed = set()

    def is_removed(serial):
        return serial and not serial.isspace() and int(serial) in removed

    for line in file_in:
        record = line[:6]
        if record == "HETATM":
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
    if sys.argv[0] < 2:
        sys.stderr.write("No filenames.\n")
        sys.exit(1)
    if sys.argv[1] == "nohet":
        remove_hetatm(sys.argv[2], sys.stdout)
    elif sys.argv[1] == "mw":
        print get_protein_mw(sys.argv[2])
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
            print("File: %s" % arg)
            print read_metadata(arg)

if __name__ == '__main__':
    main()
