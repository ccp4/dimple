import gzip
import sys

from c4.cell import Cell


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


if __name__ == '__main__':
    if sys.argv[0] < 2:
        sys.stderr.write("No filenames.\n")
        sys.exit(1)
    if sys.argv[1] == "nohet":
        remove_hetatm(sys.argv[2], sys.stdout)
        sys.exit(0)
    for arg in sys.argv[1:]:
        print("File: %s" % arg)
        print read_metadata(arg)
