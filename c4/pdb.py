
class PdbMeta:
    def __init__(self, cryst1_line):
        assert cryst1_line.startswith("CRYST1")
        self.a = float(cryst1_line[6:15])
        self.b = float(cryst1_line[15:24])
        self.c = float(cryst1_line[24:33])
        self.alpha = float(cryst1_line[33:40])
        self.beta = float(cryst1_line[40:47])
        self.gamma = float(cryst1_line[47:54])
        self.cell = (self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
        self.symmetry = cryst1_line[55:66].strip()
        self.z = int(cryst1_line[66:70])
    def __str__(self):
        return """\
cell: %(cell)s
symmetry: "%(symmetry)s"
z: %(z)d""" % self.__dict__


def read_metadata(pdb):
    with open(pdb) as f:
        for line in f:
            if line.startswith("CRYST1"):
                return PdbMeta(line)


if __name__ == '__main__':
    import sys
    if sys.argv[0] < 2:
        sys.stderr.write("No filenames.\n")
        sys.exit(1)
    for arg in sys.argv[1:]:
        print("File: %s" % arg)
        print read_metadata(arg)
