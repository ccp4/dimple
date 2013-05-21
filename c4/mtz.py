
import subprocess
from collections import OrderedDict

class MtzMeta:
    def __init__(self, cell, symmetry, sg_number, dmin, dmax, columns):
        assert type(columns[0]) == tuple
        self.cell = cell
        if cell:
            self.a, self.b, self.c = cell[0:3]
            self.alpha, self.beta, self.gamma = cell[3:]
        self.symmetry = symmetry
        self.sg_number = sg_number
        self.dmin = dmin
        self.dmax = dmax
        assert dmin >= dmax # yes, min > max here
        self.columns = OrderedDict(columns)
    def __str__(self):
        return """\
cell: %(cell)s
symmetry: "%(symmetry)s" (symmetry group no. %(sg_number)d)
resolution range: %(dmin)s - %(dmax)s
columns: %(columns)s""" % self.__dict__


def read_metadata(hklin):
    "for now using mtzdump, directly calling libccp4/mtzlib would be better"
    p = subprocess.Popen(["mtzdump", "HKLIN", hklin],
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    stdoutdata, stderrdata = p.communicate(input="HEAD\nEND")
    retcode = p.poll()
    if retcode:
        raise RuntimeError("mtzdump of %s failed." % hklin)
    lines = stdoutdata.splitlines()
    for n, line in enumerate(lines):
        if not line.startswith(" * "):
            continue
        if line.startswith(" * Dataset ID, project/crystal/dataset names, ce"):
            try:
                cell = tuple(float(x) for x in lines[n + 5].split())
            except ValueError:
                cell = None
        elif line.startswith(" * Space group = "):
            symmetry = line.split("'")[1].strip()
            sg_number = int(line.split()[-1].rstrip(")"))
        elif line.startswith(" *  Resolution Range"):
            res_line = lines[n+2]
            #    0.00125    0.48742     (     28.339 -      1.432 A )
            lower_resol = float(res_line[28:40])
            upper_resol = float(res_line[41:53])
        elif line.startswith(" * Column Labels"):
            column_names = lines[n+2].split()
        elif line.startswith(" * Column Types"):
            column_types = lines[n+2].split()
    columns = zip(column_names, column_types)
    return MtzMeta(cell, symmetry=symmetry, sg_number=sg_number,
                   dmin=lower_resol, dmax=upper_resol, columns=columns)


def check_freerflags_column(free_mtz, data_mtz_meta):
    rfree_meta = read_metadata(free_mtz)
    if rfree_meta.dmax > data_mtz_meta.dmax:
        raise ValueError("free-R-flags dmax: %g (should be < %g)" %
                                  (rfree_meta.dmax, data_mtz_meta.dmax))
    if rfree_meta.dmin < data_mtz_meta.dmin:
        raise ValueError("free-R-flags dmin: %g (should be >= %g)" %
                                  (rfree_meta.dmin, data_mtz_meta.dmin))
    for col_label, col_type in rfree_meta.columns.iteritems():
        if col_label.lower().startswith('free') and col_type == 'I':
            return col_label
    raise ValueError("free-R column not found in %s" % free_mtz)


if __name__ == '__main__':
    import sys
    if sys.argv[0] < 2:
        sys.stderr.write("No filenames.\n")
        sys.exit(1)
    for arg in sys.argv[1:]:
        print("File: %s" % arg)
        print read_metadata(arg)
