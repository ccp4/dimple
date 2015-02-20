
import subprocess
from collections import OrderedDict
import sys
from c4.utils import put_error
from c4.cell import Cell

class MtzMeta(Cell):
    def __init__(self, cell, symmetry, sg_number, dmin, dmax, columns,
                 filename):
        assert type(columns[0]) == tuple
        Cell.__init__(self, cell)
        self.symmetry = symmetry
        self.sg_number = sg_number
        self.dmin = dmin
        self.dmax = dmax
        assert dmin >= dmax # yes, min > max here
        self.columns = OrderedDict(columns)
        self.filename = filename

    def __str__(self):
        return """\
cell: %(cell)s
symmetry: "%(symmetry)s" (symmetry group no. %(sg_number)d)
resolution range: %(dmin)s - %(dmax)s
columns: %(columns)s""" % self.__dict__

    def check_col_type(self, label, expected_type):
        if label not in self.columns:
            put_error("Column '%s' not found in %s" % (label, self.filename))
            sys.exit(1)
        col_type = self.columns[label]
        if col_type != expected_type:
            put_error("Column '%s' in %s has type '%s' (expected '%s')" %
                      (label, self.filename, col_type, expected_type))
            return False
        return True


def _run_mtzdump(hklin, keys):
    p = subprocess.Popen(["mtzdump", "HKLIN", hklin],
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    keys.append("END")
    stdoutdata, stderrdata = p.communicate(input="\n".join(keys))
    retcode = p.poll()
    if retcode:
        raise RuntimeError("mtzdump of %s failed." % hklin)
    return stdoutdata

def read_metadata(hklin):
    "for now using mtzdump, directly calling libccp4/mtzlib would be better"
    lines = _run_mtzdump(hklin, ["HEAD"]).splitlines()
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
                   dmin=lower_resol, dmax=upper_resol, columns=columns,
                   filename=hklin)

def check_freerflags_column(free_mtz, data_mtz_meta):
    names = ['FreeR_flag', 'FREE']
    rfree_meta = read_metadata(free_mtz)
    if rfree_meta.dmax > data_mtz_meta.dmax:
        raise ValueError("free-R-flags dmax: %g (should be < %g)" %
                                  (rfree_meta.dmax, data_mtz_meta.dmax))
    if rfree_meta.dmin < data_mtz_meta.dmin:
        raise ValueError("free-R-flags dmin: %g (should be >= %g)" %
                                  (rfree_meta.dmin, data_mtz_meta.dmin))
    for name in names:
        if name in rfree_meta.columns:
            rfree_meta.check_col_type(name, 'I')
            return name
    raise ValueError("free-R column not found in %s" % free_mtz)

def get_num_missing(hklin, col):
    # for now using mtzdump
    out = _run_mtzdump(hklin, ["NREF 0"])
    try:
        start = out.index('\n Col Sort')
        end = out.index('\n No. of reflections')
        lines = out[start+1:end].splitlines()[3:-1]
        for line in lines:
            sp = line.split()
            if sp[-1] == col:
                return int(sp[4])
    except ValueError:
        pass

# for testing only
def main():
    if len(sys.argv) < 2:
        sys.stderr.write("No filenames.\n")
        sys.exit(1)
    if sys.argv[1] == '-m' and len(sys.argv) >= 3:
        col = sys.argv[2]
        for arg in sys.argv[3:]:
            print arg, get_num_missing(arg, col)
    else:
        for arg in sys.argv[1:]:
            print("File: %s" % arg)
            print read_metadata(arg)

if __name__ == '__main__':
    main()
