#!/usr/bin/env cctbx.python
"""
Checks if two mtz files have the same free set for the common subset
of reflections. Compares only FreeR_flag (or FREE) columns.

In Dimple we generate the same free set when we can. This script was used
to verify it.
"""

import sys
import iotbx.mtz  # cctbx is handy for this

def read_integer_column(mtz_file, label):
    group = mtz_file.extract_integers(label)
    assert len(group.data) == len(group.indices)
    #print "#free refl. %4d of %5d" % (group.data.count(0), len(group.data))
    return dict(zip(group.indices, group.data))

def get_free_flags(filename):
    mtz_file = iotbx.mtz.object(filename)
    for label in ['FreeR_flag', 'FREE']:
        if mtz_file.has_column(label):
            return read_integer_column(mtz_file, label)
    sys.exit("Free set not found in %s" % filename)

def compare(f1, f2):
    common_hkl = set(f1) & set(f2)
    differ = sum(f1[hkl] != f2[hkl] for hkl in common_hkl)
    if differ == 0:
        print "SAME flags (%d)" % len(common_hkl)
    else:
        print "NOT same flags: %d of %d" % (differ, len(common_hkl))

def main():
    if len(sys.argv) < 3:
        sys.exit("Usage: check_same_free.py file1 file2...")
    f1 = get_free_flags(sys.argv[1])
    for arg in sys.argv[2:]:
        f2 = get_free_flags(arg)
        compare(f1, f2)

if __name__ == '__main__':
    main()

