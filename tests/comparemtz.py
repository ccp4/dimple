#!/usr/bin/env cctbx.python

USAGE = "Usage: comparemtz.py [-r] file1.mtz..."

import os
import sys
import iotbx.mtz

# Columns in refmac HKLOUT:
#  - FC_ALL is maximum likelihood scaled and is used for map calculation
#  - FC_ALL_LS is least-squares scaled and is used Rfactor calculations
#  - F (or FP or anything else, depending on LABIN keyword) -- this is
#    input Fobs scaled by Refmac. Garib plans to change it to FP_refmac.
FO_COL = 'F'
FC_COL = 'FC_ALL_LS'
#FC_COL = 'FC_ALL'

def get_arr(miller_arrays, label):
    for ma in miller_arrays:
        if ma.info().labels[0] == label:
            return ma
    info = ', '.join(ma.info().labels[0] for ma in miller_arrays)
    sys.exit("ERROR: no column %s, only: %s" % (label, info))

def get_common_reflections(mtz_objects):
    common = get_arr(mtz_objects[0].as_miller_arrays(), FO_COL)
    for mtz_obj in mtz_objects[1:]:
        fs = get_arr(mtz_obj.as_miller_arrays(), FO_COL)
        common = common.common_set(fs)
    return common

def read_files(mtz_files):
    mtz_objects = []
    for n, path in enumerate(mtz_files):
        if not os.path.exists(path):
            print 'Not found: %s' % path
            continue
        print 'Reading [%d] %-30s' % (n+1, path),
        sys.stdout.flush()
        obj = iotbx.mtz.object(path)
        res = obj.max_min_resolution()
        print ' %s   %d refl.  %.1f - %.2f' % (obj.space_group_name(),
                obj.n_reflections(), res[0], res[1])
        mtz_objects.append(obj)
        obj.name = str(n+1)
        # special handling for DLS pipelines
        if path.endswith('/dimple/final.mtz'):
            obj.name = path.split('/')[-3]
            if obj.name.endswith('-run'):
                obj.name = obj.name[:-4]
    return mtz_objects

def main():
    use_cc = True
    if '-r' in sys.argv:
        use_cc = False
        sys.argv.remove('-r')
    if len(sys.argv) == 1:
        sys.exit(USAGE)
    if len(sys.argv) == 2 and not sys.argv[1].endswith('.mtz'):
        # special case to easily compare data in "processed" in DLS
        mtz_files = [os.path.join(sys.argv[1], ap, 'dimple/final.mtz')
                     for ap in ['fast_dp', 'autoPROC/ap-run', 'xia2/3dii-run',
                                'xia2/3d-run', 'xia2/dials-run']]
    else:
        mtz_files = sys.argv[1:]

    mtz_objects = read_files(mtz_files)
    if not mtz_objects:
        sys.exit("No mtz files")

    common_refl = get_common_reflections(mtz_objects)
    n_common = common_refl.size()
    cres = common_refl.d_max_min()

    print "Common res: %.2f - %.2f." % cres,
    print "Common reflections: %d" % n_common
    print """\
                 #extra refl            %s             free/work gap
input     res.  all same_res   all same_res same_refl  all same_res same_refl\
""" % (' CC  ' if use_cc else 'Rfree')
    for mtz_obj in mtz_objects:
        miller_arrays = mtz_obj.as_miller_arrays()
        fobs = get_arr(miller_arrays, FO_COL)
        assert fobs.is_in_asu()
        # we remove Fc's missing reflections (based on Fobs)
        fcalc = get_arr(miller_arrays, FC_COL).amplitudes().common_set(fobs)
        assert fcalc.is_in_asu()
        flags = get_arr(miller_arrays, 'FreeR_flag').common_set(fobs)
        n_same = fobs.resolution_filter(cres[0], cres[1]).size()
        #cc = fobs.correlation(fcalc).coefficient()
        def stat(fo, fc):
            if use_cc:
                return fo.correlation(fc).coefficient()
            else:
                return fo.r1_factor(fc)
        def calc(fo, fc):
            return (stat(fo, fc),
                    stat(fo.resolution_filter(cres[0], cres[1]),
                         fc.resolution_filter(cres[0], cres[1])),
                    stat(fo.common_set(common_refl),
                         fc.common_set(common_refl)))
        free_flag = (flags.data() == 0)
        free_data = calc(fobs.select(free_flag), fcalc.select(free_flag))
        work_data = calc(fobs.select(~free_flag), fcalc.select(~free_flag))
        if use_cc:
            gap_data = [(w-f)*100./f for f, w in zip(free_data, work_data)]
        else:
            gap_data = [(f-w)*100./f for f, w in zip(free_data, work_data)]
        print "%-9s %4.2f %6d %5d  " % (mtz_obj.name, fobs.d_min(),
                                        fobs.size()-n_common, n_same-n_common),
        print "%.4f  %.4f  %.4f  " % free_data,
        print "%.1f%%  %.1f%%  %.1f%%" % tuple(gap_data)

main()
