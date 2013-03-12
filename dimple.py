#!/usr/bin/env python

"""\
Usage: dimple input.mtz input.pdb output_dir
"""

import os
import sys
from c4.workflow import Workflow, JobError, put, put_error, \
                        parse_workflow_commands

RFREE_FOR_MOLREP = 0.4

class Options:
    def __init__(self):
        self.pdb = None
        self.mtz = None
        self.final_pdb = "final.pdb"
        self.final_mtz = "final.mtz"

def run_pipeline(wf, opt):
    put("Change mtz symmetry if needed. Use pdb as reference.\n")
    wf.pointless(hklin=opt.mtz, xyzin=opt.pdb, hklout="pointless.mtz").run()
    pdb_meta = wf.read_pdb_metadata(opt.pdb)
    mtz_meta = wf.read_mtz_metadata(opt.mtz)
    if mtz_meta.symmetry == pdb_meta.symmetry:
        put(" Same symmetry in pdb and mtz (%s).\n" % pdb_meta.symmetry)
        truncate_hklin = "pointless.mtz"
    else:
        put(" Different symmetry in pdb (%s) and mtz (%s), reindexing mtz\n" %
                (pdb_meta.symmetry, mtz_meta.symmetry))
        wf.reindex(hklin="pointless.mtz", hklout="spacegroup.mtz",
                   symmetry=pdb_meta.symmetry).run()
        truncate_hklin = "spacegroup.mtz"

    put("Obtain structure factor amplitudes\n")
    wf.truncate(hklin=truncate_hklin, hklout="truncate.mtz",
                labin="IMEAN=IMEAN SIGIMEAN=SIGIMEAN",
                labout="F=F SIGF=SIGF").run()

    put("Add missing reflections and flag for cross-validation\n")
    wf.unique(hklout="unique.mtz",
              cell=mtz_meta.cell, symmetry=pdb_meta.symmetry,
              resolution=mtz_meta.resolution_range[1],
              labout="F=F_UNIQUE SIGF=SIGF_UNIQUE").run()
    wf.freerflag(hklin="unique.mtz", hklout="free.mtz").run()
    wf.cad(hklin=["truncate.mtz", "free.mtz"], hklout="prepared.mtz",
           keys="""labin file 1 ALL
                   labin file 2 E1=FreeR_flag
                   """).run()

    if all(pdb_meta.cell[i] - mtz_meta.cell[i] < 1e-3 for i in range(6)):
        put("Cell dimensions in pdb and mtz are the same.\n")
        correct_cell_pdb = opt.pdb
    else:
        put("Different cell in pdb %s ...\n" % str(pdb_meta.cell))
        put("              and mtz %s, changing pdb\n" % str(mtz_meta.cell))
        wf.change_pdb_cell(xyzin=opt.pdb, xyzout="prepared.pdb",
                           cell=mtz_meta.cell)
        correct_cell_pdb = "prepared.pdb"

    if False:
        rb_xyzin = "prepared_nohet.pdb"
        n_het = wf.remove_hetatm(xyzin=correct_cell_pdb, xyzout=rb_xyzin)
        put("Removed %s atoms marked as HETATM in pdb.\n" % n_het)
    else:
        rb_xyzin = correct_cell_pdb

    refmac_labin = "FP=F SIGFP=SIGF FREE=FreeR_flag"
    refmac_labout = ("FC=FC PHIC=PHIC FWT=2FOFCWT PHWT=PH2FOFCWT "
                     "DELFWT=FOFCWT PHDELWT=PHFOFCWT")
    put("Rigid-body refinement.\n")
    wf.refmac5(hklin="prepared.mtz", xyzin=rb_xyzin,
               hklout="refmacRB.mtz", xyzout="refmacRB.pdb",
               labin=refmac_labin, labout=refmac_labout,
               keys="""refinement type rigidbody resolution 15 3.5
                       scale type simple lssc anisotropic experimental
                       solvent yes vdwprob 1.4 ionprob 0.8 mshrink 0.8
                       rigidbody ncycle 10""").run()

    if "free_r" not in wf.jobs[-1].data:
        put("WARNING: unknown free_r, something went wrong.\n")
        refmac_xyzin = "refmacRB.pdb"
    elif wf.jobs[-1].data["free_r"] > RFREE_FOR_MOLREP:
        put("Run MR for R_free > %g\n" % RFREE_FOR_MOLREP)
        wf.molrep(f="prepared.mtz", m="refmacRB.pdb").run()
        refmac_xyzin = "molrep.pdb"
        # we may add refmac (restr,ncyc=5) and findwaters here
    else:
        put("No MR for R_free < %g\n" % RFREE_FOR_MOLREP)
        refmac_xyzin = "refmacRB.pdb"

    put("Final restrained refinement.\n")
    wf.refmac5(hklin="prepared.mtz", xyzin=refmac_xyzin,
               hklout=opt.final_mtz, xyzout=opt.final_pdb,
               labin=refmac_labin, labout=refmac_labout,
               keys="""make hydrogen all hout no cispeptide yes ssbridge yes
                       refinement type restrained
                       weight matrix 0.2
                       scale type simple lssc anisotropic experimental
                       solvent yes vdwprob 1.4 ionprob 0.8 mshrink 0.8
                       ncycle 8""").run()
    fb_job = wf.find_blobs(opt.final_mtz, opt.final_pdb, sigma=0.8).run()
    blobs = fb_job.data["blobs"]
    with open("coot.py", "wb") as f:
        f.write(coot_py_template % (opt.final_pdb, opt.final_mtz))
        if blobs:
            f.write("set_rotation_centre(%g,%g,%g)" % blobs[0])


def parse_dimple_commands(args):
    opt = Options()
    # special mode for compatibility with ccp4i
    if len(args) == 8 and args[::2] == ["HKLIN", "XYZIN", "HKLOUT", "XYZOUT"]:
        mtz, pdb, opt.final_mtz, opt.final_pdb = args[1::2]
        output_dir = os.path.join(os.environ["CCP4_SCR"], "dimple_out")
    else:
        assert len(args) == 3, "3 arguments expected"
        mtz, pdb, output_dir = args
        if mtz.endswith(".pdb") and pdb.endswith(".mtz"):  # no problem
            mtz, pdb = pdb, mtz
        assert mtz.endswith(".mtz"), "1st arg should be mtz file"
        assert pdb.endswith(".pdb"), "2nd arg should be pdb file"
        if os.path.exists(output_dir):
            assert os.path.isdir(output_dir), "Not a directory: " + output_dir
        for filename in [mtz, pdb]:
            assert os.path.isfile(filename), "File not found: " + filename
    opt.mtz = os.path.abspath(mtz)
    opt.pdb = os.path.abspath(pdb)
    return opt, output_dir

coot_py_template = """
pdb = "%s"
mtz = "%s"
center = (34.28,  31.86,  57.20)

set_nomenclature_errors_on_read("ignore")
molecule = read_pdb(pdb)
map21 = make_and_draw_map(mtz, "2FOFCWT", "PH2FOFCWT", "", 0, 0)
map11 = make_and_draw_map(mtz, "FOFCWT", "PHFOFCWT", "", 0, 1)
"""

def main():
    args = sys.argv[1:]
    if parse_workflow_commands(args):
        return

    for necessary_var in ("CCP4", "CCP4_SCR"):
        if necessary_var not in os.environ:
            put_error('$%s not found, giving up' % necessary_var)
            sys.exit(1)
    if not os.path.isdir(os.environ["CCP4_SCR"]):
        put_error('No such directory: $CCP4_SCR, refmac shall not work!')

    try:
        opt, output_dir = parse_dimple_commands(args)
    except AssertionError as e:
        put_error(e, __doc__.rstrip())
        sys.exit(1)

    wf = Workflow(output_dir)
    try:
        run_pipeline(wf=wf, opt=opt)
    except JobError as e:
        put_error(e.msg, comment=e.note)
        wf.pickle_jobs()
        sys.exit(1)
    wf.pickle_jobs()

if __name__ == "__main__":
    main()
