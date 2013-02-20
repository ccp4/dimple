#!/usr/bin/env python

"""\
Usage: dimple input.mtz input.pdb output_dir
"""

import os
import sys
from c4.workflow import Workflow, JobError, put, put_error, \
                        show_info, open_pickled_workflow

RFREE_FOR_MOLREP = 0.4

def run_pipeline(wf, pdb, mtz):
    put("Change mtz symmetry if needed. Use pdb as reference.\n")
    wf.pointless(hklin=mtz, xyzin=pdb, hklout="pointless.mtz").run()
    pdb_meta = wf.read_pdb_metadata(pdb)
    mtz_meta = wf.read_mtz_metadata(mtz)
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
        refmac_rigid_xyzin = pdb
    else:
        put("Different cell in pdb %s ...\n" % str(pdb_meta.cell))
        put("              and mtz %s, changing pdb\n" % str(mtz_meta.cell))
        wf.change_pdb_cell(xyzin=pdb, xyzout="prepared.pdb",
                           cell=mtz_meta.cell)
        refmac_rigid_xyzin = "prepared.pdb"

    refmac_labin = "FP=F SIGFP=SIGF FREE=FreeR_flag"
    refmac_labout = ("FC=FC PHIC=PHIC FWT=2FOFCWT PHWT=PH2FOFCWT "
                     "DELFWT=FOFCWT PHDELWT=PHFOFCWT")
    put("Rigid-body refinement.\n")
    wf.refmac5(hklin="prepared.mtz", xyzin=refmac_rigid_xyzin,
               hklout="refmacRB.mtz", xyzout="refmacRB.pdb",
               labin=refmac_labin, labout=refmac_labout,
               keys="""refinement type rigidbody resolution 15 3.5
                       scale type simple lssc anisotropic experimental
                       solvent yes vdwprob 1.4 ionprob 0.8 mshrink 0.8
                       rigidbody ncycle 10""").run()

    if "free_r" not in wf.jobs[-1].data:
        put("WARNING: unknown free_r, something went wrong.\n")
        refmac_xyzin="refmacRB.pdb"
    elif wf.jobs[-1].data["free_r"] > RFREE_FOR_MOLREP:
        put("Run MR for R_free > %g\n" % RFREE_FOR_MOLREP)
        wf.molrep(f="prepared.mtz", m="refmacRB.pdb").run()
        refmac_xyzin="molrep.pdb"
        # we may add refmac (restr,ncyc=5) and findwaters here
    else:
        put("No MR for R_free < %g\n" % RFREE_FOR_MOLREP)
        refmac_xyzin="refmacRB.pdb"

    put("Final restrained refinement.\n")
    wf.refmac5(hklin="prepared.mtz", xyzin=refmac_xyzin,
               hklout="final.mtz", xyzout="final.pdb",
               labin=refmac_labin, labout=refmac_labout,
               keys="""make hydrogen all hout no cispeptide yes ssbridge yes
                       refinement type restrained
                       weight matrix 0.2
                       scale type simple lssc anisotropic experimental
                       solvent yes vdwprob 1.4 ionprob 0.8 mshrink 0.8
                       ncycle 8""").run()
    wf.find_blobs("final.mtz", "final.pdb", sigma=0.8).run()


def main():
    if "CCP4" not in os.environ:
        sys.stderr.write('$CCP4 not found, giving up\n')
        sys.exit(1)

    if len(sys.argv) >= 3 and sys.argv[1] == "info":
        wf = open_pickled_workflow(sys.argv[2])
        job_numbers = [int(job_str)-1 for job_str in sys.argv[3:]]
        show_info(wf, job_numbers)
        return

    try:
        assert len(sys.argv) == 4, "3 arguments expected"
        mtz, pdb, output_dir = sys.argv[1:]
        if mtz.endswith(".pdb") and pdb.endswith(".mtz"): # no problem
            mtz, pdb = pdb, mtz
        assert mtz.endswith(".mtz"), "1st arg should be mtz file"
        assert pdb.endswith(".pdb"), "2nd arg should be pdb file"
        if os.path.exists(output_dir):
            assert os.path.isdir(output_dir), "Not a directory: " + output_dir
        for filename in [mtz, pdb]:
            assert os.path.isfile(filename), "File not found: " + filename
    except AssertionError as e:
        put_error(e, __doc__.rstrip())
        sys.exit(1)
    wf = Workflow(output_dir)
    try:
        run_pipeline(wf=wf, mtz=os.path.abspath(mtz), pdb=os.path.abspath(pdb))
    except JobError as e:
        put_error(e.msg, comment=e.note)
        wf.pickle_jobs()
        sys.exit(1)
    wf.pickle_jobs()

if __name__ == "__main__":
    main()

