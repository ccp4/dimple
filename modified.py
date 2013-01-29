#!/usr/bin/env python

"""\
Usage: dimple input.mtz input.pdb output_dir
"""

import os
import sys
import c4.workflow
import c4.mtz
import c4.pdb
import c4.utils

def run_pipeline(wf, pdb, mtz):
    wf.pointless(hklin=mtz, xyzin=pdb, hklout="pointless.mtz").run()
    pdb_meta = c4.pdb.read_metadata(pdb)
    mtz_meta = c4.mtz.read_metadata(mtz)
    wf.unique(hklout="unique.mtz",
              cell=mtz_meta.cell, symmetry=pdb_meta.symmetry,
              resolution=mtz_meta.resolution_range[1],
              labout="F=F_UNIQUE SIGF=SIGF_UNIQUE").run()
    wf.freerflag(hklin="unique.mtz", hklout="free.mtz").run()
    if mtz_meta.symmetry == pdb_meta.symmetry:
        print "Same symmetry in pdb and mtz (%s), good." % pdb_meta.symmetry
        os.symlink("pointless.mtz", "spacegroup.mtz")
    else:
        print "Different symmetry in pdb (%s) and mtz (%s), reindexing mtz" % (
                pdb_meta.symmetry, mtz_meta.symmetry)
        wf.reindex(hklin="pointless.mtz", hklout="spacegroup.mtz",
                   symmetry=pdb_meta.symmetry).run()
    wf.truncate(hklin="spacegroup.mtz", hklout="truncate.mtz",
                labin="IMEAN=IMEAN SIGIMEAN=SIGIMEAN",
                labout="F=F SIGF=SIGF").run()
    wf.cad(hklin=["truncate.mtz", "free.mtz"], hklout="prepared.mtz",
           keys="""labin file 1 ALL
                   labin file 2 ALL""").run()
    if all(pdb_meta.cell[i] - mtz_meta.cell[i] < 1e-3 for i in range(6)):
        print "Cell dimensions in pdb and mtz are the same, good."
        os.symlink(pdb, "prepared.pdb")
    else:
        print "Different cell in pdb %s ..." % str(pdb_meta.cell)
        print "              and mtz %s, changing pdb" % str(mtz_meta.cell)
        wf.change_pdb_cell(xyzin=pdb, xyzout="prepared.pdb",
                           cell=mtz_meta.cell)
    refmac_labin = "FP=F SIGFP=SIGF FREE=FreeR_flag"
    refmac_labout = ("FC=FC PHIC=PHIC FWT=2FOFCWT PHWT=PH2FOFCWT "
                     "DELFWT=FOFCWT PHDELWT=PHFOFCWT")
    wf.refmac5(hklin="prepared.mtz", xyzin="prepared.pdb",
               hklout="refmacRB.mtz", xyzout="refmacRB.pdb",
               labin=refmac_labin, labout=refmac_labout,
               keys="""refinement type rigidbody resolution 15 3.5
                       scale type simple lssc anisotropic experimental
                       solvent yes vdwprob 1.4 ionprob 0.8 mshrink 0.8
                       rigidbody ncycle 10""").run_and_parse()
    wf.refmac5(hklin="prepared.mtz", xyzin="refmacRB.pdb",
               hklout="final.mtz", xyzout="final.pdb",
               labin=refmac_labin, labout=refmac_labout,
               keys="""make hydrogen all hout no cispeptide yes ssbridge yes
                       refinement type restrained
                       weight matrix 0.2
                       scale type simple lssc anisotropic experimental
                       solvent yes vdwprob 1.4 ionprob 0.8 mshrink 0.8
                       ncycle 8""").run_and_parse()


def main():
    if "CCP4" not in os.environ:
        sys.stderr.write('$CCP4 not found, giving up\n')
        sys.exit(1)
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
    except AssertionError, e:
        sys.stderr.write(c4.utils.red("Error: %s." % e) + "\n%s" % __doc__)
        sys.exit(1)
    wf = c4.workflow.Workflow(output_dir)
    run_pipeline(wf=wf, mtz=os.path.abspath(mtz), pdb=os.path.abspath(pdb))

if __name__ == "__main__":
    main()

