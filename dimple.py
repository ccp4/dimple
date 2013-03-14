#!/usr/bin/env python

import os
import sys
import argparse
from c4.workflow import Workflow, JobError, put, put_error, \
                        parse_workflow_commands

RFREE_FOR_MOLREP = 0.4

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
               hklout=opt.hklout, xyzout=opt.xyzout,
               labin=refmac_labin, labout=refmac_labout,
               keys="""make hydrogen all hout no cispeptide yes ssbridge yes
                       refinement type restrained
                       weight matrix 0.2
                       scale type simple lssc anisotropic experimental
                       solvent yes vdwprob 1.4 ionprob 0.8 mshrink 0.8
                       ncycle 8""").run()
    fb_job = wf.find_blobs(opt.hklout, opt.xyzout, sigma=0.8).run()
    blobs = fb_job.data["blobs"]
    wf.write_coot_script("coot.py", pdb=opt.xyzout, mtz=opt.hklout,
                         center=(blobs[0] if blobs else None))
    # For now not more than two blobs, in future better blob/ligand scoring
    for n, b in enumerate(blobs[:2]):
        wf.make_png("blob%s" % (n+1), pdb=opt.xyzout, mtz=opt.hklout, center=b)



def parse_dimple_commands():
    parser = argparse.ArgumentParser(
                              usage="dimple input.mtz input.pdb output_dir")
    parser.add_argument('mtz', metavar='input.mtz')
    parser.add_argument('pdb', metavar='input.pdb')
    parser.add_argument('output_dir')
    parser.add_argument('--from-job', metavar='N', type=int, default=0)
    parser.add_argument('--hklout', metavar='out.mtz', default='final.mtz')
    parser.add_argument('--xyzout', metavar='out.pdb', default='final.pdb')
    # get rid of 'positional arguments' in the usage method
    parser._action_groups[:1] = []

    args = sys.argv[1:]

    # special mode for compatibility with ccp4i
    legacy_args = { "HKLIN": "", "XYZIN": "",
                    "HKLOUT": "--hklout", "XYZOUT": "--xyzout" }
    if len(args) == 8 and args[0] in legacy_args:
        args = [legacy_args.get(a) or a
                for a in args if legacy_args.get(a) != ""]
        output_dir = os.path.join(os.environ["CCP4_SCR"], "dimple_out")
        args.append(output_dir)

    opt = parser.parse_args(args)
    if opt.mtz.endswith(".pdb") and opt.pdb.endswith(".mtz"):  # no problem
        opt.mtz, opt.pdb = opt.pdb, opt.mtz

    # extra checks
    if not opt.mtz.endswith(".mtz"):
        put_error("1st arg should be mtz file")
        sys.exit(1)
    if not opt.pdb.endswith(".pdb"):
        put_error("2nd arg should be pdb file")
        sys.exit(1)
    for filename in [opt.mtz, opt.pdb]:
        if not os.path.isfile(filename):
            put_error("File not found: " + filename)
            sys.exit(1)
    if os.path.exists(opt.output_dir) and not os.path.isdir(opt.output_dir):
        put_error("Not a directory: " + opt.output_dir)
        sys.exit(1)

    opt.mtz = os.path.relpath(opt.mtz, opt.output_dir)
    opt.pdb = os.path.relpath(opt.pdb, opt.output_dir)

    return opt


def main():
    if parse_workflow_commands():
        return

    for necessary_var in ("CCP4", "CCP4_SCR"):
        if necessary_var not in os.environ:
            put_error('$%s not found, giving up' % necessary_var)
            sys.exit(1)
    if not os.path.isdir(os.environ["CCP4_SCR"]):
        put_error('No such directory: $CCP4_SCR, refmac shall not work!')

    options = parse_dimple_commands()

    wf = Workflow(options.output_dir)
    wf.from_job = options.from_job
    try:
        run_pipeline(wf=wf, opt=options)
    except JobError as e:
        put_error(e.msg, comment=e.note)
        wf.pickle_jobs()
        sys.exit(1)
    wf.pickle_jobs()

if __name__ == "__main__":
    main()
