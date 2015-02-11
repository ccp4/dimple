#!/usr/bin/env python

import os
import sys
if sys.version_info[:2] != (2, 7):
    sys.stderr.write("Error. Python 2.7 is required.\n")
    sys.exit(1)
import argparse
from c4.utils import comment, put_error, syspath, adjust_path, start_log
from c4.mtz import check_freerflags_column
import c4.workflow
from c4 import coot

__version__ = '1.4'

def dimple(wf, opt):
    mtz_meta = wf.read_mtz_metadata(opt.mtz)
    wf.file_info[opt.mtz] = mtz_meta
    mtz_meta.check_col_type(opt.icolumn, 'J')
    mtz_meta.check_col_type(opt.sigicolumn, 'Q')
    _comment_summary_line("MTZ", mtz_meta)
    for p in opt.pdbs:
        wf.file_info[p] = wf.read_pdb_metadata(p)
    if len(opt.pdbs) > 1:
        comment("PDBs in order of similarity (using the first one):\n")
        opt.pdbs.sort(key=lambda x: calculate_difference_metric(wf.file_info[x],
                                                                mtz_meta))
    for p in opt.pdbs:
        _comment_summary_line(os.path.basename(p), wf.file_info[p])
    ini_pdb = opt.pdbs[0]
    pdb_meta = wf.file_info[ini_pdb]

    wf.pointless(hklin=opt.mtz, xyzin=ini_pdb, hklout="pointless.mtz",
                 keys="TOLERANCE 5").run()
    alt_reindex = wf.jobs[-1].data.get('alt_reindex')
    if alt_reindex:
        for ar in alt_reindex:
            comment("    %-10s CC: %-8s cell_deviation: %s\n" % (
                    ar['op'], ar['cc'], ar['cell_deviat']))
    else:
        comment("    ooops, no good indexing\n")
    #comment("Calculate structure factor amplitudes\n")
    wf.truncate(hklin="pointless.mtz", hklout="truncate.mtz",
                labin="IMEAN=%s SIGIMEAN=%s" % (opt.icolumn, opt.sigicolumn),
                labout="F=F SIGF=SIGF").run()

    if opt.free_r_flags:
        free_mtz = opt.free_r_flags
        try:
            free_col = check_freerflags_column(free_mtz, mtz_meta)
        except ValueError, e: # avoiding "as e" syntax for the sake of Py2.4
            put_error(e)
            sys.exit(1)
        comment("Free R flags from given file, col. %s.\n" % free_col)
    else:
        comment("Add missing reflections and free-R flags\n")
        free_mtz = "free.mtz"
        wf.unique(hklout="unique.mtz",
                  cell=mtz_meta.cell, symmetry=pdb_meta.symmetry,
                  resolution=mtz_meta.dmax,
                  labout="F=F_UNIQUE SIGF=SIGF_UNIQUE").run()
        wf.freerflag(hklin="unique.mtz", hklout=free_mtz).run()
        free_col = 'FreeR_flag'

    wf.cad(hklin=["truncate.mtz", free_mtz], hklout="prepared.mtz",
           keys="""labin file 1 ALL
                   labin file 2 E1=%s
                   reso file 2 1000.0 %g
                   """ % (free_col, mtz_meta.dmax)).run()
    if False:
        rb_xyzin = "prepared_nohet.pdb"
        n_het = wf.remove_hetatm(xyzin=ini_pdb, xyzout=rb_xyzin)
        comment("Removed %s atoms marked as HETATM in pdb.\n" % n_het)
    else:
        rb_xyzin = ini_pdb

    refmac_labin = "FP=F SIGFP=SIGF FREE=%s" % free_col
    refmac_labout = ("FC=FC PHIC=PHIC FWT=2FOFCWT PHWT=PH2FOFCWT "
                     "DELFWT=FOFCWT PHDELWT=PHFOFCWT")
    comment("Rigid-body refinement with resolution 3.5 A, 10 cycles.\n")
    wf.refmac5(hklin="prepared.mtz", xyzin=rb_xyzin,
               hklout="refmacRB.mtz", xyzout="refmacRB.pdb",
               labin=refmac_labin, labout=refmac_labout, libin=None,
               keys="""refinement type rigidbody resolution 15 3.5
                       scale type simple lssc anisotropic experimental
                       solvent yes vdwprob 1.4 ionprob 0.8 mshrink 0.8
                       rigidbody ncycle 10""").run()

    if "free_r" not in wf.jobs[-1].data:
        comment("WARNING: unknown free_r, something went wrong.\n")
        refmac_xyzin = "refmacRB.pdb"
    elif wf.jobs[-1].data["free_r"] > opt.mr_when_rfree:
        comment("Run MR for R_free > %g\n" % opt.mr_when_rfree)
        wf.molrep(f="prepared.mtz", m="refmacRB.pdb").run()
        refmac_xyzin = "molrep.pdb"
    else:
        comment("No MR for R_free < %g\n" % opt.mr_when_rfree)
        refmac_xyzin = "refmacRB.pdb"

    if False:
        wf.findwaters(pdbin=refmac_xyzin, hklin="refmacRB.mtz",
                      f="FC", phi="PHIC", pdbout="prepared_wat.pdb", sigma=2)
        refmac_xyzin = "prepared_wat.pdb"

    comment("Final restrained refinement, 8 cycles.\n")
    if opt.weight:
        refmac_weight = "matrix 0.2"
    else:
        refmac_weight = "auto"
    restr_job = wf.refmac5(hklin="prepared.mtz", xyzin=refmac_xyzin,
                 hklout=opt.hklout, xyzout=opt.xyzout,
                 labin=refmac_labin, labout=refmac_labout, libin=opt.libin,
                 keys="""make hydrogen all hout no cispeptide yes ssbridge yes
                         refinement type restrained
                         weight %s
                         scale type simple lssc anisotropic experimental
                         solvent yes vdwprob 1.4 ionprob 0.8 mshrink 0.8
                         ncycle 8""" % refmac_weight).run()
    if opt.summary:
        comment("".join(restr_job.data["selected_lines"]))
    # if that run is repeated with --from-job it's useful to compare Rfree
    if wf.from_job > 0 and wf.from_job <= len(wf.jobs): # from_job is 1-based
        prev = [j for j in wf.repl_jobs if j.name == restr_job.name]
        if prev and prev[0].data and "free_r" in prev[0].data:
            comment("Previously:  R/Rfree %.4f/%.4f  Rfree change: %+.4f\n" % (
                    prev[0].data["overall_r"], prev[0].data["free_r"],
                    restr_job.data["free_r"] - prev[0].data["free_r"]))

    fb_job = wf.find_blobs(opt.hklout, opt.xyzout, sigma=0.8).run()
    if opt.img_format == 'none':
        return
    blobs = fb_job.data["blobs"]
    if blobs:
        if len(blobs) == 1:
            comment("Rendering density blob at (%.1f, %.1f, %.1f)\n" % blobs[0])
        else:
            comment("Rendering 2 largest blobs: at (%.1f, %.1f, %.1f) and at "
                    "(%.1f, %.1f, %.1f)\n" % (blobs[0]+blobs[1]))
        if _check_picture_tools():
            _generate_pictures(wf, opt, fb_job)
    else:
        comment("Unmodelled blobs not found.\n")


def _comment_summary_line(name, meta):
    def angle(x):
        if x == 90.: return '90'
        else:        return str(x)
    if meta:
        line = '%-21s %-12s (%.2f, %.2f, %.2f,  %s, %s, %s)\n' % (
                name, meta.symmetry, meta.a, meta.b, meta.c,
                angle(meta.alpha), angle(meta.beta), angle(meta.gamma))
    else:
        line = '%-21s ???\n' % name
    comment(line)


def match_symmetry(meta1, meta2):
    return (meta1 and meta2 and
            [a[0] for a in meta1.symmetry.split()] ==
            [a[0] for a in meta2.symmetry.split()])

def calculate_difference_metric(meta1, meta2):
    if not match_symmetry(meta1, meta2):
        return sys.float_info.max
    return sum(abs(a-b) for a,b in zip(meta1.cell, meta2.cell))


def _check_picture_tools():
    ok = True
    coot_path, coot_ver = coot.find_path_and_version()
    if coot_path:
        if coot_ver is None:
            put_error("coot not working(?), no pictures")
            ok = False
        else:
            if "with python" not in coot_ver:
                put_error("coot with Python support is needed")
                ok = False
            if "\n0.6." in coot_ver:
                put_error("coot 0.7+ is needed (0.6 would crash)")
                ok = False
    else:
        put_error("No coot, no pictures")
        ok = False
    if not syspath("render"):
        put_error("No Raster3d, no pictures")
        ok = False
    return ok


def _generate_pictures(wf, opt, fb_job):
    blobs = fb_job.data["blobs"]
    com = fb_job.data["center"]

    # run-coot.py centers on the biggest blob. It uses relative paths -
    # it can be run only from the output directory, but is not affected
    # by moving that directory to different location.
    # There are blobN-coot.py scripts generated below with absolute paths.
    # write coot script (apart from pictures) that centers on the biggest blob
    script_path = os.path.join(wf.output_dir, "run-coot.py")
    script = coot.basic_script(pdb=opt.xyzout, mtz=opt.hklout,
                               center=blobs[0], toward=com)
    open(script_path, "w").write(script)

    # blob images, for now for not more than two blobs
    basenames = []
    for n, b in enumerate(blobs[:2]):
        py_path = os.path.join(wf.output_dir, "blob%d-coot.py" % (n+1))
        with open(py_path, "w") as blob_py:
            d = os.path.abspath(wf.output_dir)
            blob_py.write(coot.basic_script(pdb=os.path.join(d, opt.xyzout),
                                            mtz=os.path.join(d, opt.hklout),
                                            center=blobs[n], toward=com))
        if n != 0:
            # workaround for buggy coot: reloading maps
            script += coot.basic_script(pdb=opt.xyzout, mtz=opt.hklout,
                                        center=b, toward=com)
        rs, names = coot.r3d_script(b, com, blobname="blob%s"%(n+1))
        script += rs
        basenames += names
    try:
        wf.coot_py(script).run()
    except c4.workflow.JobError:
        # check for a possible cause to hint the user
        # (possible workaround: change $HOME to non-existing directory)
        retcode = wf.silently_run_job(wf.coot_py(script_text=""))
        if retcode != 0:
            put_error("coot fails with options: --no-graphics --python",
                      comment="It happens when scripts in .coot or "
                              ".coot-preferences are not compatible\n"
                              "with the --no-graphics mode.")
        raise
    for basename in basenames:
        wf.render_r3d(basename, img_format=opt.img_format).run()
    wf.delete_files([name+".r3d" for name in basenames])


def parse_dimple_commands():
    dstr = ' (default: %(default)s)'
    parser = argparse.ArgumentParser(
                usage='%(prog)s [options...] input.mtz input.pdb output_dir',
                epilog=c4.workflow.commands_help, prog="dimple",
                formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('all_args', nargs='*')
    parser.add_argument('--hklout', metavar='out.mtz', default='final.mtz',
                        help='output mtz file'+dstr)
    parser.add_argument('--xyzout', metavar='out.pdb', default='final.pdb',
                        help='output pdb file'+dstr)
    parser.add_argument('-s', '--summary', action='store_true',
                        help='show refmac summary')
    parser.add_argument('-f', choices=['png', 'jpeg', 'tiff', 'none'],
                        default='png', dest='img_format',
                        help='format of generated images'+dstr)
    parser.add_argument('--weight', metavar='VALUE', type=float,
                        help='refmac matrix weight (default: auto-weight)')
    parser.add_argument('--libin', metavar='CIF',
                        help='ligand descriptions for refmac (LIBIN)')
    parser.add_argument('-R', '--free-r-flags', metavar='MTZ_FILE',
                    help='reference file with all reflections and freeR flags')
    parser.add_argument('-M', '--mr-when-rfree', type=float, default=0.4,
                        metavar='NUM',
                        help='threshold for Molecular Replacement')
    parser.add_argument('-I', '--icolumn', metavar='ICOL',
                        default='IMEAN', help='I column label'+dstr)
    parser.add_argument('--sigicolumn', metavar='SIGICOL',
                        default='SIG<ICOL>', help='SIGI column label'+dstr)
    parser.add_argument('--from-job', metavar='N', type=int, default=0,
                        help=argparse.SUPPRESS)
    parser.add_argument('--version', action='version',
                        version='%(prog)s '+__version__)
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
    # opt.all_args should be one mtz, one or more pdbs and output_dir
    if len(opt.all_args) < 3:
        put_error("At least 3 arguments expected.")
        sys.exit(1)
    opt.output_dir = opt.all_args.pop()
    if opt.output_dir.endswith('.mtz') or opt.output_dir.endswith('.pdb'):
        put_error('The last argument should be output directory')
        sys.exit(1)
    mtz_args = [a for a in opt.all_args if a.lower().endswith('.mtz')]
    if len(mtz_args) != 1:
        put_error("One mtz file should be given.")
        sys.exit(1)
    opt.mtz = mtz_args[0]
    opt.all_args.remove(opt.mtz)
    opt.pdbs = opt.all_args
    for a in opt.pdbs:
        if not a.lower().endswith('.pdb'):
            put_error("unexpected arg (neither mtz nor pdb): %s" % a)
            sys.exit(1)
    if len(opt.pdbs) == 0:
        put_error("At least one pdb file should be given.")
        sys.exit(1)

    # extra checks
    for filename in opt.pdbs + [opt.mtz, opt.free_r_flags, opt.libin]:
        if filename and not os.path.isfile(filename):
            put_error("File not found: " + filename)
            sys.exit(1)
    if os.path.exists(opt.output_dir) and not os.path.isdir(opt.output_dir):
        put_error("Not a directory: " + opt.output_dir)
        sys.exit(1)

    # Since we'll execute programs from opt.output_dir, adjust paths.
    opt.mtz = adjust_path(opt.mtz, opt.output_dir)
    opt.pdbs = [adjust_path(a, opt.output_dir) for a in opt.pdbs]
    if opt.free_r_flags:
        opt.free_r_flags = adjust_path(opt.free_r_flags, opt.output_dir)
    if opt.libin:
        opt.libin = adjust_path(opt.libin, opt.output_dir)

    # the default value of sigicolumn ('SIG<ICOL>') needs substitution
    opt.sigicolumn = opt.sigicolumn.replace('<ICOL>', opt.icolumn)

    return opt


def main():
    if c4.workflow.parse_workflow_commands():
        return

    for necessary_var in ("CCP4", "CCP4_SCR"):
        if necessary_var not in os.environ:
            put_error('$%s not found, giving up' % necessary_var)
            sys.exit(1)
    if not os.path.isdir(os.environ["CCP4_SCR"]):
        put_error('No such directory: $CCP4_SCR, refmac shall not work!')

    options = parse_dimple_commands()

    wf = c4.workflow.Workflow(options.output_dir, from_job=options.from_job)
    start_log(os.path.join(options.output_dir, "dimple.log"),
              output_dir=options.output_dir)
    try:
        dimple(wf=wf, opt=options)
    except c4.workflow.JobError, e: # avoid "as e" for the sake of Py2.4
        put_error(e.msg, comment=e.note)
        wf.pickle_jobs()
        sys.exit(1)
    wf.pickle_jobs()

if __name__ == "__main__":
    main()
