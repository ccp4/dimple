#!/usr/bin/env python

import os
import sys
if sys.version_info[:2] != (2, 7):
    sys.stderr.write("Error. Python 2.7 is required.\n")
    sys.exit(1)
import argparse
import c4.utils
from c4.utils import comment, put_error, adjust_path
from c4.mtz import check_freerflags_column, get_num_missing
from c4.pdb import is_pdb_id, download_pdb
import c4.workflow
from c4 import coot

__version__ = '2.1.1'

def dimple(wf, opt):
    comment("%8s### Dimple v%s. Problems and suggestions:"
            " ccp4@ccp4.ac.uk ###" % ('', __version__))
    mtz_meta = wf.read_mtz_metadata(opt.mtz)
    wf.file_info[opt.mtz] = mtz_meta
    mtz_meta.check_col_type(opt.icolumn, 'J')
    mtz_meta.check_col_type(opt.sigicolumn, 'Q')
    _comment_summary_line("MTZ", mtz_meta)
    if opt.dls_filename_matching:
        # the same filtering as in solve_o_matic's select_pdb.py in Diamond LS
        match_pdbs = [arg for arg in opt.pdbs
                      if os.path.basename(arg).split('.')[0] in os.getcwd()]
        if match_pdbs != opt.pdbs:
            comment("\n%d of %d PDBs have filenames matching data directory"
                    % (len(match_pdbs), len(opt.pdbs)))
            if match_pdbs:
                opt.pdbs = match_pdbs
            else:
                comment("\nWARNING: using PDBs with non-matching filenames.")
    for p in opt.pdbs:
        wf.file_info[p] = wf.read_pdb_metadata(p)
    if len(opt.pdbs) > 1:
        comment("\nPDBs in order of similarity (using the first one):")
        opt.pdbs.sort(key=lambda x: calculate_difference_metric(wf.file_info[x],
                                                                mtz_meta))
    for p in opt.pdbs:
        _comment_summary_line(os.path.basename(p), wf.file_info[p])
    ini_pdb = "ini.pdb"
    wf.copy_uncompressed(opt.pdbs[0], ini_pdb)
    pdb_meta = wf.file_info[opt.pdbs[0]]
    wf.pointless(hklin=opt.mtz, xyzin=ini_pdb, hklout="pointless.mtz",
                 keys="TOLERANCE 2").run()
    pointless_data = wf.jobs[-1].data
    alt_reindex = pointless_data.get('alt_reindex')
    if alt_reindex:
        for ar in alt_reindex:
            comment("\n    %-10s CC: %-8s cell_deviation: %s" % (
                    ar['op'], ar['cc'], ar['cell_deviat']))
    else:
        # until recently (2015) pointless did't print CC for non-ambiguous
        # spacegroups (e.g. C2), but now it always prints (?)
        comment("\n    ooops, no good indexing")
    pointless_mtz_meta = wf.read_mtz_metadata("pointless.mtz")
    if pointless_mtz_meta.symmetry != mtz_meta.symmetry:
        _comment_summary_line('reindexed MTZ', pointless_mtz_meta)
    #comment("\nCalculate structure factor amplitudes")
    if opt.ItoF_prog == 'truncate':
        wf.truncate(hklin="pointless.mtz", hklout="truncate.mtz",
                  labin="IMEAN=%s SIGIMEAN=%s" % (opt.icolumn, opt.sigicolumn),
                  labout="F=F SIGF=SIGF").run()
    else:
        wf.ctruncate(hklin="pointless.mtz", hklout="truncate.mtz",
                     colin="/*/*/[%s,%s]" % (opt.icolumn, opt.sigicolumn)).run()

    if opt.free_r_flags:
        free_mtz = opt.free_r_flags
        free_col = check_freerflags_column(wf.path(free_mtz),
                                           expected_symmetry=pdb_meta.symmetry)
        comment("\nFree-R flags from the reference file, column %s." % free_col)
    else:
        comment("\nGenerate free-R flags")
        free_mtz = "free.mtz"
        # CCP4 freerflag uses always the same pseudo-random sequence by default
        if opt.seed_freerflag:
            wf.unique(hklout="unique.mtz",
                      cell=pointless_data['output_cell'],
                      symmetry=pdb_meta.symmetry,
                      resolution=mtz_meta.dmax,
                      labout="F=F_UNIQUE SIGF=SIGF_UNIQUE").run()
            wf.freerflag(hklin="unique.mtz", hklout=free_mtz, keys="SEED").run()
        else:
            # here we'd like to have always the same set of free-r flags
            # for given PDB file. That's why it's MTZ-agnostic.
            wf.unique(hklout="unique.mtz",
                      cell=pdb_meta.cell,
                      symmetry=pdb_meta.symmetry,
                      resolution=1.0, # somewhat arbitrary limit
                      labout="F=F_UNIQUE SIGF=SIGF_UNIQUE").run()
            wf.freerflag(hklin="unique.mtz", hklout=free_mtz).run()
        free_col = 'FreeR_flag'

    prepared_mtz = "prepared.mtz"
    # TODO: add SYSAB_KEEP key if spacegroup is to be searched
    wf.cad(hklin=["truncate.mtz", free_mtz], hklout=prepared_mtz,
           keys="""labin file 1 ALL
                   labin file 2 E1=%s
                   reso file 2 1000.0 %g
                   """ % (free_col, mtz_meta.dmax)).run()
    freerflag_missing = get_num_missing(wf.path(prepared_mtz), free_col)
    if freerflag_missing:
        comment("\nHmmm, missing free-R flags for %d reflections. Adding."
                % freerflag_missing)
        wf.freerflag(hklin=prepared_mtz, hklout="prepared2.mtz",
                     keys="COMPLETE FREE="+free_col).run()
        prepared_mtz = "prepared2.mtz"
    if False:
        rb_xyzin = "prepared_nohet.pdb"
        n_het = wf.remove_hetatm(xyzin=ini_pdb, xyzout=rb_xyzin)
        comment("\nRemoved %s atoms marked as HETATM in pdb." % n_het)
    else:
        rb_xyzin = ini_pdb

    refmac_labin = "FP=F SIGFP=SIGF FREE=%s" % free_col
    refmac_labout = ("FC=FC PHIC=PHIC FWT=2FOFCWT PHWT=PH2FOFCWT "
                     "DELFWT=FOFCWT PHDELWT=PHFOFCWT")

    refmac_xyzin = None
    cell_diff = calculate_difference_metric(pdb_meta, pointless_mtz_meta)
    if cell_diff > 0.25 and opt.mr_when_rfree < 1:
        comment("\nQuite different unit cells, start from MR.")
    else:
        comment("\nRigid-body refinement with resolution 3.5 A, 10 cycles.")
        wf.refmac5(hklin=prepared_mtz, xyzin=rb_xyzin,
                   hklout="refmacRB.mtz", xyzout="refmacRB.pdb",
                   labin=refmac_labin, labout=refmac_labout, libin=None,
                   keys="""refinement type rigidbody resolution 15 3.5
                           scale type simple lssc anisotropic experimental
                           solvent yes vdwprob 1.4 ionprob 0.8 mshrink 0.8
                           rigidbody ncycle 10""").run()

        if "free_r" not in wf.jobs[-1].data:
            comment("\nWARNING: unknown free_r, something went wrong.")
            refmac_xyzin = "refmacRB.pdb"
        elif wf.jobs[-1].data["free_r"] > opt.mr_when_rfree:
            comment("\nRun MR for R_free > %g" % opt.mr_when_rfree)
        else:
            comment("\nNo MR for R_free < %g" % opt.mr_when_rfree)
            refmac_xyzin = "refmacRB.pdb"

    if refmac_xyzin is None:
        if opt.MR_prog == 'molrep':
            wf.molrep(f=prepared_mtz, m=rb_xyzin).run()
            refmac_xyzin = "molrep.pdb"
        else:
            # FIXME: should be the ratio of ASU
            vol_ratio = mtz_meta.get_volume() / pdb_meta.get_volume()
            # FIXME account for strict NCS (MTRIX records without iGiven)
            num = max(int(round(vol_ratio)), 1)
            if num != 1:
                comment("\nSearching %d molecules, mtz cell %.1f x larger "
                        "than model" % (num, vol_ratio))
            wf.phaser_auto(hklin=prepared_mtz,
                      #labin="I = %s SIGI = %s" % (opt.icolumn, opt.sigicolumn),
                      labin="F = F SIGF = SIGF",
                      sg_alt="ALL",
                      model=dict(pdb=rb_xyzin, identity=100, num=num),
                      solvent_percent=_get_solvent_percent(wf, rb_xyzin),
                      root='phaser').run()
            phaser_data = wf.jobs[-1].data
            if phaser_data['status'].startswith('Sorry'):
                return
            if phaser_data['SG'] != pointless_mtz_meta.symmetry:
                comment("\nSpacegroup changed to %s" % phaser_data['SG'])
            refmac_xyzin = "phaser.1.pdb"
            prepared_mtz = "phaser.1.mtz"

    if False:
        wf.findwaters(pdbin=refmac_xyzin, hklin="refmacRB.mtz",
                      f="FC", phi="PHIC", pdbout="prepared_wat.pdb", sigma=2)
        refmac_xyzin = "prepared_wat.pdb"

    comment("\nFinal restrained refinement, 8 cycles.")
    if opt.weight:
        refmac_weight = "matrix 0.2"
    else:
        refmac_weight = "auto"
    restr_job = wf.refmac5(hklin=prepared_mtz, xyzin=refmac_xyzin,
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
            comment("\nPreviously:  R/Rfree %.4f/%.4f  Rfree change: %+.4f" % (
                    prev[0].data["overall_r"], prev[0].data["free_r"],
                    restr_job.data["free_r"] - prev[0].data["free_r"]))

    fb_job = wf.find_blobs(opt.hklout, opt.xyzout, sigma=0.8).run()
    if opt.img_format != 'none':
        blobs = fb_job.data["blobs"]
        if blobs:
            if len(blobs) == 1:
                comment("\nRendering density blob at (%.1f, %.1f, %.1f)" %
                        blobs[0])
            else:
                comment("\nRendering 2 largest blobs: at (%.1f, %.1f, %.1f) "
                        "and at (%.1f, %.1f, %.1f)" % (blobs[0]+blobs[1]))
            if _check_picture_tools():
                _generate_pictures(wf, opt, fb_job)
        else:
            comment("\nUnmodelled blobs not found.")


def _comment_summary_line(name, meta):
    def angle(x):
        if x == 90.: return '90'
        else:        return str(x)
    if meta:
        line = '\n%-21s %-12s (%.2f, %.2f, %.2f,  %s, %s, %s)' % (
                name, meta.symmetry, meta.a, meta.b, meta.c,
                angle(meta.alpha), angle(meta.beta), angle(meta.gamma))
    else:
        line = '\n%-21s ???' % name
    comment(line)


# returns only part of the script
def _get_solvent_percent(wf, pdb):
    rw = wf.rwcontents(pdb)
    if wf.silently_run_job(rw) != 0:
        raise RuntimeError("rwcontents of %s failed." % pdb)
    def val(key):
        for line in rw.out.lines:
            if key in line:
                return float(line[line.find(key)+len(key):])
    Vm = val('The Matthews Coefficient is :')
    #mw = val('Molecular Weight of protein:')
    #vol = val('Cell volume:')
    #data_num = int(round(mtz_meta.get_volume() / mw / Vm))
    if not Vm:
        raise RuntimeError("rwcontents could not interpret %s." % pdb)
    # 1.23 is used in phaser/src/Composition.cc
    return (1 - 1.23/Vm) * 100


def match_symmetry(meta1, meta2):
    if not meta1 or not meta2:
        return None
    def sig(sym):
        first_chars = [a[0] for a in sym.split()]
        s = first_chars[0] + ''.join(sorted(first_chars[1:]))
        if s == 'I112': # I2 is equivalent to C2
            return 'C112'
        return s
    return sig(meta1.symmetry) == sig(meta2.symmetry)

def calculate_difference_metric(meta1, meta2):
    if not match_symmetry(meta1, meta2):
        return sys.float_info.max
    #return sum(abs(a-b) for a,b in zip(meta1.cell, meta2.cell))
    return meta1.to_standard().max_shift_in_mapping(meta2.to_standard())


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
    if not c4.utils.syspath("render"):
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
    for n, basename in enumerate(basenames):
        job = wf.render_r3d(basename, img_format=opt.img_format)
        if n % 3 == 0:
            job.run()
        else: # minimal output
            wf.run_job(job, show_progress=False, new_line=False)
    wf.delete_files([name+".r3d" for name in basenames])


def parse_dimple_commands(args):
    dstr = ' (default: %(default)s)'
    parser = argparse.ArgumentParser(
                usage='%(prog)s [options...] input.mtz input.pdb output_dir',
                epilog=c4.workflow.commands_help, prog="dimple",
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # positional args can be separated by options, but not after the 3rd one
    # see http://bugs.python.org/issue15112 , http://bugs.python.org/issue14191
    parser.add_argument('pos_arg1')
    parser.add_argument('pos_arg2')
    parser.add_argument('pos_arg3')
    parser.add_argument('more_args', nargs='*')
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
                        help='threshold for Molecular Replacement'+dstr)
    parser.add_argument('--MR-prog', choices=['phaser', 'molrep'],
                        default='phaser',
                        help='Molecular Replacement program'+dstr)
    parser.add_argument('-I', '--icolumn', metavar='ICOL',
                        default='IMEAN', help='I column label'+dstr)
    parser.add_argument('--sigicolumn', metavar='SIGICOL',
                        default='SIG<ICOL>', help='SIGI column label'+dstr)
    parser.add_argument('--ItoF-prog', choices=['truncate', 'ctruncate'],
                        default='truncate',
                        help='program to calculate amplitudes'+dstr)
    parser.add_argument('--cleanup', action='store_true',
                        help='remove intermediate files on exit')
    parser.add_argument('--seed-freerflag', action='store_true',
                        help=argparse.SUPPRESS)
    parser.add_argument('--dls-filename-matching', action='store_true',
                        help=argparse.SUPPRESS)
    parser.add_argument('--from-job', metavar='N', type=int, default=0,
                        help=argparse.SUPPRESS)
    parser.add_argument('--version', action='version',
                        version='%(prog)s '+__version__)
    # get rid of 'positional arguments' in the usage method
    parser._action_groups[:1] = []  # pylint: disable=protected-access

    # special mode for compatibility with ccp4i
    legacy_args = { "HKLIN": "", "XYZIN": "",
                    "HKLOUT": "--hklout", "XYZOUT": "--xyzout" }
    if len(args) == 8 and args[0] in legacy_args:
        args = [legacy_args.get(a) or a
                for a in args if legacy_args.get(a) != ""]
        output_dir = os.path.join(os.environ["CCP4_SCR"], "dimple_out")
        args.append(output_dir)

    opt = parser.parse_args(args)
    all_args = [opt.pos_arg1, opt.pos_arg2, opt.pos_arg3] + opt.more_args
    # all_args should be one mtz, one or more pdbs and output_dir
    opt.output_dir = all_args.pop()
    if (opt.output_dir.endswith('.mtz') or opt.output_dir.endswith('.pdb')
            or opt.output_dir.endswith('.gz')):
        put_error('The last argument should be output directory')
        sys.exit(1)
    mtz_args = [a for a in all_args if a.lower().endswith('.mtz')]
    if len(mtz_args) != 1:
        put_error("One mtz file should be given.")
        sys.exit(1)
    opt.mtz = mtz_args[0]
    all_args.remove(opt.mtz)
    opt.pdbs = all_args
    for n, a in enumerate(opt.pdbs):
        if is_pdb_id(a):
            opt.pdbs[n] = download_pdb(a, opt.output_dir)
        elif not (a.lower().endswith('.pdb') or a.lower().endswith('.pdb.gz')):
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


def main(args):
    if c4.workflow.parse_workflow_commands():
        return

    for necessary_var in ("CCP4", "CCP4_SCR"):
        if necessary_var not in os.environ:
            put_error('$%s not found, giving up' % necessary_var)
            sys.exit(1)
    if not os.path.isdir(os.environ["CCP4_SCR"]):
        put_error('No such directory: $CCP4_SCR, refmac shall not work!')

    options = parse_dimple_commands(args)

    wf = c4.workflow.Workflow(options.output_dir, from_job=options.from_job)
    c4.utils.start_log(os.path.join(options.output_dir, "dimple.log"),
                       output_dir=options.output_dir)
    c4.utils.log_value("version", __version__)
    c4.utils.start_log_screen(os.path.join(options.output_dir, "screen.log"))
    try:
        dimple(wf=wf, opt=options)
        comment("\n")
        if options.cleanup:
            wf.delete_files(["pointless.mtz", "truncate.mtz", "unique.mtz",
                             "free.mtz", "prepared.mtz", "prepared2.mtz",
                             "refmacRB.mtz", "refmacRB.pdb",
                             "molrep.pdb", "molrep_dimer.pdb", "molrep.crd",
                             "phaser.1.pdb", "phaser.1.mtz"])
    except c4.workflow.JobError, e: # avoid "as e" for the sake of Py2.4
        put_error(e.msg, comment=e.note)
        wf.pickle_jobs()
        sys.exit(1)
    except RuntimeError, e:
        put_error(e)
        wf.pickle_jobs()
        sys.exit(1)
    wf.options = options
    wf.pickle_jobs()

if __name__ == "__main__":
    main(sys.argv[1:])
