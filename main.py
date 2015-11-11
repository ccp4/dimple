#!/usr/bin/env python

import os
import sys
if sys.version_info[:2] != (2, 7):
    sys.stderr.write("Error. Python 2.7 is required.\n")
    sys.exit(1)
import argparse
if __name__ == "__main__" and __package__ is None:
    sys.path.insert(1,
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dimple import utils
from dimple.utils import comment, put_error
from dimple.cell import match_symmetry
from dimple.mtz import check_freerflags_column, MtzMeta
from dimple.pdb import is_pdb_id, download_pdb, check_hetatm_x
from dimple import workflow
from dimple import coots

__version__ = '2.3.9'


def dimple(wf, opt):
    comment("%8s### Dimple v%s. Problems and suggestions:"
            " ccp4@ccp4.ac.uk ###" % ('', __version__))
    mtz_meta = wf.read_mtz_metadata(opt.mtz)
    mtz_meta.check_col_type(opt.icolumn or 'IMEAN', 'J')
    mtz_meta.check_col_type(opt.sigicolumn, 'Q')
    _comment_summary_line("MTZ (%.1fA)" % mtz_meta.dmax, mtz_meta)
    if opt.dls_naming:
        opt.pdbs = dls_name_filter(opt.pdbs)
    for p in opt.pdbs:
        wf.read_pdb_metadata(p)
    if len(opt.pdbs) > 1:
        comment("\nPDBs in order of similarity (using the first one):")
        opt.pdbs.sort(key=lambda x: calculate_difference_metric(wf.file_info[x],
                                                                mtz_meta))
    utils.log_value("data_file", opt.mtz)
    utils.log_value("pdb_files", opt.pdbs)
    for p in opt.pdbs:
        _comment_summary_line(os.path.basename(p), wf.file_info[p])
    ini_pdb = "ini.pdb"
    wf.copy_uncompressed(opt.pdbs[0], ini_pdb)
    pdb_meta = wf.file_info[opt.pdbs[0]]
    if pdb_meta is None:
        put_error("Failed to read CRYST1 record from the pdb file")
        return
    remove_hetatm = False
    if remove_hetatm or check_hetatm_x(wf.path(ini_pdb), pdb_meta):
        if not remove_hetatm:
            comment("\nHETATM marked as element X would choke many programs.")
        rb_xyzin = "prepared.pdb"
        wf.temporary_files.add(rb_xyzin)
        n_het = wf.remove_hetatm(xyzin=ini_pdb, xyzout=rb_xyzin,
                                 remove_all=remove_hetatm)
        comment("\nRemoved %d HETATM atoms" % n_het)
    else:
        rb_xyzin = ini_pdb
    wf.rwcontents(xyzin=rb_xyzin).run()
    solvent_pct = wf.jobs[-1].data.get('solvent_percent')
    pdb_num_mol = wf.jobs[-1].data.get('num_mol')
    if solvent_pct is None:
        raise RuntimeError("rwcontents could not interpret %s." % rb_xyzin)
    if solvent_pct > 70:
        comment("\nHmm... %.1f%% of solvent looks suspicious" % solvent_pct)
    if abs(wf.jobs[-1].data.get('volume', 0) - pdb_meta.get_volume()) > 10:
        comment("\ndebug: problem when calculating volume?")

    ####### pointless - reindexing #######
    if match_symmetry(mtz_meta, pdb_meta):
        reindexed_mtz = "pointless.mtz"
        wf.temporary_files.add(reindexed_mtz)
        wf.pointless(hklin=opt.mtz, xyzin=rb_xyzin, hklout=reindexed_mtz,
                     keys="TOLERANCE 5").run()
        pointless_data = wf.jobs[-1].data
        alt_reindex = pointless_data.get('alt_reindex')
        if alt_reindex:
            for ar in alt_reindex:
                comment("\n    %-10s CC: %-8.3f cell diff: %.1fA" % (
                        ar['op'], ar['cc'], ar['cell_deviat']))
        else:
            # until recently (2015) pointless didn't print CC for non-ambiguous
            # spacegroups (e.g. C2), but now it always prints
            comment("\n    no good indexing")
            reindexed_mtz = opt.mtz
    else:
        reindexed_mtz = opt.mtz
    reindexed_mtz_meta = wf.read_mtz_metadata(reindexed_mtz)
    if reindexed_mtz_meta.symmetry != mtz_meta.symmetry:
        _comment_summary_line('reindexed MTZ', reindexed_mtz_meta)

    ####### (c)truncate - calculate amplitudes #######
    if (opt.ItoF_prog or opt.icolumn or 'F' not in mtz_meta.columns
                                  or 'SIGF' not in mtz_meta.columns):
        f_mtz = "amplit.mtz"
        wf.temporary_files.add(f_mtz)
        i_sigi_cols = (opt.icolumn or 'IMEAN', opt.sigicolumn)
        if opt.ItoF_prog == 'ctruncate' or (opt.ItoF_prog is None and opt.slow):
            wf.ctruncate(hklin=reindexed_mtz, hklout=f_mtz,
                         colin="/*/*/[%s,%s]" % i_sigi_cols).run()
        else:
            wf.truncate(hklin=reindexed_mtz, hklout=f_mtz,
                        labin="IMEAN=%s SIGIMEAN=%s" % i_sigi_cols,
                        labout="F=F SIGF=SIGF").run()
    else:
        f_mtz = reindexed_mtz

    ####### rigid body - check if model is good for refinement? #######
    refmac_labin_nofree = "FP=F SIGFP=SIGF"
    refmac_labout = ("FC=FC PHIC=PHIC FWT=2FOFCWT PHWT=PH2FOFCWT "
                     "DELFWT=FOFCWT PHDELWT=PHFOFCWT")
    refmac_xyzin = None
    cell_diff = calculate_difference_metric(pdb_meta, reindexed_mtz_meta)
    if cell_diff > 0.1 and opt.mr_when_r < 1:
        comment("\nQuite different unit cells, start from MR.")
    else:
        comment("\nRigid-body refinement with resolution 3.5 A, 10 cycles.")
        wf.temporary_files |= {"refmacRB.pdb", "refmacRB.mtz"}
        # it may fail because of "Disagreement between mtz and pdb"
        wf.refmac5(hklin=f_mtz, xyzin=rb_xyzin,
                   hklout="refmacRB.mtz", xyzout="refmacRB.pdb",
                   labin=refmac_labin_nofree, labout=refmac_labout,
                   libin=None,
                   keys="""refinement type rigidbody resolution 15 3.5
                           rigidbody ncycle 10""").run(may_fail=True)
        # if the error is caused by mtz/pdb disagreement, continue with MR
        if wf.jobs[-1].exit_status != 0:
            if opt.mr_when_r >= 1:
                comment("\nFailed. Would try MR, but it is disabled.")
                return
            comment("\nTry MR.")
        elif not wf.jobs[-1].data.get("overall_r"):
            comment("\nWARNING: unknown R factor, something went wrong.")
            refmac_xyzin = "refmacRB.pdb"
        elif wf.jobs[-1].data["overall_r"] > opt.mr_when_r:
            comment("\nRun MR for R > %g" % opt.mr_when_r)
        else:
            comment("\nNo MR for R < %g" % opt.mr_when_r)
            refmac_xyzin = "refmacRB.pdb"

    ####### phaser/molrep - molecular replacement #######
    if refmac_xyzin is None:
        # pdb_num_mol accounts for strict NCS (MTRIX without iGiven)
        vol_ratio = mtz_meta.asu_volume() / pdb_meta.asu_volume(pdb_num_mol)
        if vol_ratio < 0.66:
            comment("\nasu %.0f%% smaller than in the model"
                    % ((1 - vol_ratio) * 100))
        if opt.MR_prog == 'molrep' or (opt.MR_prog is None and
                                       solvent_pct > 70):
            wf.temporary_files |= {"molrep.pdb", "molrep_dimer.pdb",
                                   "molrep.crd"}
            wf.molrep(f=f_mtz, m=rb_xyzin).run()
            refmac_xyzin = "molrep.pdb"
        else:
            num = 1
            if vol_ratio > 1.5:
                num = int(round(vol_ratio))
                comment("\nSearching %d molecules, asu is %.1f x larger "
                        "than in the model" % (num, vol_ratio))
            wf.temporary_files |= {"phaser.1.pdb", "phaser.1.mtz"}
            wf.phaser_auto(hklin=f_mtz,
                           labin="F = F SIGF = SIGF",
                           model=dict(pdb=rb_xyzin, identity=100, num=num),
                           solvent_percent=solvent_pct,
                           sg_alt="ALL",
                           root='phaser').run(may_fail=True)
            phaser_data = wf.jobs[-1].data
            if (wf.jobs[-1].exit_status != 0 or
                    phaser_data['status'] == 'Sorry - No solution'):
                comment("\nGiving up.")
                return
            if phaser_data['status'].endswith('...'):
                solu_set = wf.get_phaser_solu_set(root='phaser')
                if solu_set:
                    comment("\n..." + solu_set[len(phaser_data['status'])-3:])
                    phaser_data['status'] = solu_set
            if phaser_data.get('partial_solution'):
                comment("\nOnly partial solution found")
            if phaser_data['SG'] != reindexed_mtz_meta.symmetry:
                comment("\nSpacegroup changed to %s" % phaser_data['SG'])
            refmac_xyzin = "phaser.1.pdb"
            f_mtz = "phaser.1.mtz"

    if False:
        wf.findwaters(pdbin=refmac_xyzin, hklin=f_mtz,
                      f="FC", phi="PHIC", pdbout="prepared_wat.pdb", sigma=2)
        refmac_xyzin = "prepared_wat.pdb"

    ####### adding free-R flags #######
    f_mtz_cols = wf.read_mtz_metadata(f_mtz).columns.keys()
    cad_reso = opt.reso or (reindexed_mtz_meta.dmax - MtzMeta.d_eps)
    if opt.free_r_flags:
        free_mtz = opt.free_r_flags
        free_col = check_freerflags_column(wf.path(free_mtz),
                                           expected_symmetry=pdb_meta)
        comment("\nFree-R flags from the %s file, column %s." %
                (("reference" if free_mtz != opt.mtz else 'input'), free_col))
    else:
        free_col = 'FreeR_flag'
        if free_col in f_mtz_cols:
            comment("\nReplace free-R flags")
        else:
            comment("\nGenerate free-R flags")
        free_mtz = "free.mtz"
        wf.temporary_files |= {"unique.mtz", free_mtz}
        # CCP4 freerflag uses always the same pseudo-random sequence by default
        if opt.seed_freerflag:
            wf.unique(hklout="unique.mtz", ref=reindexed_mtz_meta,
                      resolution=cad_reso).run()
            wf.freerflag(hklin="unique.mtz", hklout=free_mtz, keys="SEED").run()
        else:
            comment(" (repeatably)")
            # here we'd like to have always the same set of free-r flags
            # for given PDB file. That's why it's using cell parameters
            # from pdb not mtz and why always the same resolution is set.
            wf.unique(hklout="unique.mtz", ref=pdb_meta, resolution=1.0).run()
            wf.freerflag(hklin="unique.mtz", hklout=free_mtz).run()

    if free_mtz == opt.mtz and opt.reso is None:
        prepared_mtz = f_mtz
    else:
        prepared_mtz = "prepared.mtz"
        wf.temporary_files.add(prepared_mtz)
        wf.cad(data_in=[(f_mtz, [c for c in f_mtz_cols if c != free_col]),
                        (free_mtz, [free_col])],
               hklout=prepared_mtz,
               keys=["sysab_keep",  # does it matter?
                     "reso overall 1000.0 %g" % cad_reso]).run()
    freerflag_missing = wf.count_mtz_missing(prepared_mtz, free_col)
    if freerflag_missing:
        wf.freerflag(hklin=prepared_mtz, hklout="prepared2.mtz",
                     keys="COMPLETE FREE="+free_col,
                     parser=" (again, for %d refl. more)" % freerflag_missing
                    ).run()
        prepared_mtz = "prepared2.mtz"
        wf.temporary_files.add(prepared_mtz)

    ####### refinement #######
    if opt.weight:
        refmac_weight = "matrix %f" % opt.weight
    else:
        refmac_weight = "auto"
    restr_ref_keys = """\
     make newligand continue
     refinement type restrained
     weight %s
     """ % refmac_weight
    refmac_labin = "%s FREE=%s" % (refmac_labin_nofree, free_col)
    if opt.jelly:
        comment("\nJelly-body refinement, %d cycles." % opt.jelly)
        wf.temporary_files |= {"jelly.pdb", "jelly.mtz"}
        wf.refmac5(hklin=prepared_mtz, xyzin=refmac_xyzin,
                   hklout="jelly.mtz", xyzout="jelly.pdb",
                   labin=refmac_labin, labout=refmac_labout, libin=opt.libin,
                   keys=restr_ref_keys+"ridge distance sigma 0.01\n"
                                       "make hydrogen no\n"
                                       "ncycle %d" % opt.jelly).run()
        comment(_refmac_rms_line(wf.jobs[-1].data))
        refmac_xyzin = "jelly.pdb"
    comment("\nFinal restrained refinement, %d cycles." % opt.restr_cycles)
    restr_job = wf.refmac5(hklin=prepared_mtz, xyzin=refmac_xyzin,
                 hklout=opt.hklout, xyzout=opt.xyzout,
                 labin=refmac_labin, labout=refmac_labout, libin=opt.libin,
                 keys=restr_ref_keys+("ncycle %d" % opt.restr_cycles)).run()
    comment(_refmac_rms_line(restr_job.data))
    # if that run is repeated with --from-step it's useful to compare Rfree
    if wf.from_job > 0 and wf.from_job <= len(wf.jobs): # from_job is 1-based
        prev = [j for j in wf.repl_jobs if j.name == restr_job.name]
        if prev and prev[0].data and "free_r" in prev[0].data:
            comment("\nPreviously:  R/Rfree %.4f/%.4f  Rfree change: %+.4f" % (
                    prev[0].data["overall_r"], prev[0].data["free_r"],
                    restr_job.data["free_r"] - prev[0].data["free_r"]))

    ####### check blobs and finish #######
    fb_job = wf.find_blobs(opt.hklout, opt.xyzout, sigma=0.8).run()
    _generate_scripts_and_pictures(wf, opt, fb_job)


def _refmac_rms_line(data):
    items = []
    for name, col in [('bond', 'rmsBOND'), ('angle', 'rmsANGL'),
                      ('chiral', 'rmsCHIRAL')]:
        d = data.get(col)
        if d:
            items.append('   %s %.2f -> %.2f' % (name, d[0], d[-1]))
    return '\n  RMS:' + ''.join(items)


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


def calculate_difference_metric(meta1, meta2):
    match = match_symmetry(meta1, meta2)
    # wrong or corrupted file (no CRYST1) is worse than non-matching file
    if match is None:
        return sys.float_info.max
    if match is False:
        return sys.float_info.max / 2
    #return sum(abs(a-b) for a,b in zip(meta1.cell, meta2.cell))
    return meta1.to_standard().max_shift_in_mapping(meta2.to_standard())


def _check_picture_tools():
    ok = True
    coot_path, coot_ver = coots.find_path_and_version()
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
    if not utils.syspath("render"):
        put_error("No Raster3d, no pictures")
        ok = False
    return ok

# chmod +x
def _make_executable(path):
    mode = os.stat(path).st_mode
    os.chmod(path, mode | ((mode & 0444) >> 2))

def _generate_scripts_and_pictures(wf, opt, fb_job):
    blobs = fb_job.data["blobs"]
    if not blobs:
        comment("\nUnmodelled blobs not found.")
    elif opt.img_format != 'none' and _check_picture_tools():
        if len(blobs) == 1:
            comment("\nRendering density blob at (%.1f, %.1f, %.1f)" %
                    blobs[0])
        else:
            comment("\nRendering 2 largest blobs: at (%.1f, %.1f, %.1f) "
                    "and at (%.1f, %.1f, %.1f)" % (blobs[0]+blobs[1]))
    com = fb_job.data["center"]

    # run-coot.py centers on the biggest blob. It uses relative paths -
    # it can be run only from the output directory, but is not affected
    # by moving that directory to different location.
    # There are blobN-coot.py scripts generated below with absolute paths.
    # write coot script (apart from pictures) that centers on the biggest blob
    script_path = os.path.join(wf.output_dir, "run-coot.py")
    script = coots.basic_script(pdb=opt.xyzout, mtz=opt.hklout,
                               center=(blobs and blobs[0]), toward=com)
    open(script_path, "w").write(script)

    # blob images, for now for not more than two blobs
    for n, b in enumerate(blobs[:2]):
        py_path = os.path.join(wf.output_dir, "blob%d-coot.py" % (n+1))
        with open(py_path, "w") as blob_py:
            d = os.path.abspath(wf.output_dir)
            blob_py.write(coots.basic_script(pdb=os.path.join(d, opt.xyzout),
                                             mtz=os.path.join(d, opt.hklout),
                                             center=blobs[n], toward=com))
    # coot.sh - one-line script for convenience
    if blobs:
        coot_sh_text = '{coot} --no-guano {out}/blob1-coot.py\n'
    else:
        coot_sh_text = '{coot} --no-guano {out}/final.mtz {out}/final.pdb\n'
    coot_sh_path = os.path.join(wf.output_dir, "coot.sh")
    try:
        with open(coot_sh_path, 'w') as f:
            f.write(coot_sh_text.format(coot=coots.find_path(),
                                        out=wf.output_dir))
        _make_executable(coot_sh_path)
    except (IOError, OSError) as e:
        put_error(e)

    if opt.img_format == 'none':
        return

    script = ''
    basenames = []
    # as a workaround for buggy coot the maps are reloaded for each blob
    for n, b in enumerate(blobs[:2]):
        script += coots.basic_script(pdb=opt.xyzout, mtz=opt.hklout,
                                     center=b, toward=com)
        rs, names = coots.r3d_script(b, com, blobname="blob%s"%(n+1))
        script += rs
        basenames += names
    coot_job = wf.coot_py(script)
    try:
        coot_job.run()
    except workflow.JobError:
        # check for a possible cause to hint the user
        # (possible workaround: change $HOME to non-existing directory)
        if utils.silently_run(coot_job.args, cwd=wf.output_dir)[0] != 0:
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
            job.run(show_progress=False, new_line=False)
    wf.delete_files([name+".r3d" for name in basenames])


def parse_dimple_commands(args):
    dstr = ' (default: %(default)s)'
    parser = argparse.ArgumentParser(
                usage='%(prog)s [options...] input.mtz input.pdb output_dir',
                epilog=workflow.commands_help, prog="dimple",
                formatter_class=argparse.RawDescriptionHelpFormatter)
    # positional args can be separated by options, but not after the 3rd one
    # see http://bugs.python.org/issue15112 , http://bugs.python.org/issue14191
    parser.add_argument('pos_arg1')
    parser.add_argument('pos_arg2')
    parser.add_argument('pos_arg3')
    parser.add_argument('more_args', nargs='*')
    group1 = parser.add_argument_group('most commonly used options')
    group1.add_argument('-s', '--slow', action='count',
                        help='more refinement, etc. (can be used 2x)')
    group1.add_argument('-M', '--mr-when-r', type=float, default=0.4,
                        metavar='NUM',
                        help='threshold for Molecular Replacement'+dstr)
    group2 = parser.add_argument_group('options contolling input/output')
    group2.add_argument('-I', '--icolumn', metavar='ICOL',
                        help='I column label (default: IMEAN)')
    group2.add_argument('--sigicolumn', metavar='SIGICOL',
                        default='SIG<ICOL>', help='SIGI column label'+dstr)
    group2.add_argument('--libin', metavar='CIF',
                        help='ligand descriptions for refmac (LIBIN)')
    group2.add_argument('-R', '--free-r-flags', metavar='MTZ_FILE',
                        help='file with freeR flags '
                             '("-" = use flags from data mtz)')
    group2.add_argument('--hklout', metavar='out.mtz', default='final.mtz',
                        help='output mtz file'+dstr)
    group2.add_argument('--xyzout', metavar='out.pdb', default='final.pdb',
                        help='output pdb file'+dstr)
    group2.add_argument('-f', choices=['png', 'jpeg', 'tiff', 'none'],
                        default='png', dest='img_format',
                        help='format of generated images'+dstr)
    group2.add_argument('--no-cleanup', dest='cleanup', action='store_false',
                        help='leave intermediate files')
    group2.add_argument('--cleanup', action='store_true',
                        help=argparse.SUPPRESS)  # obsolete

    group3 = parser.add_argument_group('options customizing the run')
    group3.add_argument('--jelly', metavar='N_ITER', type=int,
                    help='run refmac jelly-body before the final refinement')
    group3.add_argument('--reso', type=float, help='limit the resolution [A]')
    group3.add_argument('--restr-cycles', metavar='N', type=int,
                        help='cycles of refmac final refinement (default: 8)')
    group3.add_argument('--weight', metavar='VALUE', type=float,
                        help='refmac matrix weight (default: auto-weight)')

    group3.add_argument('--MR-prog', choices=['phaser', 'molrep'],
                        help='Molecular Replacement program')
    group3.add_argument('--ItoF-prog', choices=['truncate', 'ctruncate'],
                        help='program to calculate amplitudes')
    group3.add_argument('--seed-freerflag', action='store_true',
                        help=argparse.SUPPRESS)
    group3.add_argument('--dls-naming', action='store_true',
                        help=argparse.SUPPRESS)
    group3.add_argument('--from-step', metavar='N', type=int, default=0,
                        help=argparse.SUPPRESS)
    parser.add_argument('--version', action='version',
                        version='%(prog)s '+__version__)
    # customize usage message: get rid of 'positional arguments',
    # rename default 'optional arguments' and shift it to the end.
    # pylint: disable=protected-access
    default_group = parser._action_groups[1]
    default_group.title = 'other options'
    parser._action_groups = parser._action_groups[2:] + [default_group]

    # special mode for compatibility with ccp4i
    legacy_args = {"HKLIN": "", "XYZIN": "",
                   "HKLOUT": "--hklout", "XYZOUT": "--xyzout"}
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
    # special mode for re-running jobs
    if all_args[0] == 'rerun':
        if os.path.isdir(all_args[1]):
            logfile = os.path.join(all_args[1], 'dimple.log')
        else:
            logfile = all_args[1]
        old_wf = utils.read_section_from_log(logfile, 'workflow')
        try:
            old_dir = os.path.join(old_wf['cwd'], old_wf['output_dir'])
            old_pdb = os.path.join(old_dir, 'ini.pdb')
            if 'data_file' not in old_wf:  # temporary, to be removed soon
                old_mtz_arg = [a for a in old_wf['args'].split()
                               if a.endswith('.mtz')][0]
                old_wf['data_file'] = os.path.join(old_wf['cwd'], old_mtz_arg)
            old_mtz = os.path.join(old_dir, old_wf['data_file'])
        except (TypeError, KeyError):
            put_error('Reading logfile failed')
            sys.exit(1)
        all_args[0:2] = [old_pdb, old_mtz]

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
    if opt.seed_freerflag and opt.free_r_flags:
        put_error("Option --seed-freerflag and --free-r-flags"
                  " don't make sense together")
        sys.exit(1)
    if opt.free_r_flags == '-':
        opt.free_r_flags = opt.mtz

    # extra checks
    for filename in opt.pdbs + [opt.mtz, opt.free_r_flags, opt.libin]:
        if filename and not os.path.isfile(filename):
            put_error("File not found: " + filename)
            sys.exit(1)
    if os.path.exists(opt.output_dir) and not os.path.isdir(opt.output_dir):
        put_error("Not a directory: " + opt.output_dir)
        sys.exit(1)

    # Since we'll execute programs from opt.output_dir, adjust paths.
    opt.mtz = utils.adjust_path(opt.mtz, opt.output_dir)
    opt.pdbs = [utils.adjust_path(a, opt.output_dir) for a in opt.pdbs]
    if opt.free_r_flags:
        opt.free_r_flags = utils.adjust_path(opt.free_r_flags, opt.output_dir)
    if opt.libin:
        opt.libin = utils.adjust_path(opt.libin, opt.output_dir)

    # the default value of sigicolumn ('SIG<ICOL>') needs substitution
    opt.sigicolumn = opt.sigicolumn.replace('<ICOL>', opt.icolumn or 'IMEAN')

    if opt.restr_cycles is None:
        opt.restr_cycles = (12 if opt.slow else 8)
    if opt.jelly is None and opt.slow:
        pass # opt.jelly = 50  # isn't it too slow?

    return opt

def dls_name_filter(pdbs):
    # Filename matching used in Diamond synchrotron. PDB filenames
    # are matched against the current (!) directory.
    # It's more relaxed than in solve_o_matic's select_pdb.py:
    # case-insensitive and ignoring non-alphanumeric characters.
    pattern = ''.join(a for a in os.getcwd().lower()
                      if a.isalnum() or a == '/')
    def token(arg):
        part = os.path.basename(arg).split('.')[0]
        return ''.join(a for a in part.lower() if a.isalnum())
    matched_pdbs = [arg for arg in pdbs if token(arg) in pattern]
    if matched_pdbs != pdbs:
        comment("\n%d of %d PDBs have filenames matching data directory"
                % (len(matched_pdbs), len(pdbs)))
    return matched_pdbs


def main(args):
    if workflow.parse_workflow_commands():
        return

    for necessary_var in ("CCP4", "CCP4_SCR"):
        if necessary_var not in os.environ:
            put_error('$%s not found, giving up' % necessary_var)
            sys.exit(1)
    if not os.path.isdir(os.environ["CCP4_SCR"]):
        put_error('No such directory: $CCP4_SCR, refmac shall not work!')

    options = parse_dimple_commands(args)

    wf = workflow.Workflow(options.output_dir, from_job=options.from_step)
    utils.start_log(os.path.join(options.output_dir, "dimple.log"),
                    output_dir=options.output_dir)
    utils.log_value("version", __version__)
    utils.start_log_screen(os.path.join(options.output_dir, "screen.log"))
    try:
        dimple(wf=wf, opt=options)
        exit_status = 0
    except workflow.JobError as e:
        put_error(e.msg, comment=e.note)
        try:
            utils.report_disk_space([wf.output_dir, os.getenv("CCP4_SCR")])
        except KeyboardInterrupt:
            comment("\nok, exiting...")
        exit_status = 1
    except RuntimeError as e:
        put_error(e)
        exit_status = 1
    finally:
        comment("\n")
    if options.cleanup:
        wf.delete_files(wf.temporary_files)
    wf.options = options
    wf.dump_pickle()
    return exit_status

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
