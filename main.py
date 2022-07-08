#!/usr/bin/env python

#  Copyright 2013-2016 Diamond Light Source Ltd
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.

import os
import sys
import argparse
if __name__ == '__main__' and __package__ is None:
    sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__
                                                                       ))))
from dimple import utils
from dimple.utils import comment, put_error
from dimple.cell import Cell, match_symmetry, calculate_difference
from dimple.mtz import check_freerflags_column, MtzMeta, DEFAULT_FREE_COLS
from dimple.pdb import is_pdb_id, download_pdb, check_hetatm_x
from dimple import workflow
from dimple import coots
from dimple import contaminants

__version__ = '2.6.2'

PROG = 'dimple'
USAGE_SHORT = '%s [options...] input.mtz input.pdb output_dir' % PROG

# sometimes provided models are incomplete, be suspicious above this solvent%
HIGH_SOLVENT_PCT = 75
# Sometimes provided models are too big - then all chains are put into
# a single ensemble for MR. The first threshold always triggers ensembling,
# the second one only if the asu volume of the data and model is different.
VERY_LOW_SOLVENT_PCT = 10
LOW_SOLVENT_PCT = 30
# do not search blobs if the model is too bad
BAD_FINAL_RFREE = 0.5
# do not check the list of contaminants if the model is good
GOOD_FINAL_RFREE = 0.4

def dimple(wf, opt):
    comment('     ### Dimple v%s. Problems and suggestions:'
            ' ccp4.github.io/dimple ###' % __version__)
    mtz_meta = wf.read_mtz_metadata(opt.mtz)
    _comment_summary_line('MTZ (%.1fA)' % mtz_meta.dmax, mtz_meta)
    if opt.dls_naming:
        opt.pdbs = dls_name_filter(opt.pdbs)
    opt.pdbs = utils.filter_out_duplicate_files(opt.pdbs, relto=opt.output_dir)
    if not opt.pdbs:
        comment('\nNo non-empty pdb files given. Nothing to do.')
        return
    for p in opt.pdbs:
        wf.read_pdb_metadata(p, print_errors=(len(opt.pdbs) > 1))
    if len(opt.pdbs) > 1:
        comment('\nPDBs in order of similarity (using the first one):')
        opt.pdbs.sort(key=lambda x: calculate_difference(wf.file_info[x],
                                                         mtz_meta))
    utils.log_value('data_file', opt.mtz)
    utils.log_value('pdb_files', opt.pdbs)
    for p in opt.pdbs:
        _comment_summary_line(os.path.basename(p), wf.file_info[p])
    ini_pdb = 'ini.pdb'
    wf.copy_uncompressed(opt.pdbs[0], ini_pdb)
    pdb_meta = wf.file_info[opt.pdbs[0]]
    if pdb_meta is None:
        put_error('PDB file missing CRYST1 record, starting from MR')
    if opt.no_hetatm or check_hetatm_x(wf.path(ini_pdb), pdb_meta):
        if not opt.no_hetatm:
            comment('\nHETATM marked as element X would choke many programs.')
        rb_xyzin = 'prepared.pdb'
        wf.temporary_files.add(rb_xyzin)
        n_het = wf.remove_hetatm(xyzin=ini_pdb, xyzout=rb_xyzin,
                                 remove_all=opt.no_hetatm)
        comment('\nRemoved %d HETATM atoms' % n_het)
    else:
        rb_xyzin = ini_pdb
    # run rwcontents even without CRYST1 - it will show mol. weight only
    wf.rwcontents(xyzin=rb_xyzin).run()
    rw_data = wf.jobs[-1].data
    if pdb_meta is None:
        pass  # we already had a warning message
    elif rw_data.get('solvent_percent') is None:
        put_error('rwcontents could not interpret %s' % rb_xyzin)
    elif rw_data['solvent_percent'] > HIGH_SOLVENT_PCT:
        comment('\nHmm... %.1f%% of solvent or incomplete model' %
                rw_data['solvent_percent'])
        if abs(wf.jobs[-1].data.get('volume', 0) - pdb_meta.get_volume()) > 10:
            comment('\ndebug: problem when calculating volume?')

    ####### pointless - reindexing #######
    if match_symmetry(mtz_meta, pdb_meta) and opt.mr_when_r > 0 and (
            0.7 < mtz_meta.get_volume() / pdb_meta.get_volume() < 1.4):
        reindexed_mtz = 'reindexed.mtz'
        # Reindexed file can be useful for further refinement, don't delete it.
        #wf.temporary_files.add(reindexed_mtz)
        wf.pointless(hklin=opt.mtz, xyzin=rb_xyzin, hklout=reindexed_mtz,
                     keys='TOLERANCE 5').run(may_fail=True)
        alt_reindex = wf.jobs[-1].data.get('alt_reindex')
        if wf.jobs[-1].exit_status == 0 and alt_reindex:
            for ar in alt_reindex:
                comment('\n    %-10s CC: %-8.3f cell diff: %.1fA' % (
                        ar['op'], ar['cc'], ar['cell_deviat']))
        else:
            # until recently (2015) pointless didn't print CC for non-ambiguous
            # spacegroups (e.g. C2), but now it always prints
            comment('\n    no good indexing')
            reindexed_mtz = opt.mtz
    else:
        reindexed_mtz = opt.mtz
    reindexed_mtz_meta = wf.read_mtz_metadata(reindexed_mtz)
    if reindexed_mtz_meta.symmetry != mtz_meta.symmetry:
        _comment_summary_line('reindexed MTZ', reindexed_mtz_meta)

    ####### (c)truncate - calculate amplitudes if needed #######
    if not opt.fcolumn:
        opt.fcolumn = 'F' if 'F' in mtz_meta.columns else 'FP'
    elif opt.icolumn or opt.ItoF_prog:
        put_error('Ignoring options --fcolumn/--sigfcolumn')
    opt.sigfcolumn = opt.sigfcolumn.replace('<FCOL>', opt.fcolumn)
    if (opt.ItoF_prog or opt.icolumn or opt.fcolumn not in mtz_meta.columns or
            opt.sigfcolumn not in mtz_meta.columns):
        f_mtz = 'amplit.mtz'
        wf.temporary_files.add(f_mtz)
        i_sigi_cols = _find_i_sigi_columns(mtz_meta, opt)
        if opt.ItoF_prog == 'ctruncate' or (opt.ItoF_prog is None and opt.slow):
            colano = None
            if opt.anode and all(col in mtz_meta.columns for col in
                                 ['I(+)', 'SIGI(+)', 'I(-)', 'SIGI(-)']):
                colano = '/*/*/[I(+),SIGI(+),I(-),SIGI(-)]'
            wf.ctruncate(hklin=reindexed_mtz, hklout=f_mtz,
                         colin='/*/*/[%s,%s]' % i_sigi_cols,
                         colano=colano).run()
        else:
            wf.truncate(hklin=reindexed_mtz, hklout=f_mtz,
                        labin='IMEAN=%s SIGIMEAN=%s' % i_sigi_cols,
                        labout='F=F SIGF=SIGF').run()
        opt.fcolumn = 'F'
        opt.sigfcolumn = 'SIGF'
    else:
        f_mtz = reindexed_mtz

    ####### rigid body - check if model is good for refinement? #######
    refmac_labin_nofree = 'FP=%s SIGFP=%s' % (opt.fcolumn, opt.sigfcolumn)
    refmac_xyzin = None
    cell_diff = calculate_difference(pdb_meta, reindexed_mtz_meta)
    if pdb_meta is None:
        pass  # the error message was already printed
    elif opt.mr_when_r <= 0:
        comment('\nMR requested unconditionally.')
    elif cell_diff > 0.1 and opt.mr_when_r < 1:
        comment('\nDifferent unit cells.')
    elif pdb_meta.symmetry != reindexed_mtz_meta.symmetry:
        comment('\nDifferent space groups.')
    else:
        comment('\nRigid-body refinement with resolution 3.5 A, %d cycles.' %
                opt.rigid_cycles)
        if 'aa_count' in rw_data and 'water_count' in rw_data:
            if rw_data['aa_count'] != 0:
                comment(' %.1f waters/aa.' % (rw_data['water_count'] /
                                              rw_data['aa_count']))
            else:
                comment(' %d/0 waters/aa.' % rw_data['water_count'])
        wf.temporary_files |= {'refmacRB.pdb', 'refmacRB.mtz'}
        # it may fail because of "Disagreement between mtz and pdb"
        wf.refmac5(hklin=f_mtz, xyzin=rb_xyzin,
                   hklout='refmacRB.mtz', xyzout='refmacRB.pdb',
                   labin=refmac_labin_nofree,
                   libin=None,
                   keys="""refinement type rigidbody resolution 15 3.5
                           rigidbody ncycle %d""" % opt.rigid_cycles
                   ).run(may_fail=True)
        # if the error is caused by mtz/pdb disagreement, continue with MR
        if wf.jobs[-1].exit_status != 0:
            comment('\nTry MR.')
        elif not wf.jobs[-1].data.get('overall_r'):
            comment('\nWARNING: unknown R factor, something went wrong.\n')
            refmac_xyzin = 'refmacRB.pdb'
        elif wf.jobs[-1].data['overall_r'] > opt.mr_when_r:
            comment('\nRun MR for R > %g.' % opt.mr_when_r)
        else:
            comment('\nNo MR for R < %g.' % opt.mr_when_r)
            refmac_xyzin = 'refmacRB.pdb'

    ####### phaser/molrep - molecular replacement #######
    if refmac_xyzin is None:
        vol_ratio = None
        if pdb_meta:
            # num_mol accounts for strict NCS (MTRIX without iGiven)
            vol_ratio = (mtz_meta.asu_volume() /
                         pdb_meta.asu_volume(rw_data['num_mol']))
            comment(' Volume of asu: %.1f%% of model asu.' % (100 * vol_ratio))
        if opt.mr_when_r >= 1:
            comment('\nWould try MR, but it is disabled.')
            return
        if opt.mr_num:
            mr_num = opt.mr_num
        else:
            mr_num = guess_number_of_molecules(mtz_meta, rw_data, vol_ratio)
        mw = rw_data.get('weight')
        if isinstance(mr_num, float):
            wf.ensembler(pdbin=rb_xyzin, root='ens').run()
            n_models = len(wf.jobs[-1].data['models'])
            mw = None
            rb_xyzin = 'ens_merged.pdb'
            mr_num = max(int(round(mr_num * n_models)), 1)
        # phaser is used by default if number of searched molecules is known
        if opt.mr_prog == 'molrep':
            wf.temporary_files |= {'molrep.pdb', 'molrep_dimer.pdb',
                                   'molrep.crd'}
            wf.molrep(f=f_mtz, m=rb_xyzin).run()
            refmac_xyzin = 'molrep.pdb'
        else:
            wf.temporary_files |= {'phaser.1.pdb', 'phaser.1.mtz'}
            wf.phaser_auto(hklin=f_mtz,
                           labin='F=%s SIGF=%s' % (opt.fcolumn, opt.sigfcolumn),
                           model=dict(pdb=rb_xyzin, identity=100, num=mr_num,
                                      mw=mw),
                           sg_alt='ALL', opt=opt,
                           root='phaser').run(may_fail=True)
            if not _after_phaser_comments(wf.jobs[-1],
                                          sg_in=reindexed_mtz_meta.symmetry):
                raise RuntimeError('No phaser solution.')
            refmac_xyzin = 'phaser.1.pdb'
            f_mtz = 'phaser.1.mtz'

    if False:
        wf.findwaters(pdbin=refmac_xyzin, hklin=f_mtz,
                      f='FC', phi='PHIC', pdbout='prepared_wat.pdb', sigma=2)
        refmac_xyzin = 'prepared_wat.pdb'

    ####### adding free-R flags #######
    f_mtz_meta = wf.read_mtz_metadata(f_mtz)
    cad_reso = opt.reso or (f_mtz_meta.dmax - MtzMeta.d_eps)
    if opt.free_r_flags:
        free_mtz = opt.free_r_flags
        free_col = check_freerflags_column(wf.path(free_mtz),
                                           expected_symmetry=pdb_meta,
                                           column=opt.freecolumn)
        comment('\nFree-R flags from the %s file, column %s.' %
                (('reference' if free_mtz != opt.mtz else 'input'), free_col))
    else:
        free_col = DEFAULT_FREE_COLS[0]
        if free_col in f_mtz_meta.columns:
            comment('\nReplace free-R flags')
        else:
            comment('\nGenerate free-R flags')
        free_mtz = 'free.mtz'
        wf.temporary_files |= {'unique.mtz', free_mtz}
        if opt.seed_freerflag or cell_diff > 1e3:  # i.e. different SG
            wf.unique(hklout='unique.mtz', ref=f_mtz_meta,
                      resolution=cad_reso).run()
        else:
            comment(' (repeatably)')
            # Here we'd like to have always the same set of free-r flags
            # for given PDB file. That's why we don't use information
            # from the data file (mtz).
            wf.unique(hklout='unique.mtz', ref=pdb_meta, resolution=1.0).run()
        # CCP4 freerflag uses always the same pseudo-random sequence by default
        wf.freerflag(hklin='unique.mtz', hklout=free_mtz,
                     keys=('SEED' if opt.seed_freerflag else '')).run()

    if free_mtz == opt.mtz and opt.reso is None:
        prepared_mtz = f_mtz
    else:
        prepared_mtz = 'prepared.mtz'
        wf.temporary_files.add(prepared_mtz)
        wf.cad(data_in=[(f_mtz,
                         [c for c in f_mtz_meta.columns if c != free_col]),
                        (free_mtz, [free_col])],
               hklout=prepared_mtz,
               keys=['sysab_keep',  # does it matter?
                     'reso overall 1000.0 %g' % cad_reso]).run()
    freerflag_missing = wf.count_mtz_missing(prepared_mtz, free_col)
    if freerflag_missing:
        wf.freerflag(hklin=prepared_mtz, hklout='prepared2.mtz',
                     keys='COMPLETE FREE='+free_col,
                     parser=' (again, for %d refl. more)' % freerflag_missing
                     ).run()
        prepared_mtz = 'prepared2.mtz'
        wf.temporary_files.add(prepared_mtz)

    ####### refinement #######
    if opt.weight:
        refmac_weight = 'matrix %f' % opt.weight
    else:
        refmac_weight = 'auto'
    restr_ref_keys = """\
     make newligand continue
     refinement type restrained
     weight %s
     """ % refmac_weight
    if opt.freecolumn_val:
        restr_ref_keys += 'free %s\n' % opt.freecolumn_val
    refmac_labin = '%s FREE=%s' % (refmac_labin_nofree, free_col)
    comment('\nRestrained refinement, %d+%d cycles.' % (opt.jelly,
                                                        opt.restr_cycles))
    if opt.jelly:
        wf.temporary_files |= {'jelly.pdb', 'jelly.mtz'}
        wf.refmac5(hklin=prepared_mtz, xyzin=refmac_xyzin,
                   hklout='jelly.mtz', xyzout='jelly.pdb',
                   labin=refmac_labin, libin=opt.libin,
                   keys=restr_ref_keys + 'ridge distance sigma 0.01\n'
                                         'make hydrogen no\n'
                                         'ncycle %d' % opt.jelly +
                                         opt.extra_ref_keys).run()
        comment(_refmac_rms_line(wf.jobs[-1].data))
        refmac_xyzin = 'jelly.pdb'
    restr_job = wf.refmac5(hklin=prepared_mtz, xyzin=refmac_xyzin,
                           hklout=opt.hklout, xyzout=opt.xyzout,
                           labin=refmac_labin, libin=opt.libin,
                           keys=(restr_ref_keys +
                                 'ncycle %d' % opt.restr_cycles +
                                 opt.extra_ref_keys)).run()
    comment(_refmac_rms_line(restr_job.data))
    # if that run is repeated with --from-step it's useful to compare Rfree
    if wf.from_job > 0 and wf.from_job <= len(wf.jobs):  # from_job is 1-based
        prev = [j for j in wf.repl_jobs if j.name == restr_job.name]
        if prev and prev[0].data and 'free_r' in prev[0].data:
            comment('\nPreviously:  R/Rfree %.4f/%.4f  Rfree change: %+.4f' % (
                    prev[0].data['overall_r'], prev[0].data['free_r'],
                    restr_job.data['free_r'] - prev[0].data['free_r']))

    ####### check blobs #######
    if opt.blob_search:
        if restr_job.data['free_r'] <= BAD_FINAL_RFREE:
            if opt.gemmi_blobs:
                fb_job = wf.gemmi_blobs(opt.hklout, opt.xyzout, sigma=0.8).run()
            else:
                fb_job = wf.find_blobs(opt.hklout, opt.xyzout, sigma=0.8).run()
            coot_script = _generate_scripts_and_pictures(wf, opt, fb_job.data)
            if coot_script:
                comment('\nTo see it in Coot run %s' % coot_script)
        else:
            comment('\nNo blob search for Rfree > %g.' % BAD_FINAL_RFREE)
            _generate_scripts_and_pictures(wf, opt, None)

    if opt.anode:
        # check if mtz contains I+/- and SIGI+/-
        column_types = list(reindexed_mtz_meta.columns.values())
        if column_types.count('K') != 2 and column_types.count('M') != 2:
            comment('\nColumns I+/- and SIG+/- not found. Skipping AnoDe.')
            return
        anode_name = 'anode'
        # convert to sca for input to shelxc
        scaout = anode_name + '.sca'
        wf.mtz2sca(prepared_mtz, scaout).run()

        wf.shelxc(scaout, reindexed_mtz_meta.cell,
                  reindexed_mtz_meta.symmetry).run()

        wf.copy_uncompressed(opt.xyzout, anode_name + '.pdb')
        anode_job = wf.anode(anode_name).run()
        wf.temporary_files |= {scaout, anode_name + '.pdb', anode_name + '.hkl',
                               anode_name + '_sad.cif', anode_name + '_fa.hkl'}
        cell = Cell(reindexed_mtz_meta.cell, reindexed_mtz_meta.symmetry)
        # need orthogonal not fractional coordinates to generate coot script
        anode_job.data['blobs'] = cell.orthogonalize(anode_job.data['xyz'])
        comment(_anode_anom_peak_lines(anode_job.data))
        coot_script = _generate_scripts_and_pictures(
            wf, opt, anode_job.data, pha=anode_name+'.pha')

def _find_i_sigi_columns(mtz_meta, opt):
    if opt.icolumn:
        icolumn = opt.icolumn
        mtz_meta.check_col_type(icolumn, 'J')
    else:
        j_columns = [k for k, v in mtz_meta.columns.items() if v == 'J']
        if len(j_columns) == 1:
            icolumn = j_columns[0]
        elif 'IMEAN' in j_columns:
            icolumn = 'IMEAN'
        elif len(j_columns) > 1:
            put_error('Multiple intensity columns: %s. '
                      'Pick one with  --icolumn' % j_columns)
            sys.exit(1)
        else:
            put_error('No intensity (IMEAN) column in the MTZ file')
            sys.exit(1)

    # the default value of sigicolumn ('SIG<ICOL>') needs substitution
    sigicolumn = opt.sigicolumn.replace('<ICOL>', icolumn)
    mtz_meta.check_col_type(sigicolumn, 'Q')
    return (icolumn, sigicolumn)

def _refmac_rms_line(data):
    rb, ra, rc = [data.get(k, (-1,))
                  for k in ('rmsBOND', 'rmsANGL', 'rmsCHIRAL')]
    return ('\n    RMS:   bond %.3f -> %.3f' % (rb[0], rb[-1]) +
            '   angle %.2f -> %.2f' % (ra[0], ra[-1]) +
            '   chiral %.2f -> %.2f' % (rc[0], rc[-1]))

def _anode_anom_peak_lines(data):
    out = ''
    for n, height in enumerate(data['height'][:6]):
        out += '\n' if n % 3 == 0 else '  '
        out += 'h=%.0f ' % height
        dist = data['distance'][n]
        out += ('at ' if dist < 0.5 else '%.0fA from ' % dist) + data['atom'][n]
    return out

def _after_phaser_comments(phaser_job, sg_in):
    phaser_data = phaser_job.data
    if 'error' in phaser_data:
        comment('\n' + phaser_data['error'])
    if (phaser_job.exit_status != 0 or
            phaser_data['info'] == 'Sorry - No solution'):
        comment('\nGiving up.')
        return False
    solu_set = phaser_data.get('status', '')
    if phaser_data['info'].endswith('...'):
        comment('\n...' + solu_set[len(phaser_data['info'])-3:])
    if phaser_data.get('partial_solution'):
        # counting TF*0 or TFZ=number, but not TFZ==number
        n_comp = (solu_set.count('TF') - solu_set.count('TFZ==') +
                  solu_set.count('+TNCS'))
        comment('\nSolution found with %d components.' % n_comp)
    if phaser_data['SG'] != sg_in:
        comment('\nSpacegroup changed to %s' % phaser_data['SG'])
    return True

def _comment_summary_line(name, meta):
    comment('\n%-21s %s' % (name, meta or '???'))

# returns either number of molecules (int) or a fraction of the molecule (float)
def guess_number_of_molecules(mtz_meta, rw_data, vol_ratio):
    # if the number of molecules seems to be 1 or 2, don't go into Matthews
    if vol_ratio and rw_data.get('solvent_percent', 100) < HIGH_SOLVENT_PCT:
        if 0.7 < vol_ratio < 1.33:
            return 1
        if 1.8 < vol_ratio < 2.2:
            return 2

    Va = mtz_meta.asu_volume()
    m = rw_data['weight']

    # Vm = Va/(n*M)
    # Vs = 1 - 1.23/Vm  => Vs = 1 - n * 1.23*M/Va
    def calc_Vs(nmol):
        return 100 * (1 - nmol * 1.23 * m / Va)

    # For our purpose, it's better to overestimate the number of molecules,
    # because we can use "partial solution" from Phaser.
    # OTOH the search with overestimated n is slower and more likely to fail.
    # We also have preference for even numbers because they are more frequent
    # and Phaser can make use of tNCS if it's present.
    # Let's pick the largest n that gives solvent content (Vs) at least 30%.
    # If n is odd, try n-1 if Vs is still above 45%.

    # Vm = Va/(n*M)  =>  n = Va/(Vm*M)
    # 1-1.23/Vm=30% => Vm=1.76
    n = max(int(Va / (1.76 * m)), 1)
    if n % 2 == 1 and calc_Vs(n-1) < 45:
        n -= 1

    Vsn = calc_Vs(n)
    if n > 1:
        # 1-1.23/Vm=50% => Vm=2.46
        other_n = min(int(round(Va / (2.46 * m))), n-1)
        comment('\n%.0f%% solvent for %d, %.0f%% for %d components.'
                % (calc_Vs(other_n), other_n, Vsn, n))
    elif Vsn > 0:
        comment('\n%.0f%% solvent for single component.' % Vsn)
    else:
        comment('\nModel too big to fit in the unit cell.')

    # if model is too big we will try to split it
    if Vsn < VERY_LOW_SOLVENT_PCT or (Vsn < LOW_SOLVENT_PCT and vol_ratio):
        comment(' Let us try to split the model.')
        return float(vol_ratio or Va / (2.4 * m))
    return n


def _write_script(path, content, executable=False):
    try:
        with open(path, 'w') as f:
            f.write(content)
        if executable:  # chmod +x
            mode = os.stat(path).st_mode
            os.chmod(path, mode | ((mode & 0o444) >> 2))
    except (IOError, OSError) as e:
        put_error(e)

def _generate_scripts_and_pictures(wf, opt, data, pha=None):
    blobs = data['blobs'] if data else []
    coot_path = coots.find_path()
    if not blobs:
        comment('\nUnmodelled blobs not found.')
    elif opt.img_format:
        if coot_path:
            coot_ver = coots.find_version(coot_path)
            if coot_ver is None:
                put_error('coot not working(?), no pictures')
                opt.img_format = None
            elif 'with python' not in coot_ver:
                put_error('coot with Python support is needed')
                opt.img_format = None
        else:
            put_error('No coot, no pictures')
            opt.img_format = None
        if not utils.syspath('render'):
            put_error('No Raster3d, no pictures')
            opt.img_format = None
        if opt.img_format:
            if len(blobs) == 1:
                comment('\nRendering density blob at (%.1f, %.1f, %.1f)' %
                        blobs[0])
            else:
                comment('\nRendering 2 largest blobs: at (%.1f, %.1f, %.1f) '
                        'and at (%.1f, %.1f, %.1f)' % (blobs[0]+blobs[1]))
    com = data and data.get('center')
    if pha:
        normal_map = False
        refl = pha
        prefix = 'anom-'
    else:
        normal_map = True
        refl = opt.hklout
        prefix = ''

    # run-coot.py centers on the biggest blob. It uses relative paths -
    # it can be run only from the output directory, but is not affected
    # by moving that directory to different location.
    # There are blobN-coot.py scripts generated below with absolute paths.
    # write coot script (apart from pictures) that centers on the biggest blob
    script_path = os.path.join(wf.output_dir, prefix + 'run-coot.py')
    script = coots.basic_script(pdb=opt.xyzout, refl=refl,
                                normal_map=normal_map,
                                center=(blobs and blobs[0]), toward=com,
                                white_bg=opt.white_bg)
    _write_script(script_path, script, executable=True)

    # blob images, for now for not more than two blobs
    d = os.path.abspath(wf.output_dir)
    for n, b in enumerate(blobs[:2]):
        py_path = os.path.join(wf.output_dir,
                               '%sblob%d-coot.py' % (prefix, n+1))
        content = coots.basic_script(pdb=os.path.join(d, opt.xyzout),
                                     refl=os.path.join(d, refl),
                                     normal_map=normal_map,
                                     center=blobs[n], toward=com,
                                     white_bg=opt.white_bg)
        _write_script(py_path, content)
    # coot.sh - one-line script for convenience
    if blobs:
        coot_sh_text = '{coot} --no-guano {out}/%sblob1-coot.py\n' % prefix
    else:
        coot_sh_text = '{coot} --no-guano {out}/final.mtz {out}/final.pdb\n'
    coot_sh_path = os.path.join(wf.output_dir, prefix + 'coot.sh')
    _write_script(coot_sh_path, coot_sh_text.format(coot=coot_path or 'coot',
                                                    out=wf.output_dir),
                  executable=True)

    if opt.img_format and blobs:
        script = ''
        basenames = []
        # as a workaround for buggy coot the maps are reloaded for each blob
        for n, b in enumerate(blobs[:2]):
            script += coots.basic_script(pdb=opt.xyzout, refl=refl,
                                         normal_map=normal_map,
                                         center=b, toward=com,
                                         white_bg=opt.white_bg)
            rs, names = coots.r3d_script(center=b, toward=com,
                                         blobname='%sblob%s' % (prefix, n+1))
            script += rs
            basenames += names
        coot_job = wf.coot_py(script)
        try:
            coot_job.run()
        except workflow.JobError:
            # check for a possible cause to hint the user
            # (possible workaround: change $HOME to non-existing directory)
            if utils.silently_run(coot_job.args, cwd=wf.output_dir)[0] != 0:
                put_error('coot fails with options: --no-graphics --python',
                          comment='It happens when scripts in .coot or '
                                  '.coot-preferences are not compatible\n'
                                  'with the --no-graphics mode.')
            raise
        for n, basename in enumerate(basenames):
            try:
                job = wf.render_r3d(basename, img_format=opt.img_format)
                if n % 3 == 0:
                    job.run()
                else:  # minimal output
                    job.run(show_progress=False, new_line=False)
                wf.delete_files([basename + '.r3d'])
            except workflow.JobError as e:
                # Raster3D may fail saying "increase MAXDET and recompile".
                # This is not critical, so Dimple doesn't stop.
                put_error('Rendering failed, no picture', comment=' ' + e.note)
    return coot_sh_path


def parse_dimple_commands(args):
    dstr = ' (default: %(default)s)'
    parser = argparse.ArgumentParser(  # noqa: E126 visual indent
                usage=USAGE_SHORT, epilog=workflow.commands_help, prog=PROG,
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
    group2.add_argument('--sigicolumn', metavar='SIGICOL', default='SIG<ICOL>',
                        help='SIGI column label'+dstr)
    group2.add_argument('--fcolumn', metavar='FCOL',
                        help='F column label (default: F)')
    group2.add_argument('--sigfcolumn', metavar='SIGFCOL', default='SIG<FCOL>',
                        help='SIGF column label'+dstr)
    group2.add_argument('--libin', metavar='CIF',
                        help='ligand descriptions for refmac (LIBIN)')
    group2.add_argument('--refmac-key', metavar='LINE', action='append',
                        help='extra Refmac keywords to be used in refinement')
    group2.add_argument('-R', '--free-r-flags', metavar='MTZ_FILE',
                        help='file with freeR flags '
                             '("-" = use flags from data mtz)')
    group2.add_argument('--freecolumn', metavar='COL[=N]',
                        help='Rfree column with optional value (default: 0)')
    group2.add_argument('--hklout', metavar='out.mtz', default='final.mtz',
                        help='output mtz file'+dstr)
    group2.add_argument('--xyzout', metavar='out.pdb', default='final.pdb',
                        help='output pdb file'+dstr)
    group2.add_argument('-f', choices=['png', 'jpeg', 'none'],
                        dest='img_format',
                        help='format of generated images'+dstr)
    group2.add_argument('--white-bg', dest='white_bg', action='store_true',
                        help='white background in Coot and in images')
    group2.add_argument('--no-cleanup', dest='cleanup', action='store_false',
                        help='leave intermediate files')
    group2.add_argument('--cleanup', action='store_true',
                        help=argparse.SUPPRESS)  # obsolete

    group_w = parser.add_argument_group('what is calculated')
    group_w.add_argument('--no-blob-search',
                         dest='blob_search', action='store_false',
                         help='do not search for unmodelled blobs')
    group_w.add_argument('--anode', action='store_true',
                         help='use SHELX/AnoDe to find peaks in anomalous map')

    group3 = parser.add_argument_group('options customizing the run')
    group3.add_argument('--no-hetatm', action='store_true',
                        help='remove HETATM atoms from the given model')
    group3.add_argument('--rigid-cycles', metavar='N', type=int,
                        help='cycles of rigid-body refinement (default: 10)')
    group3.add_argument('--jelly', metavar='N', type=int,
                        help='cycles of jelly-body refinement (default: 4)')
    group3.add_argument('--restr-cycles', metavar='N', type=int,
                        help='cycles of refmac final refinement (default: 8)')
    group3.add_argument('--reso', type=float, help='limit the resolution [A]')
    group3.add_argument('--weight', metavar='VALUE', type=float,
                        help='refmac matrix weight (default: auto-weight)')
    group3.add_argument('--mr-prog', choices=['phaser', 'molrep'],
                        default='phaser',
                        help='Molecular Replacement program' + dstr)
    group3.add_argument('--mr-num', type=int,
                        help='number of molecules for MR (default: auto)')
    group3.add_argument('--mr-reso', type=float, default=3.25,
                        help='high resolution for MR '
                             '(if >10 interpreted as eLLG)' + dstr)
    group3.add_argument('--ItoF-prog', choices=['truncate', 'ctruncate'],
                        help='program to calculate amplitudes')
    group3.add_argument('--gemmi-blobs', action='store_true',
                        help=argparse.SUPPRESS)
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
    legacy_args = {'HKLIN': '', 'XYZIN': '',
                   'HKLOUT': '--hklout', 'XYZOUT': '--xyzout'}
    if len(args) == 8 and args[0] in legacy_args:
        args = [legacy_args.get(a) or a
                for a in args if legacy_args.get(a) != '']
        output_dir = os.path.join(os.environ.get('CCP4_SCR', ''), 'dimple_out')
        args.append(output_dir)

    # special mode for checking pdb file[s]
    if len(args) >= 1 and all(arg.endswith('.pdb') for arg in args):
        special_pdb_mode(args)
        sys.exit(0)

    # special mode for checking mtz file
    if len(args) == 1 and args[0].endswith('.mtz'):
        special_mtz_mode(args)
        sys.exit(1)

    opt = parser.parse_args(args)
    all_args = [opt.pos_arg1, opt.pos_arg2, opt.pos_arg3] + opt.more_args
    # all_args should be one mtz, one or more pdbs and output_dir
    opt.output_dir = all_args.pop()
    if opt.img_format == 'none':  # this option is kept for compatibility only
        opt.img_format = None
    if (opt.output_dir.endswith('.mtz') or opt.output_dir.endswith('.pdb') or
            opt.output_dir.endswith('.gz')):
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
            if not os.path.exists(old_pdb):
                if not old_wf.get('pdb_files'):
                    put_error('No pdb files in the original run?')
                    sys.exit(1)
                old_pdb = os.path.join(old_dir, old_wf['pdb_files'][0])
            if 'data_file' not in old_wf:  # temporary, to be removed soon
                old_mtz_arg = [a for a in old_wf['args'].split()
                               if a.endswith('.mtz')][0]
                old_wf['data_file'] = os.path.join(old_wf['cwd'], old_mtz_arg)
            old_mtz = os.path.join(old_dir, old_wf['data_file'])
        except (TypeError, KeyError):
            put_error('Reading logfile failed', comment='is it dimple.log?')
            sys.exit(1)
        all_args[0:2] = [old_pdb, old_mtz]

    mtz_args = [a for a in all_args if a.lower().endswith('.mtz')]
    if len(mtz_args) != 1:
        put_error('One mtz file should be given')
        sys.exit(1)
    opt.mtz = mtz_args[0]
    all_args.remove(opt.mtz)
    opt.pdbs = all_args
    for n, a in enumerate(opt.pdbs):
        if is_pdb_id(a):
            opt.pdbs[n] = download_pdb(a, opt.output_dir)
        elif not any(a.lower().endswith(ext) for ext in ['.pdb', '.pdb.gz',
                                                         '.ent', '.ent.gz']):
            put_error('unexpected arg (neither mtz nor pdb): %s' % a)
            sys.exit(1)
    if len(opt.pdbs) == 0:
        put_error('At least one pdb file should be given')
        sys.exit(1)
    if opt.seed_freerflag and opt.free_r_flags:
        put_error('Option --seed-freerflag and --free-r-flags'
                  ' do not make sense together')
        sys.exit(1)
    if opt.free_r_flags == '-':
        opt.free_r_flags = opt.mtz
    opt.freecolumn_val = None
    if opt.freecolumn and '=' in opt.freecolumn:
        opt.freecolumn, opt.freecolumn_val = opt.freecolumn.rsplit('=', 1)
    if opt.freecolumn and not opt.free_r_flags:
        if opt.freecolumn == DEFAULT_FREE_COLS[0] and opt.freecolumn_val:
            pass  # this may be useful for excluding different set
        else:
            put_error('--freecolumn suggests that you want to use existing free'
                      ' flags.\nFor this you need also option --free-r-flags')
            sys.exit(1)
    opt.extra_ref_keys = ''.join('\n'+key for key in opt.refmac_key or [])

    # extra checks
    for filename in opt.pdbs + [opt.mtz, opt.free_r_flags, opt.libin]:
        if filename and not os.path.isfile(filename):
            put_error('File not found: ' + filename)
            sys.exit(1)
    if os.path.exists(opt.output_dir) and not os.path.isdir(opt.output_dir):
        put_error('Not a directory: ' + opt.output_dir)
        sys.exit(1)

    # Since we'll execute programs from opt.output_dir, adjust paths.
    opt.mtz = utils.adjust_path(opt.mtz, opt.output_dir)
    opt.pdbs = [utils.adjust_path(a, opt.output_dir) for a in opt.pdbs]
    if opt.free_r_flags:
        opt.free_r_flags = utils.adjust_path(opt.free_r_flags, opt.output_dir)
    if opt.libin:
        opt.libin = utils.adjust_path(opt.libin, opt.output_dir)

    # set defaults that depend on the 'slow' level
    if opt.slow is None:
        opt.slow = 0
    elif opt.slow > 2:
        opt.slow = 2
    if opt.rigid_cycles is None:
        opt.rigid_cycles = 10
    if opt.restr_cycles is None:
        opt.restr_cycles = [8, 10, 12][opt.slow]
    if opt.jelly is None:
        opt.jelly = [4, 10, 100][opt.slow]

    return opt

def special_pdb_mode(args):
    print('Proper usage: %s' % USAGE_SHORT)
    check_ccp4_envvars()
    print('...actually we can run rwcontents for you')
    wf = workflow.Workflow('')
    wf.enable_logs = False
    for p in args:
        try:
            wf.read_pdb_metadata(p, print_errors=True)
            _comment_summary_line(os.path.basename(p), wf.file_info[p])
            wf.rwcontents(xyzin=p).run()
        except (IOError, RuntimeError) as e:
            put_error(e)
        except workflow.JobError as e:
            put_error(e.msg, comment=e.note)
    print('\n\n...but this is NOT how dimple is supposed to be run.')

def special_mtz_mode(args):
    print('Usage: %s' % USAGE_SHORT)
    check_ccp4_envvars()
    wf = workflow.Workflow('')
    wf.enable_logs = False
    try:
        mtz_meta = wf.read_mtz_metadata(args[0])
        print('Basic MTZ file info:')
        print(mtz_meta.info())
        contam_info = contaminants.get_info(mtz_meta)
        if contam_info:
            print(contam_info)
    except (IOError, RuntimeError) as e:
        put_error(e)

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
        comment('\n%d of %d PDBs have filenames matching data directory'
                % (len(matched_pdbs), len(pdbs)))
    return matched_pdbs

def check_ccp4_envvars():
    for necessary_var in ('CCP4', 'CCP4_SCR'):
        if necessary_var not in os.environ:
            put_error('$%s not found, giving up' % necessary_var)
            sys.exit(1)
    if not os.path.isdir(os.environ['CCP4_SCR']):
        put_error('No such directory: $CCP4_SCR, refmac shall not work!')

def check_contaminants_if_bad(wf, mtz):
    ref_job = wf.get_final_refinement_job()
    if not ref_job or ref_job.data.get('free_r', 1) > GOOD_FINAL_RFREE:
        mtz_meta = wf.read_mtz_metadata(mtz)  # it's cached
        info = contaminants.get_info(mtz_meta)
        if info:
            comment('\n' + info)

def main(args):
    if workflow.parse_workflow_commands():
        return

    options = parse_dimple_commands(args)
    check_ccp4_envvars()
    try:
        wf = workflow.Workflow(options.output_dir, from_job=options.from_step)
        utils.start_log(os.path.join(options.output_dir, 'dimple.log'),
                        output_dir=options.output_dir)
        utils.log_value('version', __version__)
        utils.start_log_screen(os.path.join(options.output_dir, 'screen.log'))

        dimple(wf=wf, opt=options)
        check_contaminants_if_bad(wf, mtz=options.mtz)
        exit_status = 0
    except workflow.JobError as e:
        put_error(e.msg, comment=e.note)
        try:
            utils.report_disk_space([wf.output_dir, os.getenv('CCP4_SCR')])
        except KeyboardInterrupt:
            comment('\nok, exiting...')
        exit_status = 1
    except (RuntimeError, IOError, OSError) as e:
        put_error(e)
        exit_status = 1
    finally:
        comment('\n')
    if options.cleanup:
        wf.delete_files(wf.temporary_files)
    wf.options = options
    try:
        wf.dump_pickle()
    except IOError as e:
        put_error(e)
        exit_status = 1
    return exit_status

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
