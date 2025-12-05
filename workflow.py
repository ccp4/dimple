import os
import sys
from subprocess import Popen, PIPE
import gzip
import errno
import shlex
import re
import threading
import queue
import time
import pickle
import shutil
if __name__ == '__main__' and __package__ is None:
    sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__
                                                                       ))))
import gemmi
from dimple import coots
from dimple import mtz
from dimple import pdb
from dimple import utils

_jobindex_fmt = '%3d '
_jobname_fmt = '%-15s'
_elapsed_fmt = '%5.1fs  '

PICKLE_FILENAME = 'workflow.pickle'


# heh, python from MSYS2 is handy, but needs some monkey-patching
if sys.platform == 'msys':
    old_opjoin = os.path.join
    def new_opjoin(*args):
        for n in range(len(args)-1, -1, -1):
            if n == 0 or os.path.isabs(args[n]):
                return old_opjoin(*args[n:])
    os.path.join = new_opjoin


class JobError(Exception):
    def __init__(self, msg, note=None):  # pylint: disable=super-init-not-called
        self.msg = msg
        self.note = note


class Output:
    "storage for Job's stdout/stderr"
    def __init__(self, role):
        self.role = role  # 'out' or 'err'
        self.file_extension = {'out': 'log', 'err': 'err'}[role]
        self.lines = []
        self.saved_to = None
        self.que = None

    def __nonzero__(self):
        return bool(self.lines or self.saved_to)

    def size_as_str(self):
        if self.lines:
            return '%d lines' % len(self.lines)
        else:
            return '-      '

    def read_line(self):
        while self.que is not None:
            try:
                line = self.que.get_nowait()
            except queue.Empty:
                break
            self.lines.append(line)
            yield line

    def finish_que(self):
        if self.que:
            while not self.que.empty():
                self.lines.append(self.que.get_nowait())
            self.que = None

    def save_output(self, output_dir, basename, remove_long_list=True):
        filename = basename + '.' + self.file_extension
        if self.lines:
            with open(os.path.join(output_dir, filename), 'wb') as f:
                for line in self.lines:
                    f.write(line)
            self.saved_to = filename
            utils.log_value('std'+self.role, filename)
            if remove_long_list and len(self.lines) > 5:
                self.lines = []

    def summary(self):
        n = len(self.lines)
        if 0 < n <= 3:
            return b''.join(self.lines)
        elif self.saved_to:
            return '-> %s' % self.saved_to
        elif n == 0:
            return ''
        else:  # n > 3
            return b''.join(self.lines[:3]) + ('%d more lines' % (n-3)).encode()


class Job:
    def __init__(self, workflow, prog):
        self.name = os.path.basename(prog) or prog  # only used to show info
        self.workflow = workflow
        self.args = [prog]
        self.std_input = ''
        self.stdin_file = None  # if set, it overwrites std_input
        # the rest is set after the job is run
        self.exit_status = None
        self.out = Output('out')
        self.err = Output('err')
        self.started = None  # will be set to time.time() at start
        self.total_time = None  # will be set when job ends
        # possible values: None (stdout preview),
        #                  string starting with space (' ') that is just shown,
        #                  or name of global function that parses output
        self.parser = None
        # job-specific data from output parsing
        self.data = {}

    def __repr__(self):
        if self.started:
            t = time.strftime(' %Y-%m-%d %H:%M', time.localtime(self.started))
        else:
            t = ''
        return '<Job %s%s>' % (self.name, t)

    def args_as_str(self):
        s = ' '.join(shlex.quote(a) for a in self.args)
        if self.stdin_file:
            s += ' < ' + self.stdin_file
        elif self.std_input:
            s += ' << EOF\n%s\nEOF' % self.std_input
        return s

    def run(self, show_progress=True, new_line=True, may_fail=False):
        self.workflow.run_job(job=self,
                              show_progress=show_progress, new_line=new_line)
        # exit_status may be None if --from-step is used
        if not may_fail and self.exit_status:
            notes = [self.args_as_str(), '']
            if self.out.saved_to:
                notes += ['stdout -> %s/%s' % (self.workflow.output_dir,
                                               self.out.saved_to)]
            if self.err:
                notes += ['stderr:', str(self.err.summary())]
            raise JobError('%s failed (exit status %d)' % (self.name,
                                                           self.exit_status),
                           note='\n'.join(notes))
        return self

    def parse(self):
        if self.parser is None:  # preview mode
            # generic non-parser
            line = b''
            for line in self.out.read_line():
                pass
            if line:
                trimmed = _PREVIEW_DISCARD_RE.sub('', line.strip().decode())
                ret = '[%d] %-44.44s' % (len(self.out.lines), trimmed)
            else:
                ret = 'stdout:%11s' % self.out.size_as_str()
                if self.err:
                    ret += ' stderr: %s' % self.err.size_as_str()
                ret = ret.ljust(50)
            return ret

        elif self.parser == '' or self.parser[0] == ' ':
            return self.parser

        else:
            p = globals()[self.parser]
            return p(self)


def _format(fmt, arg):
    return (fmt % arg) if arg else ''

# parsers for various programs
def _find_blobs_parser(job):
    if 'blobs' not in job.data:
        job.data['blobs'] = []
        job.data['scores'] = []
    for line in job.out.read_line():
        if line.startswith(b'#'):
            sp = line.split(None, 6)
            score = float(sp[5])
            if True:  # score > 150: XXX: better scoring may be needed
                x, y, z = sp[-1].strip(b'() \t\r\n').split(b',')
                job.data['blobs'].append((float(x), float(y), float(z)))
                job.data['scores'].append(score)
        elif line.startswith(b'Protein mass center:'):
            sp = line.split(b'(')[1].rstrip(b'\r\n )').split(b',')
            ctr = tuple(float(x) for x in sp)
            job.data['center'] = ctr
        elif line.startswith(b'Density std.dev'):
            job.data['density_info'] = line.strip().decode()
    scores = job.data['scores']
    if scores:
        return 'Blob scores: ' + ' '.join('%.0f' % sc for sc in scores[:8])
    elif 'density_info' in job.data:
        return 'searching with' + job.data['density_info'].split(',')[1]
    else:
        return ''

def _gemmi_blobs_parser(job):
    if 'blobs' not in job.data:
        job.data['blobs'] = []
        job.data['scores'] = []
    for line in job.out.read_line():
        if line.startswith(b'#'):
            score = float(line.split()[1])
            x, y, z = line.split(b'(')[1].split(b')')[0].split(b',')
            job.data['blobs'].append((float(x), float(y), float(z)))
            job.data['scores'].append(score)
        elif line.startswith(b'Center of mass:'):
            ctr = tuple(float(x) for x in line.split()[3:])
            job.data['center'] = ctr
    scores = job.data['scores']
    if scores:
        return 'Blob scores: ' + ' '.join('%.1f' % sc for sc in scores[:8])
    else:
        return ''

def _anode_parser(job):
    if 'xyz' not in job.data:
        job.data['xyz'] = []
        job.data['height'] = []
        job.data['sof'] = []
        job.data['distance'] = []
        job.data['atom'] = []
    found_strongest_peaks = False
    for line in job.out.read_line():
        if b'Strongest unique anomalous peaks' in line:
            found_strongest_peaks = True
            continue
        if found_strongest_peaks:
            tokens = line.split()
            if len(tokens) == 8:
                job.data['xyz'].append(tuple(float(t) for t in tokens[1:4]))
                job.data['height'].append(float(tokens[4]))
                job.data['sof'].append(float(tokens[5]))
                job.data['distance'].append(float(tokens[6]))
                job.data['atom'].append(tokens[7].decode())
    if found_strongest_peaks:
        return ('%s anomalous peaks with height h>4 sigma'
                % len(job.data['height']))
    else:
        return ''

def _rwcontents_parser(job):
    d = job.data
    for line in job.out.read_line():
        if line.startswith(b' Cell volume:'):
            vol = float(line.split(b':')[-1])
            if vol != 0:
                d['volume'] = vol
        elif line.startswith(b' Molecular Weight of protein:'):
            d['weight'] = float(line.split(b':')[-1])
        elif line.startswith(b' Molecular Weight of all atoms:'):
            d['total_weight'] = float(line.split(b':')[-1])
        elif line.startswith(b' Number of amino-acids residues ='):
            d['aa_count'] = int(line.split(b'=')[-1])
        elif line.startswith(b'                      - number of waters'):
            d['water_count'] = float(line.split()[-1])
        elif line.startswith(b' The Matthews Coefficient is :'):
            Vm = float(line.split(b':')[-1])
            if Vm != 0:
                d['Vm'] = Vm
                # 1.23 is used in Rupp's papers and in Phaser
                d['solvent_percent'] = (1 - 1.23/Vm) * 100
    if 'volume' in d and 'weight' in d and 'Vm' in d and 'num_mol' not in d:
        d['num_mol'] = int(round(d['volume'] / d['weight'] / d['Vm']))
    protein_kDa = d.get('weight', 0) / 1000.  # Da -> kDa
    total_kDa = d.get('total_weight', 0) / 1000.
    msg = '%s x %.0fkDa (+ %.0fkDa het)' % (d.get('num_mol', '??'),
                                            protein_kDa,
                                            total_kDa - protein_kDa)
    if 'volume' in d:
        msg += ' in %.fnm3, %.0f%% solvent' % (d['volume'] / 1000,
                                               d.get('solvent_percent', 0))
    return msg

def _cad_parser(job):
    d = job.data
    for line in job.out.read_line():
        # for now we're only interested in number of reflections from HKLIN1
        if b'* Number of Reflections =' in line and b'refl_in1' not in d:
            d['refl_in1'] = int(line.split(b'=')[1])
        elif b' Final Total of Unique records to HKLOUT =' in line:
            d['refl_out'] = int(line.split(b'=')[1])
    return '#refl  %s -> %s' % (d.get('refl_in1', ''), d.get('refl_out', ''))


class Ccp4LogTable(object):
    def __init__(self, title_line):
        assert b'$TABLE:' in title_line
        self.title = title_line.split(b':')[1].decode()
        self.columns = []
        self.data = []
        self.section = 0  # 0=header, 1=columns, 2=data, 3=end

    def send_line(self, line):
        line = line.strip()
        if line.startswith(b'$$'):
            self.section += 1
            if len(line) > 2:
                self.send_line(line[2:])
            return self.section != 3

        if self.section == 1:
            self.columns += line.rstrip(b'$').decode().split()
        elif self.section == 2:
            self.data.append(line.decode().split())
        return True

    def column(self, name):
        try:
            idx = self.columns.index(name)
            return [float(a[idx]) for a in self.data]
        except (ValueError, IndexError):
            return


def _refmac_parser(job):
    if 'cycle' not in job.data:
        # ini_free_r, free_r and iter_free_r are set optionally
        job.data['cycle'] = 0
        # iter_*_r values were added in dimple 1.5
        job.data['iter_overall_r'] = []
    for line in job.out.read_line():
        if 'sink' in job.data:
            more = job.data['sink'].send_line(line)
            if not more:
                for name in ['rmsBOND', 'rmsANGL', 'rmsCHIRAL']:
                    col_data = job.data['sink'].column(name)
                    if col_data and any(x != 0 for x in col_data):
                        job.data[name] = col_data
                del job.data['sink']
        elif line.startswith(b'Free R factor'):
            job.data['free_r'] = float(line.split(b'=')[-1])
            if 'ini_free_r' not in job.data:
                job.data['ini_free_r'] = job.data['free_r']
                job.data['iter_free_r'] = []
            job.data['iter_free_r'].append(job.data['free_r'])
        elif line.startswith(b'Overall R factor'):
            job.data['overall_r'] = float(line.split(b'=')[-1])
            if 'ini_overall_r' not in job.data:
                job.data['ini_overall_r'] = job.data['overall_r']
            job.data['iter_overall_r'].append(job.data['overall_r'])
        elif (line.startswith(b'     Rigid body cycle =') or
              line.startswith(b'     CGMAT cycle number =')):
            job.data['cycle'] = int(line.split(b'=')[-1])
        elif line.startswith(b'$TABLE: Rfactor analysis, stats vs cycle'):
            job.data['sink'] = Ccp4LogTable(line)
    cycle_str = '%2d/%d' % (job.data['cycle'], job.data.get('ncyc', -1))
    if 'ini_overall_r' in job.data:
        if 'ini_free_r' in job.data:
            return '%s   R/Rfree  %.4f/%.4f  ->  %.4f/%.4f' % (
                   cycle_str,
                   job.data['ini_overall_r'], job.data['ini_free_r'],
                   job.data['overall_r'], job.data['free_r'])
        return '%s   R  %.4f  ->  %.4f' % (
               cycle_str, job.data['ini_overall_r'], job.data['overall_r'])
    return cycle_str

# example:
# noqa
#     Alternative reindexing     Lklhd     CC    R(E^2)   Number Cell_deviation
#    [-k,-l,h+k+l]           0.079    0.029    0.506     61253      2.99
# in pointless 1.10.22 Phil added numbers in the first column
# 1              [h,k,l]              0.499   ...
_POINTLESS_ALTREINDEX_RE = re.compile(br'^\s*\d*\s+(\[[hkl+, -]+\]'
                                      br'[ \t\r0-9+.eE-]+)$')

_PREVIEW_DISCARD_RE = re.compile(r'[^\w!"#$%&\'()*+,./:;<=>?@[\\]^_`{|}~ -]')

def _float_or_nan(s):
    try:
        x = float(s)
    except ValueError:
        x = float('nan')
    return x

def _pointless_parser(job):
    # the line with 'reflections copied' is the last we read
    # to avoid reading 'Alternative reindexing' duplicated in the summary
    if 'refl_out' not in job.data:
        for line in job.out.read_line():
            if line.startswith(b'Maximum resolution used:'):
                job.data['resol'] = float(line.split(b':')[1])
            elif line.startswith(b'Number of reflections:'):
                job.data['refl_ref'] = int(line.split(b':')[1])
            elif _POINTLESS_ALTREINDEX_RE.match(line):
                s = _POINTLESS_ALTREINDEX_RE.match(line).group(1).split()
                job.data.setdefault('alt_reindex', []).append(
                        {'op': s[0].decode(),  # noqa E126 - indentation
                         'cc': _float_or_nan(s[2]),
                         'cell_deviat': _float_or_nan(s[-1])})
            elif line.startswith(b'   Cell:') and 'output_cell' not in job.data:
                s = line.split()[1:]
                job.data['output_cell'] = tuple(float(i) for i in s)
            elif b'reflections copied to output file' in line:
                job.data['refl_out'] = int(line.split()[0])
                break
    resol_txt = _format('%.2f', job.data.get('resol'))
    refl_out = job.data.get('refl_out', '')
    return 'resol. %4s A   #refl: %5s' % (resol_txt, refl_out)

def _phaser_parser(job):
    d = job.data
    for line in job.out.read_line():
        line = line.decode()
        if line.startswith('*** Phaser Module:'):
            d['info'] = '[%s]' % line[19:70].strip().lower()
        elif 'written to PDB file:' in line:
            d['expect_solu'] = 1
        elif 'expect_solu' in d:
            if line.startswith('   SOLU SET '):
                d['status'] = line[12:].strip()
                if len(d['status']) > 52:
                    d['info'] = d['status'][:50].rsplit(' ', 1)[0] + '...'
                else:
                    d['info'] = d['status']
            elif 'status' in d and ' SOLU ' not in line:  # continuation
                d['status'] += ' ' + line.strip()
            elif line.startswith('   SOLU SPAC '):
                d['SG'] = line[13:].strip()
                del d['expect_solu']
        elif 'Sorry - No solution' in line:
            d['info'] = line.strip('* \t\r\n')
            if 'No solution with all components' in line:
                d['partial_solution'] = 'yes'
        elif ' ERROR:' in line:
            # the error about impossible content has two lines, let's reword it
            d['error'] = line.strip().replace('a protein/nucleic acid',
                                              'impossible content')
    return '%-48s' % d.get('info', '')

def _ensembler_parser(job):
    job.data.setdefault('models', [])
    for line in job.out.read_line():
        if b'Model ' in line:
            job.data['models'].append(line.split()[-1].strip(b"')").decode())
    return 'ensembled chains: ' + ' '.join(job.data['models'])

def _truncate_parser(job):
    for line in job.out.read_line():
        if line.startswith(b' Least squares straight line gives:'):
            b_str = line.partition(b':')[-1].strip(b' \tB=').split()[0]
            job.data['B-factor'] = float(b_str)
    return 'B=%4s' % _format('%.1f', job.data.get('B-factor'))

def _ctruncate_parser(job):
    d = job.data
    for line in job.out.read_line():
        if line.startswith(b'$GRAPHS: Wilson plot - estimated B factor ='):
            d['B-factor'] = float(line.partition(b'=')[-1].split()[0])
        elif line.startswith(b'Eigenvalue ratios:'):
            d['eigval-ratios'] = tuple(float(i) for i in line.split()[-3:])
        elif line.startswith(b'L statistic ='):
            d['L-test'] = float(line.partition(b'=')[-1].split()[0])
    return 'B=%4s   aniso %14s   L-test:%s' % (
           _format('%.1f', d.get('B-factor')),
           _format('%.2f:%.2f:%.2f', d.get('eigval-ratios')),
           _format('%.2f', d.get('L-test')))


def ccp4_job(workflow, prog, logical=None, ki='', parser=None, add_end=True):
    """Handle traditional convention for arguments of CCP4 programs.
    logical - dictionary with where keys are so-called logical names.
    ki (string or list of lines) - Keyworded Input to be passed through stdin.
    add_end - adds "end" as the last line of stdin

    Note: "colin" and "labin" mean the same (column label),
    but different programs use different keywords.
    """
    job = Job(workflow, utils.cbin(prog))
    if logical:
        for a in ['hklin', 'hklout', 'hklref', 'xyzin', 'xyzout', 'libin']:
            if logical.get(a):
                job.args += [a.upper(), logical[a]]
    lines = ki.splitlines() if isinstance(ki, str) else ki
    stripped = [a.strip() for a in lines if a and not a.isspace()]
    if add_end and not (stripped and stripped[-1].lower() == 'end'):
        stripped.append('end')
    job.std_input = '\n'.join(stripped)
    job.parser = parser
    return job


def _print_progress(job, event):
    while not event.wait(0.5):
        p = job.parse()
        if p is not None:
            text = (_elapsed_fmt % (time.time() - job.started)) + p
            utils.put_temporarily(text)
            utils.reset_color()


def _start_enqueue_thread(file_obj):
    def enqueue_lines(f, q):
        for line in iter(f.readline, b''):
            q.put(line)
        f.close()
    que = queue.Queue()
    thr = threading.Thread(target=enqueue_lines, args=(file_obj, que))
    thr.daemon = True
    thr.start()
    return thr, que

def _get_input_as_string(job):
    if job.stdin_file:
        path = job.workflow.path(job.stdin_file)
        try:
            return open(path, 'rb').read()
        except IOError:
            raise JobError('cannot read input from: %s' % job.stdin_file)
    else:
        return job.std_input.encode()

def _run_and_parse(process, job):
    try:
        # job.*.que can be used by parsers (via Output.read_line() or directly)
        out_t, job.out.que = _start_enqueue_thread(process.stdout)
        err_t, job.err.que = _start_enqueue_thread(process.stderr)
        try:
            job_input = _get_input_as_string(job)
            process.stdin.write(job_input)
        except IOError as e:
            utils.put('\nWarning: passing input to %s failed.\n' % job.name)
            if e.errno not in (errno.EPIPE, errno.EINVAL):
                raise
        process.stdin.close()
        out_t.join()
        err_t.join()
        process.wait()
        # nothing is written to the queues at this point
        # parse what's left in the queues
        job.parse()
    finally:
        # take care of what is left by the parser
        job.out.finish_que()
        job.err.finish_que()


class Workflow:
    def __init__(self, output_dir, from_job=0):
        self.output_dir = os.path.abspath(output_dir)
        self.jobs = []
        self.file_info = {}
        self.temporary_files = set()
        self.from_job = from_job  # skip jobs before from_job (for testing)
        if from_job >= 1:
            try:
                _pkl = self.load_pickle()
                self.repl_jobs = _pkl.jobs
                self.file_info = _pkl.file_info
            except:
                self.repl_jobs = None
        self.dry_run = False
        self.enable_logs = True
        self.argv = sys.argv
        if not os.path.isdir(self.output_dir):
            try:
                os.makedirs(self.output_dir)
            except OSError as e:
                utils.put_error(e)
                sys.exit(1)
        # this can seriously affect Refmac compiled with GFortran
        bad_var = os.getenv('GFORTRAN_UNBUFFERED_ALL')
        if bad_var and bad_var[0] not in ('0', 'n', 'N'):
            utils.put_error(  # noqa
                      '$GFORTRAN_UNBUFFERED_ALL may terribly slow down Refmac',
                      comment='It is unset internally in dimple.')
            del os.environ['GFORTRAN_UNBUFFERED_ALL']
        # avoid html-like crap in the output of CCP4 program
        os.environ['CCP_SUPPRESS_HTML'] = '1'

    def __str__(self):
        return 'Workflow with %d jobs @ %s' % (len(self.jobs), self.output_dir)

    def path(self, rel_path):
        return os.path.join(self.output_dir, rel_path)

    def dump_pickle(self):
        with open(self.path(PICKLE_FILENAME), 'wb') as f:
            pickle.dump(self, f, -1)

    def load_pickle(self):
        with open(self.path(PICKLE_FILENAME), 'rb') as f:
            return pickle.load(f)

    def run_job(self, job, show_progress, new_line=True):
        if not hasattr(sys.stdout, 'isatty') or not sys.stdout.isatty():
            show_progress = False
        self.jobs.append(job)
        job_num = len(self.jobs)
        if new_line:
            utils.put('\n' + _jobindex_fmt % job_num)
            utils.put_green(_jobname_fmt % job.name)
        else:
            utils.put(' / %d' % job_num)
        sys.stdout.flush()
        utils.log_section(job.name)

        job_idx = len(self.jobs) - 1
        if job_idx < self.from_job - 1:  # from_job is 1-based
            # unpickle or skip
            if self.repl_jobs and len(self.repl_jobs) > job_idx:
                old_job = self.repl_jobs[job_idx]
                if old_job.name == job.name:
                    job.data = old_job.data
                    job = old_job
                    utils.put('unpickled')
                    utils.log_value('not_run', 'unpickled')
                    self.jobs[-1] = job
                else:
                    utils.put('skipped (mismatch)')
                    utils.log_value('not_run', 'unpickled/mismatch')
            else:
                utils.put('skipped')
                utils.log_value('not_run', 'skipped')
            return

        job.started = time.time()
        utils.log_time('start_time', job.started)
        if job.stdin_file:
            utils.log_value('stdin', job.stdin_file)
        elif job.std_input:
            utils.log_value('input', job.std_input)
        utils.log_value('prog', job.args[0])
        utils.log_value('args', ' '.join(shlex.quote(a) for a in job.args[1:]))
        utils.log_flush()
        # job.args[0] = 'true'  # for debugging
        try:
            process = Popen(job.args, stdin=PIPE, stdout=PIPE, stderr=PIPE,
                            cwd=self.output_dir)
        except OSError as e:
            if e.errno == errno.ENOENT:
                raise JobError('Program not found: %s' % job.args[0])
            else:
                raise

        if self.dry_run:
            return

        if show_progress:
            event = threading.Event()
            progress_thread = threading.Thread(target=_print_progress,
                                               args=(job, event))
            progress_thread.daemon = True
            progress_thread.start()

        try:
            if job.parser is not None or show_progress:
                _run_and_parse(process, job)
            else:
                job_input = _get_input_as_string(job)
                out, err = process.communicate(input=job_input)
                job.out.lines = out.splitlines(True)
                job.err.lines = err.splitlines(True)
        except KeyboardInterrupt:
            raise JobError('KeyboardInterrupt while running %s' % job.name,
                           note=job.args_as_str())
        finally:
            if show_progress:
                event.set()
                progress_thread.join()
            end_time = time.time()
            job.total_time = end_time - job.started
            utils.log_time('end_time', end_time)
            job.exit_status = process.poll()
            if new_line:
                utils.put(_elapsed_fmt % job.total_time)
            parse_output = job.parse()
            utils.put('%s' % (parse_output or ''))
            if parse_output:
                utils.log_value('info', parse_output)
            if self.enable_logs:
                self._write_logs(job)
            for k, v in job.data.items():
                if k != 'info':
                    utils.log_value(k, v)
        if job.exit_status != 0:
            utils.log_value('exit_status', job.exit_status)

    def _write_logs(self, job):
        log_basename = '%02d-%s' % (len(self.jobs), job.name.replace(' ', '_'))
        for output in (job.out, job.err):
            output.save_output(self.output_dir, log_basename)

    def remove_hetatm(self, xyzin, xyzout, remove_all):
        with open(self.path(xyzout), 'w') as out:
            return pdb.remove_hetatm(self.path(xyzin), out, remove_all)

    def read_pdb_metadata(self, xyzin, print_errors):
        if xyzin not in self.file_info:
            self.file_info[xyzin] = pdb.read_metadata(self.path(xyzin),
                                                      print_errors)
        return self.file_info[xyzin]

    def read_mtz_metadata(self, hklin):
        if hklin not in self.file_info:
            self.file_info[hklin] = mtz.read_metadata(self.path(hklin))
        return self.file_info[hklin]

    def count_mtz_missing(self, hklin, col):
        key = (hklin, 'missing', col)
        if key not in self.file_info:
            self.file_info[key] = mtz.get_num_missing(self.path(hklin), col)
        return self.file_info[key]

    def molrep(self, f, m, keys=''):
        job = Job(self, utils.cbin('molrep'))
        job.args += ['-f', f, '-m', m]
        if keys:
            job.args.append('-i')
            job.std_input = keys.strip() + '\nend'
        return job

    def phaser_auto(self, hklin, labin, model, root, sg_alt, opt):
        lines = [  # noqa
          'MODE MR_AUTO',
          'SEARCH METHOD FAST',
          'SEARCH DEEP OFF',
          'ENSEMBLE p PDBFILE "%(pdb)s" IDENTITY %(identity)g' % model,
          # if --no-hetatm was used HETATM records are already removed
          'ENSEMBLE p HETATOM ON',
          'SEARCH ENSEMBLE p NUM %(num)d' % model,
          'HKLIN "%s"' % hklin,
          'LABIN %s' % labin,
          'SGALTERNATIVE SELECT %s' % sg_alt,
          'ROOT %s' % root,
          # Since Phaser 2.5.6 template matched solutions are moved
          # to the template solution origin. Which is better than
          # getting a solution one cell away, so we set template here.
          'SOLUTION TEMPLATE original_model',
          'SOLUTION 6DIM ENSE p EULER 0 0 0 FRAC 0 0 0',
          # Final refinement with too high resolution crashes
          # with std::bad_alloc, even with 8GB of memory.
          'RESOLUTION AUTO HIGH 2.0',  # ignored by Phaser if > RESO HIGH
          ]
        if model['mw']:
            lines += ['COMPOSITION PROTEIN MW %(mw)f NUMBER %(num)d' % model]
        else:
            lines += ['COMPOSITION BY AVERAGE']
        if opt.slow < 2:
            lines += ['KILL TIME 180',  # 3h is much more than we want
                      # 'MACANO PROTOCOL OFF',
                      'PURGE ROT NUM 7',
                      'PURGE TRA NUM 20',
                      'PURGE RNP NUM 10']
        if opt.mr_reso > 10:  # in this case it's not a reso
            lines += ['ELLG TARGET %g' % opt.mr_reso]
        else:
            lines += ['RESOLUTION HIGH %g' % opt.mr_reso]
        # tNCS: we go with what phaser does by default -- tNCS of order 2
        # are handled automatically. While we could specify tNCS for
        # pseudo-tripling/quadrupling of the cell (TNCS NMOL 3) I don't know
        # if it'd do more good or bad.
        job = ccp4_job(self, 'phaser', ki=lines, parser='_phaser_parser')
        return job

    def ensembler(self, pdbin, root):
        job = Job(self, 'phaser.ensembler')
        job.name = 'ensembler'
        job.args += ['root=%s' % root, pdbin]
        job.parser = '_ensembler_parser'
        return job

    # functions below use logical=locals()
    # pylint: disable=unused-argument

    def pointless(self, hklin, xyzin, hklref=None, hklout=None, keys=''):
        return ccp4_job(self, 'pointless', logical=locals(), ki=keys,
                        parser='_pointless_parser')

    def unique(self, hklout, ref, resolution,
               labout='F=F_UNIQUE SIGF=SIGF_UNIQUE'):
        # Include also reflections that may be present in other spacegroups
        # belonging to the same pointgroup and with the same reflections.
        # (ignore systematic absences from screw axes by removing screw axes.)
        return ccp4_job(self, 'unique', logical=locals(),
                        ki=['cell %g %g %g %g %g %g' % tuple(ref.cell),
                            "symmetry '%s'" % ref.unscrewed_symmetry(),
                            'resolution %.3f' % resolution,
                            'labout %s' % labout],
                        parser='')

    def freerflag(self, hklin, hklout, keys='', parser=''):
        return ccp4_job(self, 'freerflag', logical=locals(), ki=keys,
                        parser=parser)

    #def reindex(self, hklin, hklout, symmetry):
    #    return ccp4_job(self, 'reindex', logical=locals(),
    #                    ki=["symmetry '%s'" % symmetry,
    #                        'reindex h,k,l'])

    def truncate(self, hklin, hklout, labin, labout):
        return ccp4_job(self, 'truncate', logical=locals(),
                        ki=['labin %s' % labin, 'labout %s' % labout,
                            'NOHARVEST'],
                        parser='_truncate_parser')

    def ctruncate(self, hklin, hklout, colin, colano):
        job = Job(self, 'ctruncate')
        job.args += ['-hklin', hklin, '-hklout', hklout, '-colin', colin]
        if colano:
            job.args += ['-colano', colano]
        job.parser = '_ctruncate_parser'
        return job

    def cad(self, data_in, hklout, keys):
        assert isinstance(data_in, list)
        hklin_args = []
        labin = []
        for n, (hklin, labels) in enumerate(data_in):
            labels = [a for a in labels if a not in ('H', 'K', 'L')]
            hklin_args += ['HKLIN%d' % (n+1), hklin]
            labin.append('labin file %d ' % (n+1) +
                         ' '.join('E%d=%s' % (k+1, label)
                                  for k, label in enumerate(labels)))
        job = ccp4_job(self, 'cad', logical={}, ki=(labin + keys),
                       parser='_cad_parser')
        job.args += hklin_args + ['HKLOUT', hklout]
        return job

    def pdbset(self, xyzin, xyzout, cell):
        return ccp4_job(self, 'pdbset', logical=locals(),
                        ki=['cell %g %g %g %g %g %g' % cell])

    def refmac5(self, hklin, xyzin, hklout, xyzout, labin, libin, keys):
        inp = ['labin %s' % labin] + keys.splitlines()
        #inp += ['free 6']  # for testing
        job = ccp4_job(self, 'refmac5', logical=locals(), ki=inp,
                       parser='_refmac_parser')
        words = keys.split()
        ref_type = '?'
        for n, w in enumerate(words[:-2]):
            if w == 'refinement' and words[n+1] == 'type':
                ref_type = words[n+2][:5]
            elif w == 'ridge':
                ref_type = 'jelly'
        job.name += ' ' + ref_type
        job.data['ncyc'] = -1
        for n, w in enumerate(words[:-1]):
            if w.startswith('ncyc'):
                job.data['ncyc'] = int(words[n+1])
        return job

    def get_final_refinement_job(self):
        for job in reversed(self.jobs):
            if job.name in ('refmac5 restr', 'refmac5 jelly'):
                return job

    def findwaters(self, pdbin, hklin, f, phi, pdbout, sigma=2.0):
        job = Job(self, 'findwaters')
        job.args += ['--pdbin', pdbin, '--hklin', hklin, '--f', f, '--phi', phi,
                     '--pdbout', pdbout, '--sigma', '%g' % sigma]
        return job

    def find_blobs(self, hklin, xyzin, sigma=1.0):
        # for now search in PATH (which normally includes CBIN)
        job = Job(self, utils.syspath('find-blobs'))
        job.args += ['-c', '-s%g' % sigma, hklin, xyzin]
        job.parser = '_find_blobs_parser'
        return job

    # alternative to find_blobs
    def gemmi_blobs(self, hklin, xyzin, sigma=1.0):
        job = Job(self, utils.cbin('gemmi'))
        job.name += ' ' + 'blobs'
        job.args += ['blobs', '--dimple', '--sigma=%g' % sigma, hklin, xyzin]
        job.parser = '_gemmi_blobs_parser'
        return job

    def mtz2sca(self, mtzin, scaout):
        job = Job(self, utils.syspath('mtz2sca'))
        job.args += [mtzin, scaout]
        return job

    def shelxc(self, scain, cell, symmetry):
        job = Job(self, utils.syspath('shelxc'))
        name = os.path.splitext(scain)[0]
        job.args += [name]

        job.std_input = '\n'.join([
            'SAD %s' % scain,
            'CELL %s %s %s %s %s %s' % cell,
            'SPAG %s' % symmetry,
        ])
        return job

    def anode(self, name):
        job = Job(self, utils.syspath('anode'))
        job.args += [name]
        job.parser = '_anode_parser'
        return job

    def shelx2map(self, hklin, mapout, symmetry):
        job = Job(self, utils.syspath('shelx2map'))
        job.args += [hklin, "-o", mapout, "-s", symmetry]
        return job

    def rwcontents(self, xyzin):
        return ccp4_job(self, 'rwcontents', logical=dict(xyzin=xyzin),
                        parser='_rwcontents_parser')

    def coot_py(self, script_text):
        job = Job(self, coots.find_path())
        job.args += ['--python', '--no-graphics', '--no-guano']
        script_text += '\ncoot_real_exit(0)'
        # On some Wincoot installations coot-real.exe is started from
        # runwincoot.bat directly, and on some as "start ... coot-real ...".
        # There is no way afaics to pipe stdin to coot-real.
        if os.name == 'nt':
            helper_path = self.path('r3d.py')
            with open(helper_path, 'w') as f:
                f.write(script_text)
            job.args.append(helper_path)
        else:
            job.std_input = script_text
        return job

    def render_r3d(self, name, img_format):
        assert img_format is not None
        job = Job(self, utils.syspath('render'))
        # render writes normal output to stderr (and nothing to stdout)
        job.out.file_extension = 'out'
        job.err.file_extension = 'log'
        job.args += ['-'+img_format, '%s.%s' % (name, img_format)]
        job.stdin_file = name + '.r3d'
        job.parser = ' %s.%s' % (name, img_format)
        return job

    def copy_uncompressed(self, src, dst):
        src_fullpath = self.path(src)
        dst_fullpath = self.path(dst)
        if src.endswith('.gz'):
            with gzip.open(src_fullpath, 'rb') as fsrc:
                content = fsrc.read()
            with open(dst_fullpath, 'wb') as fdst:
                fdst.write(content)
        else:
            try:
                shutil.copyfile(src_fullpath, dst_fullpath)
            except shutil.Error:  # == SameFileError in Python 3.4+
                pass

    def delete_files(self, filenames):
        for f in filenames:
            path = self.path(f)
            if os.path.exists(path):
                try:
                    os.remove(path)
                except OSError as e:
                    utils.put_error(e)


def open_pickled_workflow(file_or_dir):
    if os.path.isdir(file_or_dir):
        pkl = os.path.join(file_or_dir, PICKLE_FILENAME)
    else:
        pkl = file_or_dir
    if not os.path.exists(pkl):
        utils.put_error('workflow data file not found',
                        'No such file or directory: %s' % pkl)
        sys.exit(1)
    f = open(pkl, 'rb')
    try:
        return pickle.load(f)
    except pickle.UnpicklingError:
        utils.put_error('"Unpickling" failed',
                        'Maybe this is not a pickle file: %s' % pkl)
        sys.exit(1)

def _write_workflow_steps(wf, output):
    for n, job in enumerate(wf.jobs):
        output.write('\n%3d %-15s' % (n+1, job.name))
        if job.started:
            started_at = time.localtime(job.started)
            output.write(time.strftime(' %Y-%m-%d %H:%M', started_at))
            output.write(' %7.1fs' % job.total_time)
    output.write('\n')

def show_workflow_info(wf, mesg_dict):
    sys.stdout.write('%s\n' % wf)
    sys.stdout.write('Command:\n' + ' '.join(shlex.quote(a) for a in wf.argv))
    _write_workflow_steps(wf, sys.stdout)
    sys.stderr.write("""
To see details, specify step(s):
%(prog)s info %(output_dir)s STEPS

To re-run selected steps (for debugging):
%(prog)s repeat %(output_dir)s [STEPS]

where STEPS is one or more numbers or a range (examples: 1,2 4-6 8-)
""" % mesg_dict)

def show_job_info(job):
    sys.stdout.write('%s\n' % job)
    sys.stdout.write(job.args_as_str() + '\n')
    if job.total_time:
        sys.stdout.write('Total time: %.1fs\n' % job.total_time)
    if job.parser and job.parse():
        sys.stdout.write('Output summary: %s\n' % job.parse())
    if job.out.saved_to:
        sys.stdout.write('stdout: %s\n' % job.out.summary())
    if job.err.saved_to:
        sys.stdout.write('stderr: %s\n' % job.err.summary())


def parse_steps(args, wf):
    jobs = []
    for arg in args:
        try:
            for s in arg.split(','):
                if '-' in s:
                    a_, b_ = s.split('-')
                    a = (int(a_) if a_ != '' else 1)
                    b = (int(b_) if b_ != '' else len(wf.jobs))
                    if a == 0 or b == 0:
                        raise ValueError()
                    jobs += [wf.jobs[n-1] for n in range(a, b+1)]
                else:
                    jobs.append(wf.jobs[int(s)-1])
        except (ValueError, IndexError) as e:
            sys.stderr.write('Invalid step number(s): %s\n(%s)\n' % (arg, e))
            sys.exit(1)
    return jobs

def parse_workflow_commands():
    prog = __package__ or os.path.basename(sys.argv[0])
    args = sys.argv[1:]
    if not args or args[0] not in ('info', 'repeat'):
        return False
    if len(args) == 1:
        sys.stderr.write('Specify output_dir.\n')
        return True

    # it's handy to use "/my/path/05-cad.log" as "/my/path" "5"
    ext = os.path.splitext(args[1])[1]
    if os.path.isfile(args[1]) and ext in ('.log', '.err'):
        dirname, basename = os.path.split(args[1])
        args[1:2] = [dirname, basename.split('-')[0]]

    wf = open_pickled_workflow(args[1])
    steps = args[2:]
    if not steps:
        show_workflow_info(wf, dict(prog=prog, output_dir=args[1]))
        return True
    for job in parse_steps(steps, wf):
        if args[0] == 'info':
            show_job_info(job)
        elif args[0] == 'repeat':
            try:
                job.data = {}  # reset data from parsing
                job.run()
                utils.comment('\n')
            except JobError as e:
                utils.put_error(e.msg, comment=e.note)
                sys.exit(1)
    return True

commands_help = """\
All files are stored in the specified output directory.
For quick summary (after running the program): %(prog)s info OUTPUT_DIR
"""

if __name__ == '__main__':
    def test_parser(name, logfile):
        parser = globals()['_%s_parser' % name]
        job = Job(None, name)
        job.out.que = queue.Queue()
        with open(logfile) as f:
            for line in f:
                job.out.que.put(line)
        parser(job)
        for k in sorted(job.data.keys()):
            print('%s  %s' % (k, job.data[k]))

    assert len(sys.argv) == 3
    test_parser(sys.argv[1], logfile=sys.argv[2])
