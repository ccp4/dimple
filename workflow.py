import os
import sys
from subprocess import Popen, PIPE
import gzip
import errno
import pipes
import re
import threading
import Queue
import time
import cPickle as pickle
import shutil
from dimple import utils
from dimple import mtz
from dimple import pdb
from dimple import coots


_jobindex_fmt = "%3d "
_jobname_fmt = "%-15s"
_elapsed_fmt = "%5.1fs  "

PICKLE_FILENAME = "workflow.pickle"


# heh, python from MSYS2 is handy, but needs some monkey-patching
if sys.platform == 'msys':
    old_opjoin = os.path.join
    def new_opjoin(*args):
        for n in range(len(args)-1, -1, -1):
            if n == 0 or os.path.isabs(args[n]):
                return old_opjoin(*args[n:])
    os.path.join = new_opjoin


class JobError(Exception):
    def __init__(self, msg, note=None): # pylint: disable=super-init-not-called
        self.msg = msg
        self.note = note


class Output:
    "storage for Job's stdout/stderr"
    def __init__(self, role):
        self.role = role  # "out" or "err"
        self.file_extension = {"out":"log", "err":"err"}[role]
        self.lines = []
        self.saved_to = None
        self.que = None

    def __nonzero__(self):
        return bool(self.lines or self.saved_to)

    def size_as_str(self):
        if self.lines:
            return "%d lines" % len(self.lines)
        else:
            return "-      "

    def read_line(self):
        while self.que is not None:
            try:
                line = self.que.get_nowait()
            except Queue.Empty:
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
            with open(os.path.join(output_dir, filename), "w") as f:
                for line in self.lines:
                    f.write(line)
            self.saved_to = filename
            utils.log_value("std"+self.role, filename)
            if remove_long_list and len(self.lines) > 5:
                self.lines = []

    def summary(self):
        n = len(self.lines)
        if 0 < n <= 3:
            return "".join(self.lines)
        elif self.saved_to:
            return "-> %s" % self.saved_to
        elif n == 0:
            return ""
        else: # n > 3
            return "".join(self.lines[:3]) + ("%s more lines" % (n-3))


class Job:
    def __init__(self, workflow, prog):
        self.name = os.path.basename(prog) or prog  # only used to show info
        self.workflow = workflow
        self.args = [prog]
        self.std_input = ""
        self.stdin_file = None # if set, it overwrites std_input
        # the rest is set after the job is run
        self.exit_status = None
        self.out = Output("out")
        self.err = Output("err")
        self.started = None  # will be set to time.time() at start
        self.total_time = None  # will be set when job ends
        # possible values: None (default),
        #                  'preview' (stdout preview),
        #                  string starting with space (' ') that is just shown,
        #                  or name of global function that parses output
        self.parser = None
        # job-specific data from output parsing
        self.data = {}

    def __repr__(self):
        if self.started:
            t = time.strftime(" %Y-%m-%d %H:%M", time.localtime(self.started))
        else:
            t = ""
        return "<Job %s%s>" % (self.name, t)

    def args_as_str(self):
        s = " ".join(pipes.quote(a) for a in self.args)
        if self.stdin_file:
            s += " < " + self.stdin_file
        elif self.std_input:
            s += " << EOF\n%s\nEOF" % self.std_input
        return s

    def run(self):
        return self.workflow.run_job(job=self, show_progress=True)

    def parse(self):
        preview_mode = (self.parser == "preview")
        if self.parser is None or preview_mode:
            # generic non-parser
            line = ""
            for line in self.out.read_line():
                pass

            if preview_mode and line:
                # we only remove \e[1m (bold), it can be generalized if needed
                trimmed = line.strip().replace('\033[1m', '')
                return "[%d] %-44.44s" % (len(self.out.lines), trimmed)

            ret = "stdout:%11s" % self.out.size_as_str()
            if self.err:
                ret += " stderr: %s" % self.err.size_as_str()
            if preview_mode:
                ret = ret.ljust(50)
            return ret

        elif self.parser == '' or self.parser[0] == ' ':
            return self.parser

        else:
            p = globals()[self.parser]
            return p(self)


def _format(fmt, arg):
    return (fmt % arg) if arg else ""

# parsers for various programs
def _find_blobs_parser(job):
    if "blobs" not in job.data:
        job.data["blobs"] = []
        job.data["scores"] = []
    for line in job.out.read_line():
        if line.startswith("#"):
            sp = line.split(None, 6)
            score = float(sp[5])
            if True: #score > 150: XXX: better scoring may be needed
                x, y, z = sp[-1].strip("() \t\r\n").split(",")
                job.data["blobs"].append((float(x), float(y), float(z)))
                job.data["scores"].append(score)
        elif line.startswith("Protein mass center:"):
            sp = line.split("(")[1].rstrip("\r\n )").split(",")
            ctr = tuple(float(x) for x in sp)
            job.data["center"] = ctr
        elif line.startswith("Density std.dev"):
            job.data["density_info"] = line.strip()
    scores = job.data["scores"]
    if scores:
        return "Blob scores: " + " ".join("%.0f" % sc for sc in scores)
    elif "density_info" in job.data:
        return 'searching with' + job.data["density_info"].split(',')[1]
    else:
        return ""

def _rwcontents_parser(job):
    d = job.data
    for line in job.out.read_line():
        if line.startswith(' Cell volume:'):
            d["volume"] = float(line.split(':')[-1])
        elif line.startswith(' Molecular Weight of protein:'):
            d["weight"] = float(line.split(':')[-1])
        if line.startswith(' The Matthews Coefficient is :'):
            Vm = float(line.split(':')[-1])
            d["Vm"] = Vm
            # 1.23 is used in phaser/src/Composition.cc
            d["solvent_percent"] = (1 - 1.23/Vm) * 100
    if 'volume' in d and 'weight' in d and 'Vm' in d and 'num_mol' not in d:
        d['num_mol'] = int(round(d['volume'] / d['weight'] / d['Vm']))
    return u"%d x %.0fkDa in %.fnm3  Vm=%.2f (%.0f%% of solvent)" % (
            d.get('num_mol', 0),
            d.get('weight', 0) / 1000, # Da -> kDa
            d.get('volume', 0) / 1000, # A^2 -> nm^3
            d.get('Vm', 0),
            d.get('solvent_percent', 0))

def _cad_parser(job):
    for line in job.out.read_line():
        if line.startswith(' Final Total of Unique records to HKLOUT ='):
            job.data["refl_out"] = int(line.split('=')[1])
    return "#refl -> %s" % job.data.get("refl_out", "")

def _refmac_parser(job):
    if "cycle" not in job.data:
        # ini_free_r, free_r and iter_free_r are set optionally
        job.data["cycle"] = 0
        job.data["selected_lines"] = []
        # iter_*_r values were added in dimple 1.5
        job.data["iter_overall_r"] = []
    selected = job.data["selected_lines"]
    for line in job.out.read_line():
        if line.startswith("Free R factor"):
            job.data['free_r'] = float(line.split('=')[-1])
            if 'ini_free_r' not in job.data:
                job.data['ini_free_r'] = job.data['free_r']
                job.data['iter_free_r'] = []
            job.data['iter_free_r'].append(job.data['free_r'])
        elif line.startswith("Overall R factor"):
            job.data['overall_r'] = float(line.split('=')[-1])
            if 'ini_overall_r' not in job.data:
                job.data['ini_overall_r'] = job.data['overall_r']
            job.data['iter_overall_r'].append(job.data['overall_r'])
        elif (line.startswith("     Rigid body cycle =") or
              line.startswith("     CGMAT cycle number =")):
            job.data['cycle'] = int(line.split('=')[-1])
        elif line.startswith(" $TEXT:Result: $$ Final results $$") or (
                selected and not selected[-1].startswith(" $$")):
            selected.append(line)
    cycle_str = "%2d/%d" % (job.data["cycle"], job.ncyc)
    if 'ini_overall_r' in job.data:
        if 'ini_free_r' in job.data:
            return "%s   R/Rfree  %.4f/%.4f  ->  %.4f/%.4f" % (
                    cycle_str,
                    job.data["ini_overall_r"], job.data["ini_free_r"],
                    job.data["overall_r"], job.data["free_r"])
        return "%s   R  %.4f  ->  %.4f" % (
                cycle_str, job.data["ini_overall_r"], job.data["overall_r"])
    return cycle_str

# example:
#     Alternative reindexing        Lklhd      CC     R(E^2)    Number Cell_deviation
#    [-k,-l,h+k+l]           0.079    0.029    0.506     61253      2.99
_POINTLESS_ALTREINDEX_MATCH = re.compile(r"^\s+\[[hkl+, -]+\][ \t0-9+.eE-]+$")

def _float_or_nan(s):
    try:
        x = float(s)
    except ValueError:
        x = float('nan')
    return x

def _pointless_parser(job):
    # the line with 'reflections copied' is the last we read
    # to avoid reading 'Alternative reindexing' duplicated in the summary
    if "refl_out" not in job.data:
        for line in job.out.read_line():
            if line.startswith("Maximum resolution used:"):
                job.data["resol"] = float(line.split(":")[1])
            elif line.startswith("Number of reflections:"):
                job.data["refl_ref"] = int(line.split(":")[1])
            elif _POINTLESS_ALTREINDEX_MATCH.match(line):
                s = line.split()
                job.data.setdefault('alt_reindex', []).append(
                        {'op': s[0],
                         'cc': _float_or_nan(s[2]),
                         'cell_deviat': _float_or_nan(s[-1])})
            elif line.startswith("   Cell:") and 'output_cell' not in job.data:
                s = line.split()[1:]
                job.data["output_cell"] = tuple(float(i) for i in s)
            elif "reflections copied to output file" in line:
                job.data["refl_out"] = int(line.split()[0])
                break
    resol_txt = _format("%.2f", job.data.get("resol"))
    refl_out = job.data.get("refl_out", "")
    #refl_ref = job.data.get("refl_ref", "")
    return "resol. %4s A   #refl: %5s" % (resol_txt, refl_out)


def _phaser_parser(job):
    d = job.data
    for line in job.out.read_line():
        if line.startswith('*** Phaser Module:'):
            d['status'] = '[%s]' % line[19:70].strip().lower()
        if line.startswith('   SOLU SET ') and 'LLG=' in line:
            d['status'] = line[12:].strip()
        if line.startswith('   Sorry - No solution'):
            d['status'] = line.strip()
        if line.startswith('   SOLU SPAC '):
            d['SG'] = line[13:].strip()
    return "%-48s" % d.get('status', '')

def _truncate_parser(job):
    for line in job.out.read_line():
        if line.startswith(' Least squares straight line gives:'):
            b_str = line.partition(':')[-1].strip(' \tB=').split()[0]
            job.data["B-factor"] = float(b_str)
    return "B=%4s" % _format("%.1f", job.data.get("B-factor"))

def _ctruncate_parser(job):
    d = job.data
    for line in job.out.read_line():
        if line.startswith("$GRAPHS: Wilson plot - estimated B factor ="):
            d["B-factor"] = float(line.partition('=')[-1].split()[0])
        elif line.startswith("Eigenvalue ratios:"):
            d["eigval-ratios"] = tuple(float(i) for i in line.split()[-3:])
        elif line.startswith("L statistic ="):
            d["L-test"] = float(line.partition('=')[-1].split()[0])
    return "B=%4s   aniso %14s   L-test:%s" % (
            _format("%.1f", d.get("B-factor")),
            _format("%.2f:%.2f:%.2f", d.get("eigval-ratios")),
            _format("%.2f", d.get("L-test")))


def ccp4_job(workflow, prog, logical=None,
             input="",  # pylint: disable=redefined-builtin
             parser=None, add_end=True):
    """Handle traditional convention for arguments of CCP4 programs.
    logical is dictionary with where keys are so-called logical names,
    input string or list of lines that are to be passed though stdin
    add_end adds "end" as the last line of stdin
    """
    job = Job(workflow, utils.cbin(prog))
    if logical:
        for a in ["hklin", "hklout", "hklref", "xyzin", "xyzout", "libin"]:
            if logical.get(a):
                job.args += [a.upper(), logical[a]]
    lines = (input.splitlines() if isinstance(input, basestring) else input)
    stripped = [a.strip() for a in lines if a and not a.isspace()]
    if add_end and not (stripped and stripped[-1].lower() == "end"):
        stripped.append("end")
    job.std_input = "\n".join(stripped)
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
    que = Queue.Queue()
    thr = threading.Thread(target=enqueue_lines, args=(file_obj, que))
    thr.daemon = True
    thr.start()
    return thr, que

def _get_input_as_string(job):
    if job.stdin_file:
        path = job.workflow.path(job.stdin_file)
        try:
            return open(path, "rb").read()
        except IOError:
            raise JobError("cannot read input from: %s" % job.stdin_file)
    else:
        return job.std_input

def _just_run(process, job):
    job_input = _get_input_as_string(job)
    out, err = process.communicate(input=job_input)
    job.out.lines = out.splitlines(True)
    job.err.lines = err.splitlines(True)

def _run_and_parse(process, job):
    try:
        # job.*.que can be used by parsers (via Output.read_line() or directly)
        out_t, job.out.que = _start_enqueue_thread(process.stdout)
        err_t, job.err.que = _start_enqueue_thread(process.stderr)
        try:
            job_input = _get_input_as_string(job)
            process.stdin.write(job_input)
        except IOError as e:
            utils.put("\nWarning: passing input to %s failed.\n" % job.name)
            if e.errno not in (errno.EPIPE, e.errno != errno.EINVAL):
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
            utils.put_error(
                    '$GFORTRAN_UNBUFFERED_ALL may terribly slow down Refmac',
                    comment='It is unset internally in dimple.')
            del os.environ['GFORTRAN_UNBUFFERED_ALL']
        # avoid html-like crap in the output of CCP4 program
        os.environ['CCP_SUPPRESS_HTML'] = '1'


    def __str__(self):
        return "Workflow with %d jobs @ %s" % (len(self.jobs), self.output_dir)

    def path(self, rel_path):
        return os.path.join(self.output_dir, rel_path)

    def dump_pickle(self):
        with open(self.path(PICKLE_FILENAME), "wb") as f:
            pickle.dump(self, f, -1)

    def load_pickle(self):
        with open(self.path(PICKLE_FILENAME), "rb") as f:
            return pickle.load(f)

    def silently_run_job(self, job):
        job.started = time.time()
        try:
            process = Popen(job.args, stdin=PIPE, stdout=PIPE, stderr=PIPE,
                            cwd=self.output_dir)
        except OSError as e:
            if e.errno == errno.ENOENT:
                raise JobError("Program not found: %s\n" % job.args[0])
            else:
                raise
        try:
            _just_run(process, job)
        except KeyboardInterrupt:
            raise JobError("KeyboardInterrupt while running %s" % job.name,
                           note=job.args_as_str())
        finally:
            job.total_time = time.time() - job.started
        return process.poll()

    def run_job(self, job, show_progress, new_line=True):
        if not hasattr(sys.stdout, 'isatty') or not sys.stdout.isatty():
            show_progress = False
        self.jobs.append(job)
        job_num = len(self.jobs)
        if new_line:
            utils.put("\n" + _jobindex_fmt % job_num)
            utils.put_green(_jobname_fmt % job.name)
        else:
            utils.put(" / %d" % job_num)
        sys.stdout.flush()
        utils.log_section(job.name)

        job_idx = len(self.jobs) - 1
        if job_idx < self.from_job - 1: # from_job is 1-based
            # unpickle or skip
            if self.repl_jobs and len(self.repl_jobs) > job_idx:
                old_job = self.repl_jobs[job_idx]
                if old_job.name == job.name:
                    job = old_job
                    utils.put("unpickled")
                    utils.log_value("not_run", "unpickled")
                    self.jobs[-1] = job
                else:
                    utils.put("skipped (mismatch)")
                    utils.log_value("not_run", "unpickled/mismatch")
            else:
                utils.put("skipped")
                utils.log_value("not_run", "skipped")
            return job

        job.started = time.time()
        utils.log_time("start_time", job.started)
        if job.stdin_file:
            utils.log_value("stdin", job.stdin_file)
        elif job.std_input:
            utils.log_value("input", job.std_input)
        utils.log_value("prog", job.args[0])
        utils.log_value("args", " ".join(pipes.quote(a) for a in job.args[1:]))
        #job.args[0] = "true"  # for debugging
        try:
            process = Popen(job.args, stdin=PIPE, stdout=PIPE, stderr=PIPE,
                            cwd=self.output_dir)
        except OSError as e:
            if e.errno == errno.ENOENT:
                raise JobError("Program not found: %s\n" % job.args[0])
            else:
                raise

        if self.dry_run:
            return job

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
                _just_run(process, job)
        except KeyboardInterrupt:
            raise JobError("KeyboardInterrupt while running %s" % job.name,
                           note=job.args_as_str())
        finally:
            if show_progress:
                event.set()
                progress_thread.join()
            end_time = time.time()
            job.total_time = end_time - job.started
            utils.log_time("end_time", end_time)
            job.exit_status = process.poll()
            if new_line:
                utils.put(_elapsed_fmt % job.total_time)
            parse_output = job.parse()
            utils.put("%s" % (parse_output or ""))
            if parse_output:
                utils.log_value("info", parse_output)
            self._write_logs(job)
            for k, v in job.data.iteritems():
                if k == "selected_lines":
                    v = "\n" + "".join(v) # selected_lines have newlines
                utils.log_value(k, v)
        if job.exit_status:
            utils.log_value("exit_status", job.exit_status)
            all_args = job.args_as_str()
            notes = [all_args, ""]
            if job.out.saved_to:
                notes += ["stdout -> %s/%s" % (self.output_dir,
                                               job.out.saved_to)]
            if job.err:
                notes += ["stderr:", job.err.summary()]
            raise JobError("%s failed (exit status %d)" % (job.name,
                                                           job.exit_status),
                           note="\n".join(notes))
        return job

    def _write_logs(self, job):
        log_basename = "%02d-%s" % (len(self.jobs), job.name.replace(" ","_"))
        for output in (job.out, job.err):
            output.save_output(self.output_dir, log_basename)

    def remove_hetatm(self, xyzin, xyzout, remove_all):
        with open(self.path(xyzout), "wb") as out:
            return pdb.remove_hetatm(self.path(xyzin), out, remove_all)

    def read_pdb_metadata(self, xyzin):
        if xyzin not in self.file_info:
            self.file_info[xyzin] = pdb.read_metadata(self.path(xyzin))
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

    def molrep(self, f, m, keys=""):
        job = Job(self, utils.cbin("molrep"))
        job.args += ["-f", f, "-m", m]
        if keys:
            job.args.append("-i")
            job.std_input = keys.strip() + "\nend"
        return job

    def phaser_auto(self, hklin, labin, model, root, sg_alt="NONE",
                    solvent_percent=None):
        lines = ['MODE MR_AUTO',
                 'SEARCH METHOD FAST',
                 'SEARCH DEEP OFF',
                 'HKLIN "%s"' % hklin,
                 'LABIN %s' % labin,
                 'SGALTERNATIVE SELECT %s' % sg_alt]
        if solvent_percent:
            lines += ['COMPOSITION BY SOLVENT',
                      'COMPOSITION PERCENTAGE %f' % solvent_percent]
        lines += ('ENSEMBLE p PDBFILE %(pdb)s IDENTITY %(identity)g\n'
                  'SEARCH ENSEMBLE p NUM %(num)d' % model).splitlines()
        # For tNCS we go with what phaser does by default -- tNCS of order 2
        # are handled automatically. While we could specify tNCS for
        # pseudo-tripling/quadrupling of the cell (TNCS NMOL 3) I don't know
        # if it'd do more good or bad.
        lines += ["ROOT %s" % root]
        job = ccp4_job(self, "phaser", input=lines, parser="_phaser_parser")
        return job

    # functions below use logical=locals()
    # pylint: disable=unused-argument

    def pointless(self, hklin, xyzin, hklref=None, hklout=None, keys=""):
        return ccp4_job(self, "pointless", logical=locals(), input=keys,
                        parser="_pointless_parser")

    def unique(self, hklout, cell, symmetry, resolution,
               labout="F=F_UNIQUE SIGF=SIGF_UNIQUE"):
        return ccp4_job(self, "unique", logical=locals(),
                        input=["cell %g %g %g %g %g %g" % tuple(cell),
                               "symmetry '%s'" % symmetry,
                               "resolution %.3f" % resolution,
                               "labout %s" % labout],
                        parser="")

    def freerflag(self, hklin, hklout, keys=""):
        return ccp4_job(self, "freerflag", logical=locals(), input=keys,
                        parser="")

    #def reindex(self, hklin, hklout, symmetry):
    #    return ccp4_job(self, "reindex", logical=locals(),
    #                    input=["symmetry '%s'" % symmetry,
    #                           "reindex h,k,l"])

    def truncate(self, hklin, hklout, labin, labout):
        return ccp4_job(self, "truncate", logical=locals(),
                        input=["labin %s" % labin, "labout %s" % labout,
                               "NOHARVEST"],
                        parser="_truncate_parser")

    def ctruncate(self, hklin, hklout, colin):
        job = Job(self, "ctruncate")
        job.args += ["-hklin", hklin, "-hklout", hklout, "-colin", colin]
        job.parser = "_ctruncate_parser"
        return job

    def cad(self, hklin, hklout, keys):
        assert isinstance(hklin, list)
        job = ccp4_job(self, "cad", logical={}, input=keys,
                       parser="_cad_parser")
        # is hklinX only for cad?
        for n, name in enumerate(hklin):
            job.args += ["HKLIN%d" % (n+1), name]
        job.args += ["HKLOUT", hklout]
        return job

    def pdbset(self, xyzin, xyzout, cell):
        return ccp4_job(self, "pdbset", logical=locals(),
                        input=["cell %g %g %g %g %g %g" % cell])

    def refmac5(self, hklin, xyzin, hklout, xyzout, labin, labout, libin, keys):
        inp = ["labin %s" % labin, "labout %s" % labout] + keys.splitlines()
        #inp += ['free 6']  # for testing
        job = ccp4_job(self, "refmac5", logical=locals(), input=inp,
                       parser="_refmac_parser")
        words = keys.split()
        ref_type = "?"
        for n, w in enumerate(words[:-2]):
            if w == "refinement" and words[n+1] == "type":
                ref_type = words[n+2][:5]
            elif w == "ridge":
                ref_type = "jelly"
        job.name += " " + ref_type
        job.ncyc = -1
        for n, w in enumerate(words[:-1]):
            if w.startswith("ncyc"):
                job.ncyc = int(words[n+1])
        return job

    def findwaters(self, pdbin, hklin, f, phi, pdbout, sigma=2.0):
        job = Job(self, "findwaters")
        job.args += ["--pdbin", pdbin, "--hklin", hklin, "--f", f, "--phi", phi,
                     "--pdbout", pdbout, "--sigma", "%g" % sigma]
        return job

    def find_blobs(self, hklin, xyzin, sigma=1.0):
        # for now search in PATH (which normally includes CBIN)
        job = Job(self, utils.syspath("find-blobs"))
        job.args += ["-c", "-s%g" % sigma, hklin, xyzin]
        job.parser = "_find_blobs_parser"
        return job

    def rwcontents(self, xyzin):
        return ccp4_job(self, "rwcontents", logical=dict(xyzin=xyzin),
                        parser="_rwcontents_parser")

    def coot_py(self, script_text):
        job = Job(self, coots.find_path())
        job.args += ["--python", "--no-graphics", "--no-guano"]
        script_text += "\ncoot_real_exit(0)"
        # On some Wincoot installations coot-real.exe is started from
        # runwincoot.bat directly, and on some as "start ... coot-real ...".
        # There is no way afaics to pipe stdin to coot-real.
        if os.name == 'nt':
            helper_path = self.path("r3d.py")
            with open(helper_path, "w") as f:
                f.write(script_text)
            job.args.append(helper_path)
        else:
            job.std_input = script_text
        job.parser = "preview"
        return job

    def render_r3d(self, name, img_format="png"):
        job = Job(self, utils.syspath("render"))
        # render writes normal output to stderr (and nothing to stdout)
        job.out.file_extension = "out"
        job.err.file_extension = "log"
        job.args += ["-"+img_format, "%s.%s" % (name, img_format)]
        job.stdin_file = name+".r3d"
        job.parser = " %s.%s" % (name, img_format)
        return job

    def copy_uncompressed(self, src, dst):
        src_fullpath = self.path(src)
        dst_fullpath = self.path(dst)
        if src.endswith(".gz"):
            with gzip.open(src_fullpath, 'rb') as fsrc:
                content = fsrc.read()
            with open(dst_fullpath, 'wb') as fdst:
                fdst.write(content)
        else:
            shutil.copy2(src_fullpath, dst_fullpath)

    def delete_files(self, filenames):
        for f in filenames:
            path = self.path(f)
            if os.path.exists(path):
                try:
                    os.remove(path)
                except OSError, e:
                    utils.put_error(e)


def open_pickled_workflow(file_or_dir):
    if os.path.isdir(file_or_dir):
        pkl = os.path.join(file_or_dir, PICKLE_FILENAME)
    else:
        pkl = file_or_dir
    if not os.path.exists(pkl):
        utils.put_error("workflow data file not found",
                        "No such file or directory: %s" % pkl)
        sys.exit(1)
    f = open(pkl, "rb")
    return pickle.load(f)

def _write_workflow_steps(wf, output):
    for n, job in enumerate(wf.jobs):
        output.write("\n%3d %-15s" % (n+1, job.name))
        if job.started:
            started_at = time.localtime(job.started)
            output.write(time.strftime(" %Y-%m-%d %H:%M", started_at))
            output.write(" %7.1fs" % job.total_time)
    output.write("\n")

def show_workflow_info(wf, mesg_dict):
    sys.stdout.write("%s\n" % wf)
    sys.stdout.write("Command: " + " ".join(pipes.quote(a) for a in wf.argv))
    _write_workflow_steps(wf, sys.stdout)
    sys.stderr.write("""
To see details, specify step(s):
%(prog)s info %(output_dir)s STEPS

To re-run selected steps (for debugging):
%(prog)s repeat %(output_dir)s [STEPS]

where STEPS is one or more numbers or a range (examples: 1,2 4-6 8-)
""" % mesg_dict)

def show_job_info(job):
    sys.stdout.write("%s\n" % job)
    sys.stdout.write(job.args_as_str() + "\n")
    if job.total_time:
        sys.stdout.write("Total time: %.1fs\n" % job.total_time)
    if job.parser and job.parse():
        sys.stdout.write("Output summary: %s\n" % job.parse())
    if job.out.saved_to:
        sys.stdout.write("stdout: %s\n" % job.out.summary())
    if job.err.saved_to:
        sys.stdout.write("stderr: %s\n" % job.err.summary())


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
            sys.stderr.write("Invalid step number(s): %s\n(%s)\n" % (arg, e))
            sys.exit(1)
    return jobs

def parse_workflow_commands():
    prog = os.path.basename(sys.argv[0])
    args = sys.argv[1:]
    if not args:
        return False
    if args[0] == 'info':
        if len(args) == 1:
            sys.stderr.write("Specify output_dir.\n")
            return True
        wf = open_pickled_workflow(args[1])
        if len(args) == 2:
            show_workflow_info(wf, dict(prog=prog, output_dir=args[1]))
        else:
            for job in parse_steps(args[2:], wf):
                show_job_info(job)
        return True

    if args[0] == 'repeat':
        if len(args) == 1:
            sys.stderr.write("Specify output_dir.\n")
        elif len(args) <= 2:
            wf = open_pickled_workflow(args[1])
            sys.stderr.write("Specify steps from the list "
                             "(you can use ranges, e.g: 1,2 4-6 8-):")
            _write_workflow_steps(wf, sys.stderr)
            sys.stderr.write("For complete re-run:\n%s\n"
                             % " ".join(pipes.quote(a) for a in wf.argv))
        else:
            wf = open_pickled_workflow(args[1])
            for job in parse_steps(args[2:], wf):
                try:
                    job.data = {}  # reset data from parsing
                    job.run()
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
        job.out.que = Queue.Queue()
        with open(logfile) as f:
            for line in f:
                job.out.que.put(line)
        parser(job)
        for k in sorted(job.data.keys()):
            print k, job.data[k]

    assert len(sys.argv) == 3
    test_parser(sys.argv[1], logfile=sys.argv[2])
