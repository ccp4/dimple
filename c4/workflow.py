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
import c4.utils
import c4.mtz
import c4.pdb
import c4.coot


_jobindex_fmt = "%3d "
_jobname_fmt = "%-15s"
_elapsed_fmt = "%5.1fs  "


# heh, python from MSYS2 is handy, but needs some monkey-patching
if sys.platform == 'msys':
    old_opjoin = os.path.join
    def new_opjoin(*args):
        for n in range(len(args)-1, -1, -1):
            if n == 0 or os.path.isabs(args[n]):
                return old_opjoin(*args[n:])
    os.path.join = new_opjoin


class JobError(Exception):
    def __init__(self, msg, note=None):
        self.msg = msg
        self.note = note


class Output:
    "storage for Job's stdout/stderr"
    def __init__(self, role):
        self.role = role  # "out" or "err"
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

    def save_output(self, output_dir, filename, remove_long_list=True):
        if self.lines:
            with open(os.path.join(output_dir, filename), "w") as f:
                for line in self.lines:
                    f.write(line)
            self.saved_to = filename
            c4.utils.log_value("std"+self.role, filename)
            if remove_long_list and len(self.lines) > 5:
                self.lines = []

    def summary(self):
        n = len(self.lines)
        if n < 3:
            return "".join(self.lines)
        elif self.saved_to:
            return "-> %s" % self.saved_to
        else:
            return "".join(self.lines[:3]) + ("%s more lines" % (n-3))


class Job:
    def __init__(self, workflow, prog):
        self.name = os.path.basename(prog) or prog  # only used to show info
        self.workflow = workflow
        self.args = [prog]
        self.std_input = ""
        self.stdin_file = None # if set, it overwrites std_input
        # the rest is set after the job is run
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
                ret = "[%d] %s" % (len(self.out.lines), line.rstrip())
                return ret[:50].ljust(50)

            ret = "stdout:%11s" % self.out.size_as_str()
            if self.err:
                ret += " stderr: %s" % self.err.size_as_str()
            if preview_mode:
                ret = ret.ljust(50)
            return ret

        elif self.parser[0] == ' ':
            return self.parser

        else:
            p = globals()[self.parser]
            return p(self)


def _format(fmt, arg):
    return (fmt % arg) if arg else ""

# parsers for various programs
def _find_blobs_parser(job):
    if "blobs" not in job.data:
        #sys.stdout.write("\n")
        job.data["blobs"] = []
        job.data["scores"] = []
    for line in job.out.read_line():
        #sys.stdout.write(line)
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
    scores = job.data["scores"]
    if scores:
        return "Blob scores: " + " ".join("%.0f" % sc for sc in scores)
    else:
        return ""

def _refmac_parser(job):
    if "cycle" not in job.data:
        job.data["cycle"] = 0
        job.data["free_r"] = job.data["overall_r"] = 0.
        job.data["ini_free_r"] = job.data["ini_overall_r"] = 0
        job.data["selected_lines"] = []
        # iter_*_r values were added in dimple 1.5
        job.data["iter_free_r"] = []
        job.data["iter_overall_r"] = []
    selected = job.data["selected_lines"]
    for line in job.out.read_line():
        if line.startswith("Free R factor"):
            job.data['free_r'] = float(line.split('=')[-1])
            job.data['iter_free_r'].append(job.data['free_r'])
            if not job.data['ini_free_r']:
                job.data['ini_free_r'] = job.data['free_r']
        elif line.startswith("Overall R factor"):
            job.data['overall_r'] = float(line.split('=')[-1])
            job.data['iter_overall_r'].append(job.data['overall_r'])
            if not job.data['ini_overall_r']:
                job.data['ini_overall_r'] = job.data['overall_r']
        elif (line.startswith("     Rigid body cycle =") or
              line.startswith("     CGMAT cycle number =")):
            job.data['cycle'] = int(line.split('=')[-1])
        elif line.startswith(" $TEXT:Result: $$ Final results $$") or (
                selected and not selected[-1].startswith(" $$")):
            selected.append(line)
    return "%2d/%d   R/Rfree  %.4f/%.4f  ->  %.4f/%.4f" % (
            job.data["cycle"], job.ncyc,
            job.data["ini_overall_r"], job.data["ini_free_r"],
            job.data["overall_r"], job.data["free_r"])

# example:
#     Alternative reindexing        Lklhd      CC     R(E^2)    Number Cell_deviation
#    [-k,-l,h+k+l]           0.079    0.029    0.506     61253      2.99
_POINTLESS_ALTREINDEX_MATCH = re.compile(r"^\s+\[[hkl+-, ]+\][ \t0-9+-.eE]+$")

def _pointless_parser(job):
    # the line with 'reflections copied' is the last we read
    # to avoid reading 'Alternative reindexing' duplicated in the summary
    if "refl_out" not in job.data:
        for line in job.out.read_line():
            if line.startswith("Maximum resolution used:"):
                job.data["resol"] = float(line.split(":")[1])
            elif line.startswith("Number of reflections:"):
                job.data["refl_in"] = int(line.split(":")[1])
            elif _POINTLESS_ALTREINDEX_MATCH.match(line):
                reindex = job.data.setdefault('alt_reindex', [])
                s = line.split()
                reindex.append({'op': s[0], 'cc': s[2], 'cell_deviat': s[-1]})
            elif line.startswith("   Cell:") and 'output_cell' not in job.data:
                s = line.split()[1:]
                job.data["output_cell"] = tuple(float(i) for i in s)
            elif "reflections copied to output file" in line:
                job.data["refl_out"] = int(line.split()[0])
                break
    return "resol. %4s A   #refl: %5s -> %s" % (
            _format("%.2f", job.data.get("resol")),
            job.data.get("refl_in", ""),
            job.data.get("refl_out", ""))


def _truncate_parser(job):
    for line in job.out.read_line():
        if line.startswith(' Least squares straight line gives:'):
            b_str = line.partition(':')[-1].strip(' \tB=').split()[0]
            job.data["B-factor"] = float(b_str)
    return " B=%4s" % _format("%.1f", job.data.get("B-factor"))

def _ctruncate_parser(job):
    d = job.data
    for line in job.out.read_line():
        if line.startswith("$GRAPHS: Wilson plot - estimated B factor ="):
            d["B-factor"] = float(line.partition('=')[-1].split()[0])
        elif line.startswith("Eigenvalue ratios:"):
            d["eigval-ratios"] = tuple(float(i) for i in line.split()[-3:])
        elif line.startswith("L statistic ="):
            d["L-test"] = float(line.partition('=')[-1].split()[0])
    return " B=%4s   aniso %14s   L-test:%s" % (
            _format("%.1f", d.get("B-factor")),
            _format("%.2f:%.2f:%.2f", d.get("eigval-ratios")),
            _format("%.2f", d.get("L-test")))


def ccp4_job(workflow, prog, logical=None, input="", parser=None, add_end=True):
    """Handle traditional convention for arguments of CCP4 programs.
    logical is dictionary with where keys are so-called logical names,
    input string or list of lines that are to be passed though stdin
    add_end adds "end" as the last line of stdin
    """
    job = Job(workflow, c4.utils.cbin(prog))
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


def _print_elapsed(job, event):
    while not event.wait(0.5):
        p = job.parse()
        if p is not None:
            text = (_elapsed_fmt % (time.time() - job.started)) + p
            c4.utils.put(text)
            sys.stdout.flush()
            c4.utils.put("\b"*len(text))
            c4.utils.reset_color()


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
            c4.utils.put("\nWarning: passing input to %s failed.\n" % job.name)
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
        self.from_job = from_job # skip jobs before from_job (for testing)
        if from_job >= 1:
            try:
                self.repl_jobs = self.unpickle_jobs().jobs
            except:
                self.repl_jobs = None
        self.dry_run = False
        self.argv = sys.argv
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)
        # this can seriously affect Refmac compiled with GFortran
        bad_var = os.getenv('GFORTRAN_UNBUFFERED_ALL')
        if bad_var and bad_var[0] not in ('0', 'n', 'N'):
            c4.utils.put_error(
                    '$GFORTRAN_UNBUFFERED_ALL may terribly slow down Refmac',
                    comment='It is unset internally in dimple.')
            del os.environ['GFORTRAN_UNBUFFERED_ALL']
        # avoid html-like crap in the output of CCP4 program
        os.environ['CCP_SUPPRESS_HTML'] = '1'


    def __str__(self):
        return "Workflow with %d jobs @ %s" % (len(self.jobs), self.output_dir)

    def path(self, rel_path):
        return os.path.join(self.output_dir, rel_path)

    def pickle_jobs(self, filename="workflow.pickle"):
        with open(self.path(filename), "wb") as f:
            pickle.dump(self, f, -1)

    def unpickle_jobs(self, filename="workflow.pickle"):
        with open(self.path(filename), "rb") as f:
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
            raise JobError("\nKeyboardInterrupt while running %s" % job.name,
                           note=job.args_as_str())
        finally:
            job.total_time = time.time() - job.started
        return process.poll()

    def run_job(self, job, show_progress):
        if not hasattr(sys.stdout, 'isatty') or not sys.stdout.isatty():
            show_progress = False
        self.jobs.append(job)
        job_num = len(self.jobs)
        c4.utils.put(_jobindex_fmt % job_num)
        c4.utils.put_green(_jobname_fmt % job.name)
        sys.stdout.flush()
        c4.utils.log_section(job.name)

        job_idx = len(self.jobs) - 1
        if job_idx < self.from_job - 1: # from_job is 1-based
            if self.repl_jobs and len(self.repl_jobs) > job_idx:
                old_job = self.repl_jobs[job_idx]
                if old_job.name == job.name:
                    job = old_job
                    c4.utils.put("unpickled\n")
                    c4.utils.log_value("not_run", "unpickled")
                    self.jobs[-1] = job
                else:
                    c4.utils.put("skipped (mismatch)\n")
                    c4.utils.log_value("not_run", "unpickled/mismatch")
            else:
                c4.utils.put("skipped\n")
                c4.utils.log_value("not_run", "skipped")
            return job

        job.started = time.time()
        c4.utils.log_time("start_time", job.started)
        if job.stdin_file:
            c4.utils.log_value("stdin", job.stdin_file)
        elif job.std_input:
            c4.utils.log_value("input", job.std_input)
        c4.utils.log_value("prog", job.args[0])
        c4.utils.log_value("args",
                           " ".join(pipes.quote(a) for a in job.args[1:]))
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
            progress_thread = threading.Thread(target=_print_elapsed,
                                               args=(job, event))
            progress_thread.daemon = True
            progress_thread.start()

        try:
            if job.parser is not None or show_progress:
                _run_and_parse(process, job)
            else:
                _just_run(process, job)
        except KeyboardInterrupt:
            raise JobError("\nKeyboardInterrupt while running %s" % job.name,
                           note=job.args_as_str())
        finally:
            if show_progress:
                event.set()
                progress_thread.join()
            end_time = time.time()
            job.total_time = end_time - job.started
            c4.utils.log_time("end_time", end_time)
            retcode = process.poll()
            c4.utils.put(_elapsed_fmt % job.total_time)
            parse_output = job.parse()
            c4.utils.put("%s\n" % (parse_output or ""))
            if parse_output:
                c4.utils.log_value("info", parse_output)
            self._write_logs(job)
            for k, v in job.data.iteritems():
                if k == "selected_lines":
                    v = "\n" + "".join(v) # selected_lines have newlines
                c4.utils.log_value(k, v)
        if retcode:
            c4.utils.log_value("exit_status", retcode)
            all_args = job.args_as_str()
            notes = [all_args, ""]
            if job.out.saved_to:
                notes += ["stdout -> %s/%s" % (self.output_dir,
                                               job.out.saved_to)]
            if job.err:
                notes += ["stderr:", job.err.summary()]
            raise JobError("%s failed (exit status %d)" % (job.name, retcode),
                           note="\n".join(notes))
        return job

    def _write_logs(self, job):
        log_basename = "%02d-%s" % (len(self.jobs), job.name.replace(" ","_"))
        job.out.save_output(self.output_dir, "%s.log" % log_basename)
        job.err.save_output(self.output_dir, "%s.err" % log_basename)

    def remove_hetatm(self, xyzin, xyzout):
        with open(xyzout, "wb") as out:
            return c4.pdb.remove_hetatm(xyzin, out)

    def read_pdb_metadata(self, xyzin):
        return c4.pdb.read_metadata(self.path(xyzin))

    def read_mtz_metadata(self, hklin):
        return c4.mtz.read_metadata(self.path(hklin))

    def get_protein_mw(self, xyzin):
        return c4.pdb.get_protein_mw(self.path(xyzin))

    def molrep(self, f, m, keys=""):
        job = Job(self, c4.utils.cbin("molrep"))
        job.args += ["-f", f, "-m", m]
        if keys:
            job.args.append("-i")
            job.std_input = keys.strip() + "\nend"
        return job

    def phaser(self, hklin, labin, mode, script, root):
        lines = ["MODE %s" % mode,
                 "HKLIN %s" % hklin,
                 "LABIN %s" % labin
                ] + script.splitlines() + [
                 "ROOT %s" % root ]
        job = ccp4_job(self, "phaser", input=lines)
        return job

    def pointless(self, hklin, xyzin, hklref=None, hklout=None, keys=""):
        return ccp4_job(self, "pointless", logical=locals(), input=keys,
                        parser="_pointless_parser")

    def unique(self, hklout, cell, symmetry, resolution,
               labout="F=F_UNIQUE SIGF=SIGF_UNIQUE"):
        return ccp4_job(self, "unique", logical=locals(),
                        input=["cell %g %g %g %g %g %g" % tuple(cell),
                               "symmetry '%s'" % symmetry,
                               "resolution %.3f" % resolution,
                               "labout %s" % labout])

    def freerflag(self, hklin, hklout, keys=""):
        return ccp4_job(self, "freerflag", logical=locals(), input=keys)

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
        assert type(hklin) is list
        job = ccp4_job(self, "cad", logical={}, input=keys)
        # is hklinX only for cad?
        for n, name in enumerate(hklin):
            job.args += ["HKLIN%d" % (n+1), name]
        job.args += ["HKLOUT", hklout]
        return job

    def pdbset(self, xyzin, xyzout, cell):
        return ccp4_job(self, "pdbset", logical=locals(),
                        input=["cell %g %g %g %g %g %g" % cell])

    def refmac5(self, hklin, xyzin, hklout, xyzout, labin, labout, libin, keys):
        job = ccp4_job(self, "refmac5", logical=locals(),
                       input=(["labin %s" % labin, "labout %s" % labout] +
                              keys.splitlines()),
                       parser="_refmac_parser")
        words = keys.split()
        for n, w in enumerate(words[:-2]):
            if w == "refinement" and words[n+1] == "type":
                job.name += " " + words[n+2][:5]
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

    def find_blobs(self, mtz, pdb, sigma=1.0):
        # for now search in PATH (which normally includes CBIN)
        job = Job(self, c4.utils.syspath("find-blobs"))
        job.args += ["-c", "-s%g" % sigma, mtz, pdb]
        job.parser = "_find_blobs_parser"
        return job

    def coot_py(self, script_text):
        job = Job(self, c4.coot.find_path())
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
        job = Job(self, c4.utils.syspath("render"))
        job.args += ["-"+img_format, "%s.%s" % (name, img_format)]
        job.stdin_file = name+".r3d"
        job.parser = " -> %s.%s" % (name, img_format)
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
                    c4.utils.put_error(e)


def open_pickled_workflow(file_or_dir):
    if os.path.isdir(file_or_dir):
        pkl = os.path.join(file_or_dir, "workflow.pickle")
    else:
        pkl = file_or_dir
    if not os.path.exists(pkl):
        c4.utils.put_error("workflow data file not found",
                           "No such file or directory: %s" % pkl)
        sys.exit(1)
    f = open(pkl, "rb")
    return pickle.load(f)


def show_workflow_info(wf, mesg_dict):
    sys.stdout.write("%s\n" % wf)
    sys.stdout.write("Command: " + " ".join(pipes.quote(a) for a in wf.argv))
    for n, job in enumerate(wf.jobs):
        sys.stdout.write("\n%3d %-15s" % (n+1, job.name))
        if job.started:
            started_at = time.localtime(job.started)
            sys.stdout.write(time.strftime(" %Y-%m-%d %H:%M", started_at))
            sys.stdout.write(" %7.1fs" % job.total_time)
    sys.stdout.write("\n")
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
        except (ValueError, IndexError):
            raise ValueError("Invalid step number(s): %s" % arg)
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
        if len(args) <= 2:
            wf = open_pickled_workflow(args[1])
            sys.stderr.write("Specify steps. For complete re-run:\n%s\n"
                             % " ".join(pipes.quote(a) for a in wf.argv))
        else:
            wf = open_pickled_workflow(args[1])
            for job in parse_steps(args[2:], wf):
                try:
                    job.data = {}  # reset data from parsing
                    job.run()
                except JobError as e:
                    c4.utils.put_error(e.msg, comment=e.note)
                    sys.exit(1)
        return True

commands_help = """\
All files are stored in the specified output directory.
For quick summary (after running the program): %(prog)s info OUTPUT_DIR
"""

