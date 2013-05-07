import os
import sys
from subprocess import Popen, PIPE
import errno
import pipes
import threading
import Queue
import time
import cPickle as pickle
import c4.utils
import c4.mtz
import c4.pdb
import c4.coot


_jobindex_fmt = "%3d "
_jobname_fmt = "%-15s"
_elapsed_fmt = "%5.1fs  "


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

    def clear(self):
        self.lines = None

    def append(self, line):
        self.lines.append(line)

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
        while not self.que.empty():
            self.lines.append(self.que.get_nowait())
        self.que = None

    def save_output(self, output_dir, filename, remove_long_list=True):
        if self.lines:
            with open(os.path.join(output_dir, filename), "w") as f:
                for line in self.lines:
                    f.write(line)
            self.saved_to = filename
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
        # the rest is set after the job is run
        self.out = Output("out")
        self.err = Output("err")
        self.started = None  # will be set to time.time() at start
        self.total_time = None  # will be set when job ends
        # output parsing helpers
        self.parser = None
        # job-specific data from output parsing
        self.data = {}

    def __str__(self):
        desc = "Job %s" % self.name
        if self.started:
            desc += time.strftime(" %Y-%m-%d %H:%M",
                                  time.localtime(self.started))
        return desc

    def args_as_str(self):
        s = " ".join(pipes.quote(a) for a in self.args)
        if self.std_input:
            s += " << EOF\n%s\nEOF" % self.std_input
        return s

    def run(self):
        return self.workflow.run_job(job=self, show_progress=True)

    def parse(self):
        if self.parser:
            p = globals()[self.parser]
            return p(self)
        else:
            # generic non-parser
            for line in self.out.read_line():
                pass
            ret = "stdout:%11s" % self.out.size_as_str()
            if self.err:
                ret += " stderr: %s" % self.err.size_as_str()
            return ret


# parsers for various programs
def _find_blobs_parser(job):
    if "blobs" not in job.data:
        #sys.stdout.write("\n")
        job.data["blobs"] = []
        job.data["scores"] = []
    for line in job.out.read_line():
        #sys.stdout.write(line)
        if line.startswith("#"):
            sp = line.split()
            score = float(sp[5])
            if True: #score > 150: XXX: better scoring may be needed
                xyz = tuple(float(x.strip(",()")) for x in sp[-3:])
                job.data["blobs"].append(xyz)
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
        job.data["summary"] = []
    summary = job.data["summary"]
    for line in job.out.read_line():
        if line.startswith("Free R factor"):
            job.data['free_r'] = float(line.split('=')[-1])
        elif line.startswith("Overall R factor"):
            job.data['overall_r'] = float(line.split('=')[-1])
        elif (line.startswith("     Rigid body cycle =") or
              line.startswith("     CGMAT cycle number =")):
            job.data['cycle'] = int(line.split('=')[-1])
        elif line.startswith(" $TEXT:Result: $$ Final results $$") or (
                summary and not summary[-1].startswith(" $$")):
            summary.append(line)
    return "cycle %2d/%d   R-free / R = %.4f / %.4f" % (
        job.data["cycle"], job.ncyc, job.data["free_r"], job.data["overall_r"])


def ccp4_job(workflow, prog, logical=None, input="", add_end=True):
    """Handle traditional convention for arguments of CCP4 programs.
    logical is dictionary with where keys are so-called logical names,
    input string or list of lines that are to be passed though stdin
    add_end adds "end" as the last line of stdin
    """
    job = Job(workflow, c4.utils.full_path_of(prog))
    if logical:
        for a in ["hklin", "hklout", "hklref", "xyzin", "xyzout"]:
            if logical.get(a):
                job.args.extend([a.upper(), logical[a]])
    lines = (input.splitlines() if isinstance(input, basestring) else input)
    stripped = [a.strip() for a in lines if a and not a.isspace()]
    if add_end and not (stripped and stripped[-1].lower() == "end"):
        stripped.append("end")
    if job.std_input:
        job.std_input += "\n"
    job.std_input += "\n".join(stripped)
    return job


def _print_elapsed(job, event):
    while not event.wait(0.5):
        p = job.parse()
        if p is not None:
            text = (_elapsed_fmt % (time.time() - job.started)) + p
            c4.utils.put(text)
            sys.stdout.flush()
            c4.utils.put("\b"*len(text))


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

def _just_run(process, job):
    out, err = process.communicate(input=job.std_input)
    job.out.lines = out.splitlines(True)
    job.err.lines = err.splitlines(True)

def _run_and_parse(process, job):
    # job.*.que can be used by parsers (via Output.read_line() or directly)
    out_t, job.out.que = _start_enqueue_thread(process.stdout)
    err_t, job.err.que = _start_enqueue_thread(process.stderr)
    try:
        process.stdin.write(job.std_input)
    except IOError as e:
        c4.utils.put("\nWarning: passing std input to %s failed.\n" % job.name)
        if e.errno not in (errno.EPIPE, e.errno != errno.EINVAL):
            raise
    process.stdin.close()
    out_t.join()
    err_t.join()
    process.wait()
    # nothing is written to the queues at this point
    # parse what's left in the queues
    job.parse()
    # take care of what is left by the parser
    job.out.finish_que()
    job.err.finish_que()


class Workflow:
    def __init__(self, output_dir):
        self.output_dir = os.path.abspath(output_dir)
        self.jobs = []
        self.from_job = 0 # skip jobs before N (useful for testing)
        self.dry_run = False
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)

    def __str__(self):
        return "Workflow with %d jobs @ %s" % (len(self.jobs), self.output_dir)

    def pickle_jobs(self, filename="workflow.pickle"):
        with open(os.path.join(self.output_dir, filename), "wb") as f:
            pickle.dump(self, f, -1)

    def run_job(self, job, show_progress):
        if not hasattr(sys.stdout, 'isatty') or not sys.stdout.isatty():
            show_progress = False
        self.jobs.append(job)
        c4.utils.put(_jobindex_fmt % len(self.jobs))
        c4.utils.put_green(_jobname_fmt % job.name)
        sys.stdout.flush()

        if self.from_job >= len(self.jobs) + 1:
            c4.utils.put("skipped\n")
            return job

        job.started = time.time()
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
            self._write_logs(job)
            # Queues could cause PicklingError, empty and delete them
            job.out.finish_que()
            job.err.finish_que()
            raise JobError("\nKeyboardInterrupt while running %s" % job.name,
                           note=job.args_as_str())

        if show_progress:
            event.set()
            progress_thread.join()

        job.total_time = time.time() - job.started
        retcode = process.poll()
        c4.utils.put(_elapsed_fmt % job.total_time)
        c4.utils.put("%s\n" % (job.parse() or ""))
        self._write_logs(job)
        if retcode:
            all_args = " ".join(pipes.quote(a) for a in job.args)
            notes = []
            if job.out.saved_to:
                notes = ["stdout -> %s/%s" % (self.output_dir,
                                              job.out.saved_to)]
            if job.err:
                notes += ["stderr:", job.err.summary()]
            raise JobError("Non-zero return value from:\n%s" % all_args,
                           note="\n".join(notes))
        return job

    def _write_logs(self, job):
        log_basename = "%02d-%s" % (len(self.jobs), job.name.replace(" ","_"))
        job.out.save_output(self.output_dir, "%s.log" % log_basename)
        job.err.save_output(self.output_dir, "%s.err" % log_basename)

    def change_pdb_cell(self, xyzin, xyzout, cell):
        #for now using pdbset
        self.pdbset(xyzin=xyzin, xyzout=xyzout, cell=cell).run()

    def remove_hetatm(self, xyzin, xyzout):
        with open(xyzout, "wb") as out:
            return c4.pdb.remove_hetatm(xyzin, out)

    def read_pdb_metadata(self, xyzin):
        return c4.pdb.read_metadata(os.path.join(self.output_dir, xyzin))

    def read_mtz_metadata(self, hklin):
        return c4.mtz.read_metadata(os.path.join(self.output_dir, hklin))

    def molrep(self, f, m):
        job = Job(self, "molrep")
        job.args.extend(["-f", f, "-m", m])
        return job

    def pointless(self, hklin, xyzin, hklref=None, hklout=None, keys=""):
        return ccp4_job(self, "pointless", logical=locals(), input=keys)

    def mtzdump(self, hklin, keys=""):
        return ccp4_job(self, "mtzdump", logical=locals())

    def unique(self, hklout, cell, symmetry, resolution,
               labout="F=F_UNIQUE SIGF=SIGF_UNIQUE"):
        return ccp4_job(self, "unique", logical=locals(),
                        input=["cell %g %g %g %g %g %g" % cell,
                               "symmetry '%s'" % symmetry,
                               "resolution %.3f" % resolution,
                               "labout %s" % labout])

    def freerflag(self, hklin, hklout):
        return ccp4_job(self, "freerflag", logical=locals())

    def reindex(self, hklin, hklout, symmetry):
        return ccp4_job(self, "reindex", logical=locals(),
                        input=["symmetry '%s'" % symmetry,
                               "reindex h,k,l"])

    def truncate(self, hklin, hklout, labin, labout):
        return ccp4_job(self, "truncate", logical=locals(),
                        input=["labin %s" % labin, "labout %s" % labout])

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

    def refmac5(self, hklin, xyzin, hklout, xyzout, labin, labout, keys):
        job = ccp4_job(self, "refmac5", logical=locals(),
                       input=(["labin %s" % labin, "labout %s" % labout] +
                              keys.splitlines()))
        words = keys.split()
        for n, w in enumerate(words[:-2]):
            if w == "refinement" and words[n+1] == "type":
                job.name += " " + words[n+2][:5]
        job.ncyc = -1
        for n, w in enumerate(words[:-1]):
            if w.startswith("ncyc"):
                job.ncyc = int(words[n+1])
        job.parser = "_refmac_parser"
        return job

    def findwaters(self, pdbin, hklin, f, phi, pdbout, sigma=2.0):
        job = Job(self, "findwaters")
        job.args += ["--pdbin", pdbin, "--hklin", hklin, "--f", f, "--phi", phi,
                     "--pdbout", pdbout, "--sigma", "%g" % sigma]
        return job

    def find_blobs(self, mtz, pdb, sigma=1.0):
        job = Job(self, c4.utils.full_path_of("find-blobs"))
        job.args += ["-c", "-s%g" % sigma, mtz, pdb]
        job.parser = "_find_blobs_parser"
        return job


    def write_coot_script(self, name, pdb=None, mtz=None,
                          center=None, toward=None):
        assert center is not None
        path = os.path.join(self.output_dir, name)
        with open(path, "w") as f:
            f.write(c4.coot.basic_script(pdb=pdb, mtz=mtz,
                                         center=center, toward=toward))


    def make_img(self, blobname, pdb=None, mtz=None, center=None, toward=None,
                 format="png"):
        names = c4.coot.generate_r3d(pdb=pdb, mtz=mtz, center=center,
                                     blobname=blobname, cwd=self.output_dir,
                                     toward=toward)
        render_path = c4.utils.find_in_path("render")
        if not render_path:
            c4.utils.put_error("No Raster3d, no pictures")
            return
        for basename in names:
            print "rendering %s/%s.%s" % (self.output_dir, basename, format)
            r3d = open(os.path.join(self.output_dir, basename+".r3d")).read()
            render = Popen([render_path, "-%s" % format,
                                         "%s.%s" % (basename, format)],
                           stdin=PIPE, stdout=PIPE, stderr=PIPE,
                           cwd=self.output_dir)
            render.communicate(input=r3d)


def open_pickled_workflow(file_or_dir):
    if os.path.isdir(file_or_dir):
        pkl = os.path.join(file_or_dir, "workflow.pickle")
    else:
        pkl = file_or_dir
    if not os.path.exists(pkl):
        c4.utils.put_error("Workflow data file not found",
                           "No such file: %s" % pkl)
        sys.exit(1)
    f = open(pkl)
    return pickle.load(f)

def show_info(wf, job_numbers):
    if not job_numbers:
        sys.stdout.write("%s\n" % wf)
        for n, job in enumerate(wf.jobs):
            sys.stdout.write("%3d %s\n" % (n+1, job))
        sys.stderr.write("To see details, add job number(s).\n")
    for job_nr in job_numbers:
        show_job_info(wf.jobs[job_nr])

def show_job_info(job):
    sys.stdout.write("%s\n" % job)
    sys.stdout.write(job.args_as_str())
    sys.stdout.write("\nTotal time: %.1fs\n" % job.total_time)
    if job.parser and job.parse():
        sys.stdout.write("Output summary: %s\n" % job.parse())
    if job.out.saved_to:
        sys.stdout.write("stdout: %s\n" % job.out.summary())
    if job.err.saved_to:
        sys.stdout.write("stderr: %s\n" % job.err.summary())


def repeat_jobs(wf, job_numbers):
    if not job_numbers:
        sys.stderr.write("Which job(s) to repeat?\n")
    for job_nr in job_numbers:
        job = wf.jobs[job_nr]
        job.data = {}  # reset data from parsing
        job.run()


def parse_workflow_commands():
    args = sys.argv[1:]
    if len(args) >= 2 and args[0] == "info":
        wf = open_pickled_workflow(args[1])
        job_numbers = [int(job_str)-1 for job_str in args[2:]]
        show_info(wf, job_numbers)
        return True

    if len(args) >= 2 and args[0] == "repeat":
        wf = open_pickled_workflow(args[1])
        job_numbers = [int(job_str)-1 for job_str in args[2:]]
        try:
            repeat_jobs(wf, job_numbers)
        except JobError as e:
            c4.utils.put_error(e.msg, comment=e.note)
            sys.exit(1)
        return True

