import os
import sys
from subprocess import Popen, PIPE
import errno
import threading
import Queue
import time
import cPickle as pickle
#import pickle
import c4.mtz
import c4.pdb

def put(text):
    sys.stdout.write(text)

def put_green(text):
    if hasattr(sys.stdout, 'isatty') and sys.stdout.isatty():
        put("\033[92m%s\033[0m" % text)
    else:
        put(text)

def put_error(err, comment=""):
    if hasattr(sys.stderr, 'isatty') and sys.stderr.isatty():
        error_line = "Error: \033[91m%s\033[0m." % err # in bold red
    else:
        error_line = "Error: %s." % err
    sys.stderr.write(error_line + "\n" + comment)


class Job:
    def __init__(self, workflow, prog):
        self.name = os.path.basename(prog) or prog # only used to show info
        self.workflow = workflow
        self.args = [prog]
        self.std_input = ""
        # the rest is set after the job is run
        self.out = [] # either string or list of lines
        self.err = [] # either string or list of lines
        self.started = None # will be set to time.time() at start
        self.total_time = None # will be set when job ends
        self.parser = None
        self.data = {} # parsing helpers

    def __str__(self):
        desc = "Job %s" % self.name
        if self.started:
            desc += time.strftime(" %Y-%m-%d %H:%M",
                                            time.localtime(self.started))
        return desc

    def run(self):
        return self.workflow.run_job(job=self, show_progress=True)

    def parse(self):
        if self.parser:
            p = getattr(self, self.parser)
            return p()
        else:
            # generic non-parser shows output size if the job is finished
            if self.total_time is None:
                return ""
            def size(out):
                if not out:               return "-"
                elif type(out) is list:   return "%dL" % len(out)
                else:                     return "%.1fkB" % (len(out) / 1024.)
            ret = "stdout:%7s" % size(self.out)
            if self.err:
                ret += " stderr: %s" % size(self.err)
            return ret

    def _read_output(self):
        while True and self.data['out_q'] is not None:
            try:
                line = self.data['out_q'].get_nowait()
            except Queue.Empty:
                break
            self.out.append(line)
            yield line

    def _find_blobs_parser(self):
        if "blobs" not in self.data:
            sys.stdout.write("\n")
            self.data["blobs"] = []
        for line in self._read_output():
            sys.stdout.write(line)
            if line.startswith("#"):
                sp = line.split()
                score = float(sp[5])
                if score > 150:
                    xyz = tuple(float(x.strip(",()")) for x in sp[-3:])
                    self.data["blobs"].append(xyz)

    def _refmac_parser(self):
        t = self.data
        if "cycle" not in t:
            t["cycle"] = 0
            t["free_r"] = t["overall_r"] = 0.
        for line in self._read_output():
            if line.startswith("Free R factor"):
                t['free_r'] = float(line.split('=')[-1])
            elif line.startswith("Overall R factor"):
                t['overall_r'] = float(line.split('=')[-1])
            elif (line.startswith("     Rigid body cycle =") or
                  line.startswith("     CGMAT cycle number =")):
                t['cycle'] = int(line.split('=')[-1])
        return "cycle %(cycle)2d/%(ncyc)d   "\
                "R-free / R = %(free_r).4f / %(overall_r).4f" % t


def ccp4_job(workflow, prog, logical=None, input="", add_end=True):
    """Handle traditional convention for arguments of CCP4 programs.
    logical is dictionary with where keys are so-called logical names,
    input string or list of lines that are to be passed though stdin
    add_end adds "end" as the last line of stdin
    """
    assert os.environ.get("CCP4")
    full_path = os.path.join(os.environ["CCP4"], "bin", prog)
    job = Job(workflow, full_path)
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


_jobindex_fmt = "%3d "
_jobname_fmt = "%-15s"
_elapsed_fmt = "%5.1fs  "

def _print_elapsed(job, event):
    while not event.wait(0.5):
        p = job.parse()
        if p is not None:
            text = (_elapsed_fmt % (time.time() - job.started)) + p
            put(text)
            sys.stdout.flush()
            put("\b"*len(text))


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

def _write_output(output, filename):
    with open(filename, "w") as f:
        if type(output) is list:
            for line in output:
                f.write(line)
        else:
            f.write(output)


class Workflow:
    def __init__(self, output_dir):
        self.output_dir = os.path.abspath(output_dir)
        self.jobs = []
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)

    def __str__(self):
        return "Workflow with %d jobs @ %s" % (len(self.jobs), self.output_dir)

    def pickle_jobs(self, filename="workflow.pickle"):
        os.chdir(self.output_dir)
        with open(filename, "wb") as f:
            pickle.dump(self, f, -1)

    def run_job(self, job, show_progress):
        if not hasattr(sys.stdout, 'isatty') or not sys.stdout.isatty():
            show_progress = False
        os.chdir(self.output_dir)
        self.jobs.append(job)
        put(_jobindex_fmt % len(self.jobs))
        put_green(_jobname_fmt % job.name)
        sys.stdout.flush()
        job.started = time.time()
        #job.args[0] = "true" # for debugging
        try:
            p = Popen(job.args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        except OSError as e:
            if e.errno == errno.ENOENT:
                raise RuntimeError("Program not found: %s\n" % job.args[0])
            else:
                raise

        if show_progress:
            event = threading.Event()
            progress_thread = threading.Thread(target=_print_elapsed,
                                               args=(job, event))
            progress_thread.daemon = True
            progress_thread.start()

        if job.parser:
            out_t, job.data['out_q'] = _start_enqueue_thread(p.stdout)
            err_t, job.data['err_q'] = _start_enqueue_thread(p.stderr)
            try:
                p.stdin.write(job.std_input)
            except IOError as e:
                put("\nWarning: passing std input to %s failed.\n" % job.name)
                if e.errno not in (errno.EPIPE, e.errno != errno.EINVAL):
                    raise
            p.stdin.close()
            out_t.join()
            err_t.join()
            p.wait()
            # nothing is written to the queues at this point
            # parse what's left in the queues
            job.parse()
            # take care of what is left by the parser
            while not job.data['out_q'].empty():
                job.out.append(job.data['out_q'].get_nowait())
            while not job.data['err_q'].empty():
                job.err.append(job.data['err_q'].get_nowait())
            job.data['out_q'] = None
            job.data['err_q'] = None
        else:
            job.out, job.err = p.communicate(input=job.std_input)

        if show_progress:
            event.set()
            progress_thread.join()

        job.total_time = time.time() - job.started
        retcode = p.poll()
        put(_elapsed_fmt % job.total_time)
        put("%s\n" % (job.parse() or ""))
        self._write_logs(job)
        if retcode:
            raise RuntimeError("`%s' failed. All args:\n%s" % (
                        job.args[0], " ".join('"%s"' % a for a in job.args)))
        return job

    def _write_logs(self, job):
        log_basename = "%02d-%s" % (len(self.jobs), job.name.replace(" ","_"))
        if job.out:
            _write_output(job.out, "%s.log" % log_basename)
            job.out = "-> %s.log" % log_basename
        if job.err:
            _write_output(job.err, "%s.err" % log_basename)
            job.err = "-> %s.err" % log_basename

    def change_pdb_cell(self, xyzin, xyzout, cell):
        #for now using pdbset
        self.pdbset(xyzin=xyzin, xyzout=xyzout, cell=cell).run()

    def read_pdb_metadata(self, xyzin):
        return c4.pdb.read_metadata(xyzin)

    def read_mtz_metadata(self, hklin):
        return c4.mtz.read_metadata(hklin)

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
        job.data['ncyc'] = -1
        for n, w in enumerate(words[:-1]):
            if w.startswith("ncyc"):
                job.data['ncyc'] = int(words[n+1])
        job.parser = "_refmac_parser"
        return job

    def find_blobs(self, mtz, pdb, sigma=1.0):
        this_dir = os.path.abspath(os.path.dirname(__file__))
        job = Job(self, os.path.join(this_dir, "find-blobs"))
        job.args += ["-s%g" % sigma, mtz, pdb]
        job.parser = "_find_blobs_parser"
        return job


if __name__ == '__main__':
    usage = "Usage: python -m c4.workflow output_dir [N]\n"
    if len(sys.argv) < 2:
        sys.stderr.write(usage)
        sys.exit(0)
    if os.path.isdir(sys.argv[1]):
        pkl = os.path.join(sys.argv[1], "workflow.pickle")
    elif os.path.exists(sys.argv[1]):
        pkl = sys.argv[1]
    else:
        sys.stderr.write("No such file: %s\n" % sys.argv[1])
        sys.exit(1)
    with open(pkl) as f:
        wf = pickle.load(f)
    if len(sys.argv) == 2:
        sys.stdout.write("%s\n" % wf)
        for n, job in enumerate(wf.jobs):
            sys.stdout.write("%3d %s\n" % (n+1, job))
        sys.stderr.write("To see details, add job number(s).\n")
    else:
        for job_str in sys.argv[2:]:
            job_nr = int(job_str) - 1
            job = wf.jobs[job_nr]
            sys.stdout.write("%s\n" % job)
            sys.stdout.write(" ".join('"%s"' % a for a in job.args))
            if job.std_input:
                sys.stdout.write(" << EOF\n%s\nEOF" % job.std_input)
            sys.stdout.write("\nTotal time: %.1fs\n" % job.total_time)
            sys.stdout.write("Parser: %s\n" % job.parser)
            if job.out and type(job.out) is str and len(job.out) < 160:
                sys.stdout.write("stdout: %s\n" % job.out)
            if job.err and type(job.err) is str and len(job.err) < 160:
                sys.stdout.write("stderr: %s\n" % job.err)
