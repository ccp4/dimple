import os
import sys
from subprocess import Popen, PIPE
import threading
import Queue
import time

from . import utils


class Job:
    def __init__(self, workflow, prog):
        self.name = prog # name is only used to show info
        self.workflow = workflow
        self.args = [prog]
        self.std_input = ""
        # the rest is set after the job is run
        self.out = [] # either string or list of lines
        self.err = [] # either string or list of lines
        self.started = None # will be set to time.time() at start
        self.elapsed = -1.0 # will be set when job ends
        self.parser = None
        self.tmp = {} # parsing helpers
        self.err_queue = None # parsing helpers

    def run(self):
        self.workflow.run_job(job=self, show_progress=True)
    def run_and_parse(self):
        self.workflow.run_job(job=self, show_progress=True)


class Ccp4Job(Job):
    def __init__(self, workflow, prog, logical=None, input="", add_end=True):
        Job.__init__(self, workflow, prog)
        if logical:
            self.add_logical(logical)
        if input or add_end:
            self.add_input(input, add_end)

    def add_logical(self, loc):
        "command line arguments that follow CCP4 logical-names convention"
        for a in ["hklin", "hklout", "hklref", "xyzin", "xyzout"]:
            if loc.get(a):
                self.args.extend([a.upper(), loc[a]])

    def add_input(self, lines, add_end):
        "passing `keyworded' arguments through stdin is also CCP4 convention"
        if isinstance(lines, basestring):
            lines = lines.splitlines()
        lines = [a.strip() for a in lines if a and not a.isspace()]
        if add_end and not (lines and lines[-1].lower() == "end"):
            lines.append("end")
        if self.std_input:
            self.std_input += "\n"
        self.std_input += "\n".join(lines)

    def _refmac_parser(self):
        t = self.tmp
        if "cycle" not in t:
            t["cycle"] = 0
            t["free_r"] = t["overall_r"] = 0.
        while True:
            try:
                line = t['out_q'].get_nowait()
                self.out.append(line)
                if line.startswith("Free R factor"):
                    t['free_r'] = float(line.split('=')[-1])
                elif line.startswith("Overall R factor"):
                    t['overall_r'] = float(line.split('=')[-1])
                elif (line.startswith("     Rigid body cycle =") or
                      line.startswith("     CGMAT cycle number =")):
                    t['cycle'] = int(line.split('=')[-1])
            except Queue.Empty:
                break
        return "#%(cycle)-2d  R-free / R = %(free_r).4f / %(overall_r).4f" % t


def out_size(out):
    if not out:
        return "-"
    elif type(out) == list:
        return "%dL" % len(out)
    else:
        return "%.1fkB" % (len(out) / 1024.)

def log(text):
    sys.stdout.write(text)

def _print_elapsed(job, event):
    while not event.wait(0.5):
        text = "elapsed:%5.1fs" % (time.time() - job.started)
        if job.parser:
            text += "  " + job.parser()
        log(text)
        sys.stdout.flush()
        log("\b"*len(text))


def start_enqueue_thread(file_obj):
    def enqueue_lines(f, q):
        for line in iter(f.readline, b''):
            q.put(line)
        f.close()
    que = Queue.Queue()
    thr = threading.Thread(target=enqueue_lines, args=(file_obj, que))
    thr.daemon = True
    thr.start()
    return thr, que

class Workflow:
    def __init__(self, output_dir):
        self.output_dir = os.path.abspath(output_dir)
        self.job_log = []
        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)

    def run_job(self, job, show_progress, dry_run=False):
        os.chdir(self.output_dir)
        self.job_log.append(job)
        log("%3d  %-30s" % (len(self.job_log), utils.green(job.name)))
        sys.stdout.flush()
        job.started = time.time()
        if dry_run:
            log("--- DRY RUN\n%s\n%s\n---\n" % (job.args, job.std_input))
        else:
            p = Popen(job.args, stdin=PIPE, stdout=PIPE, stderr=PIPE)
            if show_progress:
                event = threading.Event()
                progress_thread = threading.Thread(target=_print_elapsed,
                                                   args=(job, event))
                progress_thread.daemon = True
                progress_thread.start()
                if job.parser:
                    out_t, job.tmp['out_q'] = start_enqueue_thread(p.stdout)
                    err_t, job.tmp['err_q'] = start_enqueue_thread(p.stderr)
                    p.stdin.write(job.std_input)
                    p.stdin.close()
                    out_t.join()
                    err_t.join()
                    p.wait()
                    # nothing is written to the queues at this point
                    while not job.tmp['out_q'].empty():
                        job.out.append(job.tmp['out_q'].get_nowait())
                    while not job.tmp['err_q'].empty():
                        job.err.append(job.tmp['err_q'].get_nowait())
                else:
                    job.out, job.err = p.communicate(input=job.std_input)
                event.set()
                progress_thread.join()
            else:
                job.out, job.err = p.communicate(input=job.std_input)
            job.elapsed = time.time() - job.started
            retcode = p.poll()
            if retcode:
                raise RuntimeError("Program failed: %s" % job.args)
        log("elapsed:%5.1fs  " % job.elapsed)
        if not job.parser:
            log("stdout:%7s stderr: %s" % (out_size(job.out),
                                           out_size(job.err)))
        log("\n")

    def change_pdb_cell(self, xyzin, xyzout, cell):
        #for now using pdbset
        self.pdbset(xyzin=xyzin, xyzout=xyzout, cell=cell).run()

    def pointless(self, hklin, xyzin, hklref=None, hklout=None, keys=""):
        return Ccp4Job(self, "pointless", logical=locals(), input=keys)

    def mtzdump(self, hklin, keys=""):
        return Ccp4Job(self, "mtzdump", logical=locals())

    def unique(self, hklout, cell, symmetry, resolution,
               labout="F=F_UNIQUE SIGF=SIGF_UNIQUE"):
        return Ccp4Job(self, "unique", logical=locals(),
                       input=["cell %g %g %g %g %g %g" % cell,
                              "symmetry '%s'" % symmetry,
                              "resolution %.3f" % resolution,
                              "labout %s" % labout])

    def freerflag(self, hklin, hklout):
        return Ccp4Job(self, "freerflag", logical=locals())

    def reindex(self, hklin, hklout, symmetry):
        return Ccp4Job(self, "reindex", logical=locals(),
                       input=["symmetry '%s'" % symmetry,
                              "reindex h,k,l"])

    def truncate(self, hklin, hklout, labin, labout):
        return Ccp4Job(self, "truncate", logical=locals(),
                       input=["labin %s" % labin, "labout %s" % labout])

    def cad(self, hklin, hklout, keys):
        assert type(hklin) is list
        job = Ccp4Job(self, "cad", logical={}, input=keys)
        # is hklinX only for cad?
        for n, name in enumerate(hklin):
            job.args += ["HKLIN%d" % (n+1), name]
        job.args += ["HKLOUT", hklout]
        return job

    def pdbset(self, xyzin, xyzout, cell):
        return Ccp4Job(self, "pdbset", logical=locals(),
                       input=["cell %g %g %g %g %g %g" % cell])

    def refmac5(self, hklin, xyzin, hklout, xyzout, labin, labout, keys):
        job = Ccp4Job(self, "refmac5", logical=locals(),
                      input=(["labin %s" % labin, "labout %s" % labout] +
                              keys.splitlines()))
        words = keys.split()
        for n, w in enumerate(words[:-1]):
            if w == "refinement" and words[n+1] == "type":
                job.name += " " + words[n+2][:5]
            elif w.startswith("ncyc"):
                job.name += "*" + words[n+1]
        job.parser = job._refmac_parser
        return job

