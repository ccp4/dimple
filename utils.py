import os
import sys
import platform
import subprocess
import time

_logfile = None
_logfile_sections = None
_screen_logfile = None

# To make testing easier, programs found in the directory of these scripts
# take precedence over programs in standard locations.
# So it's possible to just copy here, say, refmac5 binary and test it.
_dimple_dir = os.path.abspath(os.path.dirname(__file__))

# start log in ini-like format,
# readable for humans and easily parsed in Python:
#   import ConfigParser
#   log = ConfigParser.RawConfigParser()
#   log.read('dimple.log')
# and then, for example:
#   ini_free_r = log.getfloat('refmac5 restr', 'ini_free_r')
#   free_r = log.getfloat('refmac5 restr', 'free_r')
# the first section is [workflow], next ones correspond to jobs:
#   start_time = log.get(log.sections()[1], 'start_time') # first job start
#   end_time = log.get(log.sections()[-1], 'end_time') # last job end
def start_log(filename, output_dir):
    global _logfile  # pylint: disable=global-statement
    global _logfile_sections  # pylint: disable=global-statement
    _logfile = open(filename, "w")
    _logfile.write("# workflow log (compatible with Python ConfigParser)\n")
    _logfile_sections = set()
    log_section("workflow")
    log_value("host", platform.node())
    log_value("platform", platform.platform())
    log_value("cwd", os.getcwd())
    log_value("prog", sys.argv[0])
    if len(sys.argv) > 1:
        _logfile.write("args:\n")
        for arg in sys.argv[1:]:
            _logfile.write(" %s\n" % arg)
    log_value("output_dir", output_dir)
    log_value("CCP4", os.getenv("CCP4", ""))
    log_value("CCP4_SCR", os.getenv("CCP4_SCR", ""))
    _logfile.flush()

def _log_comment(text):
    global _logfile  # pylint: disable=global-variable-not-assigned
    if _logfile:
        _logfile.write("# ")
        _logfile.write(text.rstrip("\n").replace("\n", "\n# "))
        _logfile.write("\n")

def log_section(name):
    global _logfile  # pylint: disable=global-variable-not-assigned
    global _logfile_sections  # pylint: disable=global-variable-not-assigned
    if _logfile:
        if name in _logfile_sections:
            counter = 2
            while ('%s %d' % (name, counter)) in _logfile_sections:
                counter += 1
            name = '%s %d' % (name, counter)
        _logfile_sections.add(name)
        _logfile.write("\n[%s]\n" % name)
        _logfile.flush()

def log_value(key, value):
    global _logfile  # pylint: disable=global-variable-not-assigned
    if _logfile:
        # TODO: store data structures with json.dumps() for easier parsing
        value = str(value).rstrip().replace("\n", "\n ")
        _logfile.write("%s: %s\n" % (key, value))

def log_time(key, timestamp):
    log_value(key, time.strftime("%Y-%m-%d %H:%M:%S",
                                 time.localtime(timestamp)))

def start_log_screen(filename):
    global _screen_logfile  # pylint: disable=global-statement
    _screen_logfile = open(filename, 'a')
    _screen_logfile.write('\n')

def _log_screen(text):
    global _screen_logfile  # pylint: disable=global-variable-not-assigned
    if _screen_logfile:
        _screen_logfile.write(text)
        _screen_logfile.flush()

def put(text, ansi_code=None):
    _log_screen(text)
    if (ansi_code is not None and os.name != 'nt' and
            hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()):
        sys.stdout.write("\033[%dm%s\033[0m" % (ansi_code, text))
    else:
        sys.stdout.write(text)

def put_temporarily(text):
    sys.stdout.write(text)
    sys.stdout.flush()
    sys.stdout.write("\b"*len(text))

def put_green(text):
    put(text, ansi_code=92)

def comment(text):
    put(text)
    sys.stdout.flush()
    _log_comment(text)

def reset_color():
    if hasattr(sys.stdout, 'isatty') and sys.stdout.isatty():
        if os.name != 'nt':
            sys.stdout.write("\033[0m")

def put_error(err, comment=None):  # pylint: disable=redefined-outer-name
    _log_comment("Error: %s." % err)
    _log_screen("\nError: %s.\n" % err)
    if hasattr(sys.stderr, 'isatty') and sys.stderr.isatty():
        if os.name != 'nt':
            err = "\033[91m%s\033[0m" % err  # in bold red
    sys.stdout.flush()
    sys.stderr.write("\nError: %s.\n" % err)
    if comment is not None:
        _log_comment(comment)
        _log_screen(comment + "\n")
        sys.stderr.write(comment + "\n")


def check_prog(dirname, prog):
    """If prog(.exe) is in dirname return the path (without extension)."""
    exe = '.exe' if os.name == 'nt' else ''
    path = os.path.join(dirname, prog)
    if os.path.exists(path + exe):
        return path


def cbin(prog):
    """$CCP4/bin unless prog or prog.exe is found in the dimple directory.
    Return value: path with filename without extension.
    """
    assert os.environ.get("CCP4")
    return check_prog(_dimple_dir, prog) or \
            os.path.join(os.environ["CCP4"], "bin", prog)


def syspath(prog):
    """Search prog(.exe) in the dimple directory and in the system $PATH.
    Return value: path with filename without extension.
    """
    dirs = [_dimple_dir] + os.environ["PATH"].split(os.pathsep)
    for d in dirs:
        path = check_prog(d, prog)
        if path:
            return path
    put_error("Program not found: %s" % prog)


# Use relpath if possible, absolute paths clutter commands and make
# moving directory harder.
def adjust_path(path, relative_to):
    if os.path.isabs(path) or os.path.isabs(relative_to):
        return os.path.abspath(path)
    elif os.path.realpath(relative_to) != os.path.abspath(relative_to):
        # symlink in relative_to could make mess
        return os.path.abspath(path)
    else:
        return os.path.relpath(path, relative_to)

def _find_mount_point(path):
    path = os.path.abspath(path)
    while not os.path.ismount(path):
        path = os.path.dirname(path)
    return path


def _report_quota(quota_prog, mount_point):
    try:
        out = subprocess.check_output([quota_prog, '-w', '-p', '-f',
                                       mount_point], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        return
    lines = out.splitlines()
    if len(lines) == 3:
        if lines[1].split()[1:3] == ['blocks', 'quota']:
            blocks, quota = lines[2].split()[1:3]
            try:
                percent = '%.0f%%' % (100. * int(blocks) / int(quota))
            except ValueError:
                percent = '???'
            comment('\nUsed quota on %s: %s / %s kB (%s)' %
                    (mount_point, blocks, quota, percent))
    else:
        try:
            out = subprocess.check_output([quota_prog],
                                          stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError:
            return
        if out:
            _log_comment(out)


def report_disk_space(paths):
    if os.name == 'nt':
        return
    previous_mount_points = []
    for path in paths:
        mount = _find_mount_point(path)
        if mount in previous_mount_points:
            continue
        previous_mount_points.append(mount)
        try:
            s = os.statvfs(mount)
        except AttributeError:
            return
        free = s.f_frsize * s.f_bavail
        comment('\nFree space on %s: %.0f MB' % (mount, free / (1024.*1024)))
        for d in os.environ['PATH'].split(os.pathsep):
            quota = os.path.join(d, 'quota')
            if os.path.exists(quota):
                _report_quota(quota, mount)
                break


if __name__ == '__main__':
    print '--- testing report_disk_space() ---'
    path_args = sys.argv[1:] if len(sys.argv) > 1 else ['.']
    report_disk_space(path_args)
    print
