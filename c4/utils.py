import os
import sys
import time

_logfile = None
_logfile_sections = None
_c4_dir = os.path.abspath(os.path.dirname(__file__))

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
    log_value("cwd", os.getcwd())
    log_value("prog", sys.argv[0])
    if len(sys.argv) > 1:
        _logfile.write("args:\n")
        for arg in sys.argv[1:]:
            _logfile.write(" %s\n" % arg)
    log_value("output_dir", output_dir)
    log_value("CCP4", os.getenv("CCP4", ""))
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
        value = str(value).rstrip().replace("\n", "\n ")
        _logfile.write("%s: %s\n" % (key, value))

def log_time(key, timestamp):
    log_value(key, time.strftime("%Y-%m-%d %H:%M:%S",
                                 time.localtime(timestamp)))

def put(text):
    sys.stdout.write(text)

def put_green(text):
    if hasattr(sys.stdout, 'isatty') and sys.stdout.isatty():
        if os.name != 'nt':
            put("\033[92m%s\033[0m" % text)
        else:
            put(text)
    else:
        put(text)

def comment(text):
    put(text)
    sys.stdout.flush()
    _log_comment(text)

def reset_color():
    if hasattr(sys.stdout, 'isatty') and sys.stdout.isatty():
        if os.name != 'nt':
            put("\033[0m")

def put_error(err, comment=None):  # pylint: disable=redefined-outer-name
    _log_comment("Error: %s" % err)
    if hasattr(sys.stderr, 'isatty') and sys.stderr.isatty():
        if os.name != 'nt':
            err = "\033[91m%s\033[0m" % err  # in bold red
    sys.stdout.flush()
    sys.stderr.write("\nError: %s.\n" % err)
    if comment is not None:
        _log_comment(comment)
        sys.stderr.write(comment + "\n")


def check_prog(dirname, prog):
    """If prog(.exe) is in dirname return the path (without extension)."""
    exe = '.exe' if os.name == 'nt' else ''
    path = os.path.join(dirname, prog)
    if os.path.exists(path + exe):
        return path


def cbin(prog):
    """If prog/prog.exe is not in c4/ then $CCP4/bin is assumed.
    Return value: path with filename without extension.
    """
    assert os.environ.get("CCP4")
    return check_prog(_c4_dir, prog) or \
            os.path.join(os.environ["CCP4"], "bin", prog)


def syspath(prog):
    """If prog/prog.exe is not in c4/ then search in the system $PATH.
    Return value: path with filename without extension.
    """
    dirs = [_c4_dir] + os.environ["PATH"].split(os.pathsep)
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

