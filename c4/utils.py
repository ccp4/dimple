import os
import sys

_c4_dir = os.path.abspath(os.path.dirname(__file__))

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

def reset_color():
    if hasattr(sys.stdout, 'isatty') and sys.stdout.isatty():
        if os.name != 'nt':
            put("\033[0m")

def put_error(err, comment=None):
    if hasattr(sys.stderr, 'isatty') and sys.stderr.isatty():
        if os.name != 'nt':
            err = "\033[91m%s\033[0m" % err  # in bold red
    sys.stderr.write("Error: %s.\n" % err)
    if comment is not None:
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


