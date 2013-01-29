import os

#if hasattr(sys.stdout, 'isatty') and sys.stdout.isatty():

def red(text):
    return ("\033[91m%s\033[0m" if os.name != 'nt' else "%s") % text
def green(text):
    return ("\033[92m%s\033[0m" if os.name != 'nt' else "%s") % text

