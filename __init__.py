
import sys
if sys.version_info[:2] != (2, 7) and sys.version_info[:2] < (3, 5):
    sys.stderr.write('Error. Python 2.7 or 3.5+ is required.\n')
    sys.exit(1)
