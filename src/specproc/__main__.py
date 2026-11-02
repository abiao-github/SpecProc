"""Allow ``python -m specproc`` to behave like the ``specproc`` command."""

import sys

from specproc.main import main

if __name__ == '__main__':
    sys.exit(main())
