"""Backward-compatible alias for :mod:`specproc.utils.identify`.

The interactive identification window used to live at
``specproc.gui.line_identifier``; it now lives with the other standalone
tools in :mod:`specproc.utils`.  Importing this module keeps old scripts
working, but new code should import the ``specproc.utils`` path (or run
``specproc identify``).
"""

import sys

from specproc.utils.identify import LineIdentifier, main  # noqa: F401

if __name__ == '__main__':  # pragma: no cover
    sys.exit(main())
