# Copyright 2018-2022 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""THAPBI *Phytophthora* ITS1 Classifier Tool (PICT).

You would typically use THAPBI PICT via the command line tool it defines::

    $ thapbi_pict --help
    ...

However, it is also possible to call functions etc from within Python. The
top level package currently only defines the tool version:

    >>> from thapbi_pict import __version__
    >>> print(__version__)
    ...

The `tool documentation <https://thapbi-pict.readthedocs.io/>`_ is
hosted by `Read The Docs <https://readthedocs.org/>`_, generated
automatically from the ``docs/`` folder of the `software repository
<https://github.com/peterjc/thapbi-pict/>`_ and the *"docstrings"*
within the source code which document the Python API.
"""

__version__ = "0.13.0"
