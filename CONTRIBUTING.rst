Development Notes
=================

Python style conventions
------------------------

The Python code follows
`PEP8 <https://www.python.org/dev/peps/pep-0008/>`__ and `PEP257
docstring <https://www.python.org/dev/peps/pep-0257/>`__ style, guided
by the `Zen of Python <https://www.python.org/dev/peps/pep-0020/>`__.

Practically, coding style is enforced with the command line tools
`black <https://github.com/python/black>`__ (which can automatically
edit your code) and `flake8 <http://flake8.pycqa.org/>`__ (which in
addition to its own style checking, has a range of plugins - including a
plugin to call ``black`` from ``flake8``).

You can install these tools using:

.. code:: console

   $ pip install black flake8 flake8-black flake8-blind-except flake8-docstrings \
   flake8-rst-docstrings flake8-bugbear flake8-pie flake8-comprehensions restructuredtext-lint

You can run the checks using:

.. code:: console

   $ flake8 .

You can ask black to edit your files with:

.. code:: console

   $ black .

We recommend enabling the flake8 git pre-commit hook as follows:

.. code:: console

   $ flake8 --install-hook git
   $ git config --bool flake8.strict true

If your editor can be configured to run flake8 and/or black
automatically, even better. These checks are done as part of the
continuous integration when changes are made on GitHub.


Continuous Integration
----------------------

Currently setup to do automated testing under Linux on two free
continuous integration services, CircleCI (using Conda for dependencies)
and TravisCI (using apt-packages and PyPI).

* CircleCI: https://circleci.com/gh/peterjc/thapbi-pict/tree/master

* TravisCI: https://travis-ci.org/peterjc/thapbi-pict/branches


Dependencies
------------

For end users, we recommend installing using the
`Conda <https://conda.io/>`__ packaging system, via the
`BioConda <https://bioconda.github.io/>`__ channel, which will handle
*all* the dependencies:

https://anaconda.org/bioconda/thapbi-pict

For development (and this is reflected in the TravisCI configuration),
we recommend installing via ``pip``, which uses
`PyPI <https://pypi.python.org/>`__ for the Python dependencies. These
are declared in the ``setup.py`` script and so will automatically be
installed, but see also the ``requirements.txt`` file.

The two requirements files (``requirements.txt`` for Python dependencies
- and ``requirements-ext.txt`` for external command line bioinformatics
tools) can contain exact pinned dependency versions, allowing us to
define a more reproducible environment for running this software if
needed.

Installing from source
----------------------

First, download the code from GitHub and decompress it if required. The
best way to do this if you are likely to contribute any changes is at
the command line with ``git``.

.. code:: console

   $ git clone https://github.com/peterjc/thapbi-pict.git
   $ cd thapbi-pict

Then load the plain text SQL dump of the default database into SQLite3,
see ``database/README.rst`` for more information on this. Make it read
only to prevent accidental edits:

.. code:: console

   $ sqlite3 thapb_pict/ITS1_DB.sqlite < database/ITS1_DB.sql
   $ chmod a-w thapbi_pict/ITS1_DB.sqlite

Assuming your default Python is at least version 3.5, to install the
tool and automatically get our Python dependencies:

.. code:: console

   $ pip install .

If your system defaults to Python 2, try ``pip3 install .`` or
``python3 -m pip install .`` instead.

Once installed, you should be able to run the tool using:

.. code:: console

   $ thapbi_pict

This should automatically find the installed copy of the Python code.
Use ``thapbi_pict -v`` to report the version, or ``thapbi_pict -h`` for
help.

Release process
---------------

For a release, start from a clean git checkout (to reduce the chance of
bundling any stray local files despite a cautious ``MANIFEST.in``).

If the DB has changed, and this was not done locally, we must update it
using the plain text dump which is under version control:

.. code:: bash

   rm -rf thapbi_pict/ITS1_DB.sqlite
   sqlite3 thapbi_pict/ITS1_DB.sqlite < database/ITS1_DB.sql
   chmod a-w thapbi_pict/ITS1_DB.sqlite

If not, skip directly to:

.. code:: bash

   python setup.py sdist --formats=gztar && python setup.py bdist_wheel
   twine upload dist/thapbi_pict-X.Y.Z*
   git tag vX.Y.Z
   git push origin master --tags

The PyPI upload should trigger an automated pull request updating the 
`THAPBI PICT BioConda recipe
<https://github.com/bioconda/bioconda-recipes/blob/master/recipes/thapbi-pict/meta.yaml>`__
which will need reviewing (e.g. new dependencies) before it is merged.
