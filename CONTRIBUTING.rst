Development Notes
=================

Python style conventions
------------------------

The Python code follows `PEP8 <https://www.python.org/dev/peps/pep-0008/>`__
and `PEP257 docstring <https://www.python.org/dev/peps/pep-0257/>`__ style,
guided by the `Zen of Python <https://www.python.org/dev/peps/pep-0020/>`__.

Practically, coding style is enforced with several command line tools
including `ruff <https://github.com/astral-sh/ruff>`__ (which implements
a faster version of `black <https://github.com/python/black>`__) and `flake8
<http://flake8.pycqa.org/>`__ (with plugins) run via the tool `pre-commit
<https://pre-commit.com/>`__.

You can install these tools using:

.. code:: console

    $ pip install pre-commit
    $ pre-commit install  # within the thapbi_pict main directory

The checks will then run automatically when you make a git commit. You can
also run the checks directly using:

.. code:: console

    $ pre-commit run -a

If your editor can be configured to run flake8 and/or ruff automatically,
even better. These checks are done as part of the continuous integration when
changes are made on GitHub.


Continuous Integration
----------------------

Currently this is setup to do automated testing under Linux using free
continuous integration services:

* CircleCI (Linux): https://circleci.com/gh/peterjc/thapbi-pict/tree/master
* AppVeyor (Windows): https://ci.appveyor.com/project/peterjc/thapbi-pict/history

Dependencies
------------

See the main installation instructions for end users. For development we need
Python, a bash shell, git, and various other command line dependencies.
Installing THAPBI PICT from source (see below), will fetch Python dependencies.

The two requirements files (``requirements.txt`` for Python dependencies, and
``requirements-ext.txt`` for external command line bioinformatics tools) are
used in the continuous integration testing. These files can contain exact
pinned dependency versions, allowing us to define a more reproducible
environment for running this software if needed.

On Linux or macOS, you should have the bash shell and standard Unix tools like
``grep`` already installed. We recommend installing our specific command line
tool dependencies with  `Conda <https://conda.io/>`__ packaging system, via
the `BioConda <https://bioconda.github.io/>`__ channel:

.. code:: console

    $ conda install --file requirements-ext.txt

On Windows, few of the dependencies are available via Conda. The `Git For Windows
<https://gitforwindows.org>`_ installer will provide ``git``, ``bash``, ``grep``,
etc. You will also need to manually install sqlite3, flash, and NCBI BLAST.

Installing from source
----------------------

First, download the code from GitHub and decompress it if required. The best
way to do this if you are likely to contribute any changes is at the command
line with ``git``.

.. code:: console

    $ git clone https://github.com/peterjc/thapbi-pict.git
    $ cd thapbi-pict

Then build the default reference database, by loading the provided FASTA files
into SQLite3, see ``database/README.rst`` for more information on this. Make it
read only to prevent accidental edits:

.. code:: console

    $ cd database
    $ ./build_ITS1_DB.sh
    $ cd ..
    $ cp database/ITS1_DB.sqlite thapbi_pict/ITS1_DB.sqlite
    $ chmod a-w thapbi_pict/ITS1_DB.sqlite

Assuming your default Python is at least version 3.7, to install the tool and
automatically get our Python dependencies:

.. code:: console

    $ pip install .

If your system defaults to Python 2, try ``pip3 install .`` or
``python3 -m pip install .`` instead.

Once installed, you should be able to run the tool using:

.. code:: console

    $ thapbi_pict

This should automatically find the installed copy of the Python code.
Use ``thapbi_pict -v`` to report the version, or ``thapbi_pict -h`` for help.

Release process
---------------

For a release, start from a clean git checkout (to reduce the chance of
bundling any stray local files despite a cautious ``MANIFEST.in``). You will
need some python tools:

.. code:: console

    $ pip install -U pip twine build

First confirm if the DB at ``thapbi_pict/ITS1_DB.sqlite`` is up to date:

.. code:: bash

    sqlite3 thapbi_pict/ITS1_DB.sqlite .dump | grep -i "Imported with" | head -n 1

If there have been changes requiring the DB be rebuilt, do this:

.. code:: bash

    cd database
    ./build_ITS1_DB.sh
    git commit ITS1_DB.fasta -m "Rebuilt DB"
    cd ..

Next confirm the ``CHANGELOG.rst`` file is up to date, including using today's
date for the new version. Then actually do the build:

.. code:: bash

    rm -rf build/
    python -m build
    git tag vX.Y.Z
    git push origin master --tags
    twine upload dist/thapbi_pict-X.Y.Z*

The PyPI upload should trigger an automated pull request updating the
`THAPBI PICT BioConda recipe
<https://github.com/bioconda/bioconda-recipes/blob/master/recipes/thapbi-pict/meta.yaml>`__
which will need reviewing (e.g. new dependencies) before it is merged.

Must also turn the git tag into a "release" on GitHub, and attach the
wheel to it. This will generate a version specific DOI on Zenodo.
https://github.com/peterjc/thapbi-pict/releases
