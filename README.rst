THAPBI *Phytophthora* ITS1 Classifier Tool (PICT)
=================================================

.. image:: https://img.shields.io/pypi/v/thapbi_pict.svg
   :alt: THAPBI PICT on the Python Package Index (PyPI)
   :target: https://pypi.org/project/thapbi-pict/
.. image:: https://img.shields.io/conda/vn/bioconda/thapbi-pict.svg
   :alt: THAPBI PICT on BioConda
   :target: https://anaconda.org/bioconda/thapbi-pict
.. image:: https://img.shields.io/travis/peterjc/thapbi-pict/master.svg?label=master&logo=travis
   :alt: THAPBI PICT TravisCI build status
   :target: https://travis-ci.org/peterjc/thapbi-pict/branches
.. image:: https://img.shields.io/circleci/project/github/peterjc/thapbi-pict/master.svg?label=master&logo=CircleCI
   :alt: THAPBI PICT CircleCI build status
   :target: https://circleci.com/gh/peterjc/thapbi-pict/tree/master
.. image:: https://readthedocs.org/projects/thapbi-pict/badge/?version=latest
   :alt: Documentation Status
   :target: https://thapbi-pict.readthedocs.io/en/latest/?badge=latest
.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :alt: Code style: black
   :target: https://github.com/python/black

THAPBI PICT an ITS1-based diagnostic/profiling tool from the UK
`BBSRC <https://www.bbsrc.ac.uk>`__ funded Tree Health and Plant
Biosecurity Initiative (THAPBI) Phyto-Threats project, focused on
identifying *Phytophthora* species present in Illumina sequenced
environmental samples.

*Phytophthora* (from Greek meaning plant-destroyer) species are
economically important plant pathogens, important in both agriculture
and forestry. ITS1 is short for Internal Transcribed Spacer one, which
is a region of eukaryotes genomes between the 18S and 5.8S rRNA genes.
This is commonly used for molecular barcoding, where sequencing this
short region can identify species.

This `software repository <https://github.com/peterjc/thapbi-pict/>`_
continues earlier work including:

- https://github.com/widdowquinn/THAPBI
- https://github.com/widdowquinn/THAPBI-pycits
- https://github.com/peterthorpe5/THAPBI-pycits
- https://github.com/peterthorpe5/public_scripts/tree/master/metapy
- https://github.com/peterthorpe5/public_scripts/tree/master/metapy_tools

The `tool documentation <https://thapbi-pict.readthedocs.io/e>`_ is
hosted by `Read The Docs <https://readthedocs.org/>`_, generated
automatically from the ``docs/`` folder of this repository.


Installation
============

We recommend installing this tool on Linux using the
`Conda <https://conda.io/>`__ packaging system, via the
`BioConda <https://bioconda.github.io/>`__ channel, which will handle
*all* the dependencies:

.. code:: bash

   $ conda install thapbi_pict

The same should work on macOS, but the command line dependencies likely
rule out using Windows directly.

Alternatively, since the software is on the `Python Package Index
(PyPI) <https://pypi.python.org/>`__, the following command will install
it along with its Python dependencies:

.. code:: bash

   $ pip install thapbi_pict

However, in this case you will still need to install various external
command line tools like ``hmmer``, and others which are only used for
some classifiers (like ``blast`` and ``swarm``). If you have BioConda
setup, use the following:

.. code:: bash

   $ conda install --file requirements-ext.txt

On a typical Linux system most of the tools required will be available
via the default distribution packages, although not always under the
same package name.

On Debian (with the efforts of DebianMed), or Ubuntu Linux, try:

.. code:: bash

   $ sudo apt-get install ncbi-blast+ cutadapt hmmer swarm trimmomatic

If you want to install the very latest unreleased code, you must
download the source code from GitHub - see the ``CONTRIBUTING.rst`` file
for more details.

Once installed, you should be able to run the tool using:

.. code:: bash

   $ thapbi_pict

This should automatically find the installed copy of the Python code.
Use ``thapbi_pict -v`` to report the version, or ``thapbi_pict -h`` for
help.


Release History
===============

See the ``CHANGELOG.rst`` file (it was getting too long to include here).


Development Notes
=================

Please see the ``CONTRIBUTING.rst`` file for details of the development
setup including Python style conventions, git pre-commit hook,
continuous integration and test coverage, and release process.
