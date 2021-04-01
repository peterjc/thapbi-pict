.. image:: https://img.shields.io/pypi/v/thapbi_pict.svg
   :alt: THAPBI PICT on the Python Package Index (PyPI)
   :target: https://pypi.org/project/thapbi-pict/
.. image:: https://img.shields.io/conda/vn/bioconda/thapbi-pict.svg
   :alt: THAPBI PICT on BioConda
   :target: https://anaconda.org/bioconda/thapbi-pict
.. image:: https://img.shields.io/circleci/project/github/peterjc/thapbi-pict/master.svg?label=master&logo=CircleCI
   :alt: THAPBI PICT CircleCI build status
   :target: https://circleci.com/gh/peterjc/thapbi-pict/tree/master
.. image:: https://img.shields.io/readthedocs/thapbi-pict.svg?logo=read-the-docs
   :alt: Documentation Status
   :target: https://readthedocs.org/projects/thapbi-pict/builds/
.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :alt: Code style: black
   :target: https://github.com/python/black
.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4529395.svg
   :alt: Zenodo DOI
   :target: https://doi.org/10.5281/zenodo.4529395

THAPBI *Phytophthora* ITS1 Classifier Tool (PICT)
=================================================

About
-----

THAPBI PICT is a sequence based diagnostic/profiling tool from the UK funded
Tree Health and Plant Biosecurity Initiative (THAPBI) `Phyto-Threats project
<https://www.forestresearch.gov.uk/research/global-threats-from-phytophthora-spp/>`_,
focused on identifying *Phytophthora* species present in Illumina sequenced
environmental samples.

*Phytophthora* (from Greek meaning plant-destroyer) species are economically
important plant pathogens, important in both agriculture and forestry. ITS1 is
short for Internal Transcribed Spacer one, which is a region of eukaryotes
genomes between the 18S and 5.8S rRNA genes. This is commonly used for
molecular barcoding, where sequencing this short region can identify species.

With appropriate primer settings and a custom database of full length markers,
THAPBI PICT can be applied to other organisms and/or barcode marker sequences
- not just *Phytophthora* ITS1. It requires overlapping paired-end Illumina
reads which can be merged to cover the *full* amplicon marker. Longer markers
or fragmented amplicons are not supported. Internally it works by tracking
amplicon sequence variants (ASVs), using MD5 checksums as identifiers.

The worked examples include oomycetes, fungi, microalgae, and bats, and cover
markers in ITS1, ITS2, 18S rRNA and COI and more. The main criteria has been
mock communities with known species composition.

Installation
------------

We recommend installing this tool on Linux or macOS using the
`Conda <https://conda.io/>`__ packaging system, via the
`BioConda <https://bioconda.github.io/>`__ channel, which will handle
*all* the dependencies:

.. code:: console

   $ conda install thapbi_pict

Sadly the command line dependencies likely rule out using Windows directly.

Alternatively, since `the software is on the Python Package Index (PyPI)
<https://pypi.org/project/thapbi-pict/>`__, the following command will install
it along with its Python dependencies:

.. code:: console

   $ pip install thapbi_pict

However, in this case you will still need to install various external
command line tools. See ``INSTALL.rst`` for more details, and if you
want to modify the software read ``CONTRIBUTING.rst`` as well.

Quick Start
-----------

Once installed, you should be able to run the tool at the command line
using:

.. code:: console

   $ thapbi_pict

This should automatically find the installed copy of the Python code.
Use ``thapbi_pict -v`` to report the version, or ``thapbi_pict -h`` for
help.

Documentation
-------------

The `tool documentation <https://thapbi-pict.readthedocs.io/>`_ is hosted by
`Read The Docs <https://readthedocs.org/>`_, generated automatically from the
``docs/`` folder.

The documentation includes more detailed discussion of the sample datasets
in the ``examples/`` folder (which are based on published datasets).

Funding
-------

The Phyto-Threats project was supported from 2016 to 2019 under the Tree
Health and Plant Biosecurity Initiative (THAPBI), jointly funded by the
Biotechnology and Biological Sciences Research Council (`BBSRC
<https://bbsrc.ukri.org/>`_), Department for Environment, Food and Rural
affairs (`DEFRA <https://www.gov.uk/government/organisations/department-for-environment-food-rural-affairs>`_),
Economic and Social Research Council (`ESRC <https://esrc.ukri.org>`_),
`Forestry Commission <https://www.gov.uk/government/organisations/forestry-commission>`_,
Natural Environment Research Council (`NERC <https://nerc.ukri.org>`_)
and `Scottish Government <https://www.gov.scot/>`_.

Background
----------

THAPBI PICT continues earlier work including:

- https://github.com/widdowquinn/THAPBI
- https://github.com/widdowquinn/THAPBI-pycits
- https://github.com/peterthorpe5/THAPBI-pycits
- https://github.com/peterthorpe5/public_scripts/tree/master/metapy
- https://github.com/peterthorpe5/public_scripts/tree/master/metapy_tools

Release History
---------------

See the ``CHANGELOG.rst`` file.

Development Notes
-----------------

Please see the ``CONTRIBUTING.rst`` file for details of the development
setup including Python style conventions, git pre-commit hook,
continuous integration and test coverage, and release process.
