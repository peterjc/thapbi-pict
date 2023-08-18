.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.4529395.svg
   :alt: Zenodo DOI
   :target: https://doi.org/10.5281/zenodo.4529395
.. image:: https://img.shields.io/github/license/peterjc/thapbi-pict.svg?label=License
   :alt: MIT License
   :target: https://github.com/peterjc/thapbi-pict/blob/master/LICENSE.rst
.. image:: https://results.pre-commit.ci/badge/github/peterjc/thapbi-pict/master.svg
   :target: https://results.pre-commit.ci/latest/github/peterjc/thapbi-pict/master
   :alt: pre-commit.ci status
.. image:: https://img.shields.io/circleci/project/github/peterjc/thapbi-pict/master.svg?label=CI&logo=CircleCI
   :alt: CircleCI status
   :target: https://circleci.com/gh/peterjc/thapbi-pict/tree/master
.. image:: https://img.shields.io/appveyor/ci/peterjc/thapbi-pict/master.svg?logo=appveyor
   :alt: AppVeyor status
   :target: https://ci.appveyor.com/project/peterjc/thapbi-pict/history
.. image:: https://img.shields.io/readthedocs/thapbi-pict.svg?label=RTD&logo=read-the-docs
   :alt: Documentation Status
   :target: https://readthedocs.org/projects/thapbi-pict/builds/
.. image:: https://img.shields.io/pypi/v/thapbi_pict.svg?label=PyPI
   :alt: THAPBI PICT on the Python Package Index (PyPI)
   :target: https://pypi.org/project/thapbi-pict/
.. image:: https://img.shields.io/conda/vn/bioconda/thapbi-pict.svg?label=Bioconda
   :alt: THAPBI PICT on BioConda
   :target: https://anaconda.org/bioconda/thapbi-pict
.. image:: https://img.shields.io/badge/Code%20style-black-000000.svg
   :alt: Code style: black
   :target: https://github.com/python/black


THAPBI *Phytophthora* ITS1 Classifier Tool (PICT)
=================================================

About
-----

THAPBI PICT is a sequence based diagnostic/profiling tool from the UK funded
Tree Health and Plant Biosecurity Initiative (THAPBI) `Phyto-Threats project
<https://www.forestresearch.gov.uk/research/global-threats-from-phytophthora-spp/>`_,
initially focused on identifying *Phytophthora* species present in Illumina
sequenced environmental samples.

*Phytophthora* (from Greek meaning plant-destroyer) species are economically
important plant pathogens, in both agriculture and forestry. ITS1 is short for
Internal Transcribed Spacer one, which is a region of eukaryotes genomes
between the 18S and 5.8S rRNA genes. This is commonly used for molecular
barcoding, where sequencing this short region can identify species.

With appropriate primer settings and a custom database of full length markers,
THAPBI PICT can be applied to other organisms and/or barcode marker sequences
- not just *Phytophthora* ITS1. It requires overlapping paired-end Illumina
reads which can be merged to cover the *full* amplicon marker. Longer markers
or fragmented amplicons are not supported. Internally it works by tracking
unique amplicon sequence variants (ASVs), using MD5 checksums as identifiers.

The worked examples include oomycetes, fungi, fish, bats, and plants, and
cover markers in ITS1, ITS2, 12S, 16S, COI, and more. The main criteria has
been mock communities with known species composition.

Installation
------------

We recommend installing this tool on Linux or macOS using the
`Conda <https://conda.io/>`__ packaging system, via the
`BioConda <https://bioconda.github.io/>`__ channel, which will handle
*all* the dependencies:

.. code:: console

   $ conda install thapbi-pict

Alternatively or on Windows, since `the software is on the Python Package
Index (PyPI) <https://pypi.org/project/thapbi-pict/>`__, the following command
will install it along with its Python dependencies:

.. code:: console

   $ pip install thapbi-pict

However, in this case you will still need to install various external command
line tools. See ``INSTALL.rst`` for more details (especially for Windows),
and if you want to modify the software read ``CONTRIBUTING.rst`` as well.

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

Citation
--------

If you use THAPBI PICT in your work, please cite our *PeerJ* paper, and give
details of the version and any non-default settings used in your methods:

    Cock *et al.* (2023) "THAPBI PICT - a fast, cautious, and accurate
    metabarcoding analysis pipeline" *PeerJ* **11**:e15648
    https://doi.org/10.7717/peerj.15648

You can also cite the software specifically via Zenodo which offers version
specific DOIs as well as https://doi.org/10.5281/zenodo.4529395 which is for
the latest version.

Funding
-------

The initial work was supported from 2016 to 2019 under the Tree Health and
Plant Biosecurity Initiative (THAPBI) Phyto-Threats project:

  This research was supported by a grant funded jointly by the
  Biotechnology and Biological Sciences Research Council (`BBSRC
  <https://bbsrc.ukri.org/>`_), Department for Environment, Food and Rural
  affairs (`DEFRA <https://www.gov.uk/government/organisations/department-for-environment-food-rural-affairs>`_),
  Economic and Social Research Council (`ESRC <https://esrc.ukri.org>`_),
  `Forestry Commission <https://www.gov.uk/government/organisations/forestry-commission>`_,
  Natural Environment Research Council (`NERC <https://nerc.ukri.org>`_)
  and `Scottish Government <https://www.gov.scot/>`_, under the Tree
  Health and Plant Biosecurity Initiative, grant number ``BB/N023463/1``.

Work from 2020 to 2021 was supported in part under the *Early detection of
Phytophthora in EU and third country nurseries and traded plants (ID-PHYT)*
Euphresco project:

  Funded by DEFRA as part of the Future Proofing Plant Health project in
  support of Euphresco ID-PHYT.

Work from 2022 to 2027 was partly funded by the Rural & Environment Science
& Analytical Services (RESAS) Division of the Scottish Government.

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

See file ``CONTRIBUTING.rst`` for details of the development setup including
Python style conventions, git pre-commit hook, continuous integration and test
coverage, and release process.
