[![THAPBI PICT TravisCI build status](https://api.travis-ci.org/peterjc/thapbi-pict.svg?branch=master)](https://travis-ci.org/peterjc/thapbi-pict/branches)
[![THAPBI PICT CircleCI build status](https://circleci.com/gh/peterjc/thapbi-pict/tree/master.svg?style=svg)](https://circleci.com/gh/peterjc/thapbi-pict/tree/master)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

# THAPBI *Phytophthora* ITS1 Classifier Tool (PICT)

*Phytophthora* (from Greek meaning plant-destroyer) species are economically
important plant pathogens, important in both agriculture and forestry. ITS1 is
short for Internal Transcribed Spacer one, which is a region of eukaryotes
genomes between the 18S and 5.8S rRNA genes. This is commonly used for
molecular barcoding, where sequencing this short region can identify species.

This repository is for development of ITS1-based diagnostic/profiling tools
for the Tree Health and Plant Biosecurity Initiative (THAPBI) Phyto-Threats
project, funded by the UK's [BBSRC](https://www.bbsrc.ac.uk).

This continues earlier work including:

* https://github.com/widdowquinn/THAPBI-pycits
* https://github.com/peterthorpe5/THAPBI-pycits
* https://github.com/peterthorpe5/public_scripts/tree/master/metapy
* https://github.com/peterthorpe5/public_scripts/tree/master/metapy_tools

# Installation

In future we intend to release this software on PyPI to allow you to install
with ``pip install thapbi_pict`` and on BioConda which would allow install
with just ``conda install thapbi_pict``.

However, for now you should download the source code from GitHub, decompress
it if required, and run ``pip3 install .`` which should automatically get
our Python dependencies.

Once installed, you should be able to run the tool using:

   thapbi_pict

This should automatically find the installed copy of the Python code.


# Release History

*Pending...*

# Development Notes

Please see the ``CONTRIBUTING.md`` file for details of the development
setup including Python style conventions, git pre-commit hook, continuous
integration and test coverage.
