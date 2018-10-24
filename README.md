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

# Development Notes

## Python style conventions

The Python code follows [PEP8](https://www.python.org/dev/peps/pep-0008/)
and [PEP257 docstring](https://www.python.org/dev/peps/pep-0257/) style,
guided by the [Zen of Python](https://www.python.org/dev/peps/pep-0020/).

## Dependencies

We will use the [Conda](https://conda.io/) packaging system, including the
[BioConda](https://bioconda.github.io/) channel for bioinformatics packages,
and [PyPI](https://pypi.python.org/) packages to define a reproducible
environment for running this software.
