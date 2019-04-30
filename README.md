[![THAPBI PICT TravisCI build status](https://img.shields.io/travis/peterjc/thapbi-pict/master.svg?label=master&logo=travis)](https://travis-ci.org/peterjc/thapbi-pict/branches)
[![THAPBI PICT CircleCI build status](https://img.shields.io/circleci/project/github/peterjc/thapbi-pict/master.svg?label=master&logo=CircleCI)](https://circleci.com/gh/peterjc/thapbi-pict/tree/master)
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

* https://github.com/widdowquinn/THAPBI
* https://github.com/widdowquinn/THAPBI-pycits
* https://github.com/peterthorpe5/THAPBI-pycits
* https://github.com/peterthorpe5/public_scripts/tree/master/metapy
* https://github.com/peterthorpe5/public_scripts/tree/master/metapy_tools

# Installation

In the near future we intend to release this software as a BioConda package,
meaning the follwing command will install it along with all the dependencies:

```bash
$ conda install thapbi_pict
```

Until then, since the software is on the Python Package Index (PyPI), the
following command will install it along with its Python dependencies:

```bash
$ pip install thapbi_pict
```

However, in this case you will still need to install various external command
line tools like ``hmmer``, and others which are only used for some classifiers
(like ``blast`` and ``swarm``). If you have BioConda setup, use the following:

```bash
$ conda install blast cutadapt flash hmmer swarm trimmomatic
```

On a typical Linux system most of the tools listed will be available via the
default distribution packages, although not always under the same package name.

If you want to install the very latest unreleased code, you should download
the source code from GitHub, and decompress it if required. Then load the
plain text SQL dump of the default database into SQLite3 by running
``sqlite3 thapb_pict/ITS1_DB.sqlite < database/ITS1_DB.sql``, and next run
``pip3 install .`` which should automatically get our Python dependencies.

Once installed, you should be able to run the tool using:

```bash
$ thapbi_pict
```

This should automatically find the installed copy of the Python code. Use
``thapbi_pict -v`` to report the version, or ``thapbi_pict -h`` for help.


# Release History

| Version | Date       | Notes                                                                        |
|---------|------------|------------------------------------------------------------------------------|
| v0.0.1  | 2019-01-17 | Initial framework with ``identity`` and ``swarm`` classifiers.               |
| v0.0.2  | 2019-01-21 | Added ``assess`` command.                                                    |
| v0.0.3  | 2019-01-22 | Simplified generated filenames.                                              |
| v0.0.4  | 2019-01-24 | Added ``seq-import`` command, ``blast`` classifier, multi-taxon predictions. |
| v0.0.5  | 2019-02-06 | Hamming Loss in assessement output.                                          |
| v0.0.6  | 2019-02-07 | Misc. cleanup and import fixes.                                              |
| v0.0.7  | 2019-02-12 | Added ``plate-summary`` command, ``onebp`` classifier.                       |
| v0.0.8  | 2019-02-21 | Fix multi-class TN under-counting. New loss metric, ``swarmid`` classifier.  |
| v0.0.9  | 2019-03-05 | Looks for expected primers, discards mismatches. Caches HMM files locally.   |
| v0.0.10 | 2019-03-06 | Replace primer code allowing only 1bp differences with ``cutadapt``.         |
| v0.0.11 | 2019-03-08 | Speed up FASTQ preparation by using ``flash`` instead of ``pear`` v0.9.6.    |
| v0.0.12 | 2019-03-11 | Fixed bug in ``swarmid`` classifier.                                         |
| v0.0.13 | 2019-03-22 | Remove conserved 32bp when primer trim. Assess at sample level by default.   |
| v0.0.14 | 2019-04-01 | MD5 in dump output. Fixed importing sequences failing taxonomic validation.  |
| v0.0.15 | 2019-04-03 | Support for genus-level only entries in the database.                        |
| v0.1.0  | 2019-04-04 | Include a bundled ITS1 database.                                             |
| v0.1.1  | 2019-04-16 | Expand default taxonomy and database from Peronosporaceae to Peronosporales. |
| v0.1.2  | 2019-04-17 | Keep searching if ``onebp`` classifier perfect match is at genus-level only. |
| v0.1.3  | 2019-04-24 | Can optionally display sample metadata from TSV file in summary reports.     |
| v0.1.4  | 2019-04-25 | Sort samples using the optional metadata fields requested in reports.        |
| v0.1.5  | 2019-04-29 | Reworked optional metadata integration and its display in summary reports.   |
| v0.1.6  | 2019-04-30 | Include ready to use binary ITS1 database in source tar-ball & wheel files.  |


# Development Notes

Please see the ``CONTRIBUTING.md`` file for details of the development
setup including Python style conventions, git pre-commit hook, continuous
integration and test coverage.

For a release, start from a clean git checkout (to reduce the chance of
bundling any stray local files despite a cautious ``MANIFEST.in``).

```bash
rm -rf thapbi_pict/ITS1_DB.sqlite
sqlite3 thapbi_pict/ITS1_DB.sqlite < database/ITS1_DB.sql
chmod a-w thapbi_pict/ITS1_DB.sqlite
python setup.py sdist --formats=gztar
python setup.py bdist_wheel
twine upload dist/thapbi_pict-X.Y.Z.tar.gz dist/thapbi_pict-X.Y.Z-py3-none-any.whl
git tag -vX.Y.Z
git push origin master --tags
```
