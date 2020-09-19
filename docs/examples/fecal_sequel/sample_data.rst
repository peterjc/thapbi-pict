.. _fecal_sequel_sample_data:

Introduction
============

Data source
-----------

Script ``setup.sh`` will download the raw FASTQ files for Walker *et al.*
(2019) from https://www.ebi.ac.uk/ena/data/view/PRJNA574765

We focus on bioproject PRJNA574765 which has 60 samples and covers the mock
communities. Additionally the paper describes PRJNA525109 (41 samples
comparing genetic efficacy vs traditional survey techniques), and PRJNA525407
(9 samples looking at bat species assemblages in archaeological sites in
Belize, with an expanded reference set).

The reference set of COI sequences is taken from Supplementary S2 in the
preceding paper:

    Walker *et al.* (2016)
    Species From Feces: Order-Wide Identification of Chiroptera From Guano and
    Other Non-Invasive Genetic Samples.
    https://doi.org/10.1371/journal.pone.0162342

That paper included bioproject PRJNA325503 with 9 samples.

Provided files
--------------

Either clone the THAPBI PICT source code repository, or decompress the
latest source code release (``.tar.gz`` file). You should find it contains
a directory ``examples/fecal_sequel/`` which is for this example.

File ``PRJNA574765`` was download from the ENA and includes the FASTQ
checksums, URLs, and the key metadata. Related file ``metadata.tsv``
contains report-ready metadata about the samples (see below).

The subdirectory ``raw_data/`` will hold the compressed FASTA files. This
contains a file named ``MD5SUM.txt`` which can be used to validate the
FASTQ files (using the checksums provided by the ENA).

File ``mock_community.known.tsv`` describes the three species of bats expected
in the mock communities (which use different ratios).

File ``COI_430_bats.fasta`` of pre-trimmed bat COI markers will be generated
by downloading the FASTA file from Walker *et al.* (2016) Supplementary S2,
and underscores replaced with spaces in the record names.

File ``observed_3_bats.fasta`` contains alternative COI markers observed
in at least 10 samples, and their assumed species source.

Shell scripts ``setup.sh`` and ``run.sh`` should reproduce the analysis
discussed in the THAPBI PICT documentation.

Setup
-----

We assume you have acquired the THAPBI PICT source code, and have your command
line terminal open in the ``examples/fecal_sequel/`` folder. First we run
the ``setup.sh`` script:

.. code:: console

   $ ./setup.sh

This will download the reference FASTA file, and the the raw gzip compressed
FASTQ files from the ENA (120 files, 60 pairs, about 750MB in total), and
setup per-sample symlinks to the expected output in the ``expected/``
directory for use with classifier assessment.

If you have the ``md5sum`` tool installed (standard on Linux), verify the FASTQ
files downloaded correctly:

.. code:: console

    $ cd raw_data/
    $ md5sum -c MD5SUM.txt
    $ cd ..

There is no need to decompress the files.

Running the pipeline
--------------------

Next, you can run the ``run.sh`` script which will call THAPBI PICT to turn
the FASTA file into a database, and then run the pipeline.

Metadata
--------

The provided file ``metadata.tsv`` is based on ``PRJNA574765`` but breaks up
the sample name into separate columns:

1. Accession, assigned by the public archive, e.g. "SRR10198789"
2. Rare, which of the 3 species is at low abundance, "COTO", "EPFU" or "TABR".
3. Ratio, either "1:64" (rare) or "1:192" (very rare)
4. Replicate, "01" to "10" (leading zero for alphabetical sorting)

The four letter appreviations are *Corynorhinus townsendii* (COTO),
*Eptesicus fuscus* (EPFU) and *Tadarida brasiliensis* (TABR).

When calling THAPBI PICT, the meta data commands are given as follows:

.. code:: console

    $ thapbi_pict ... -t metadata.tsv -x 1 -c 2,3,4

Argument ``-t metadata.tsv`` says to use this file for the metadata.

The ``-x 1`` argument indicates the filename stem can be found in column 1,
Accession.

Argument ``-c 2,3,4`` says which columns to display and sort by (do not
include the indexed column again). i.e. Rare species, ratio, replicate.

We have not given a ``-g`` argument to assign colour bands in the Excel
reports, so it will default to the first column in ``-c``, meaning we get
three coloured bands for "COTO", "EPFU" and "TABR".
