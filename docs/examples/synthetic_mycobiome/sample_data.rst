.. _synthetic_mycobiome_sample_data:

Marker data
===========

Either clone the THAPBI PICT source code repository, or decompress the latest
source code release (``.tar.gz`` file). You should find it contains a
directory ``examples/synthetic_mycobiome/`` which is for this example.

Shell scripts ``setup.sh`` and ``run.sh`` should reproduce the analysis
discussed.

FASTQ data
----------

File ``PRJNA305924.tsv`` was download from the ENA and includes the FASTQ
checksums, URLs, and sample metadata (not just for the files we will be using,
but additional Illumina MiSeq runs, and Ion Torrent data too).

Script ``setup.sh`` will download the raw FASTQ files for two of the Illumina
MiSeq runs described in Palmer *et al.* (2018)
from https://www.ebi.ac.uk/ena/data/view/PRJNA305924

It will download 82 raw FASTQ files (41 pairs), taking about 6.0 GB on disk.

If you have the ``md5sum`` tool installed (standard on Linux), verify the
FASTQ files downloaded correctly:

.. code:: console

    $ cd raw_data/m4A
    $ md5sum -c MD5SUM.txt
    ...
    $ cd ../m6
    $ md5sum -c MD5SUM.txt
    ...
    $ cd ../..

There is no need to decompress the files.

Amplicon primers & reference sequences
--------------------------------------

A region of ITS2 was amplified using the fITS7/ITS4 primer pair
(``GTGARTCATCGAATCTTTG`` and ``TCCTCCGCTTATTGATATGC``) with an average
product length of 264bp using public fungal sequences.

The file ``references.fasta`` we provide is based on ``amptk_mock2.fa`` and
``amptk_mock3.fa`` from the authors' GitHub repository
<https://github.com/nextgenusfs/amptk/tree/master/amptk/DB>, but formatted
suitable for direct import into our tool with primer-trimmed sequences.

Additional file ``environment.fasta`` contains selected close matches to
sequences from the environmental samples in the NCBI found with BLASTN
against the NT database.

Metadata
--------

File ``metadata.tsv`` is based on the ENA metadata and the paper text. It has
four columns:

1. run_accession, assigned by the public archive, e.g. "SRR7109326"
2. library_name, with sequencing run as a prefix, e.g. "m6-stds" or "m6-301-1"
3. plate_name, the sequencing run, one of "m4A" or "m6"
4. sample_alias, as used in the paper, e.g. "BioMockStds" or "301-1"
5. group, human readable sample type, e.g. "Biological Mock" or "Environment"

When calling THAPBI PICT, the meta data commands are given as follows:

.. code:: console

    $ thapbi_pict ... -t metadata.tsv -x 1 -c 3,4,5

Argument ``-t metadata.tsv`` says to use this file for the metadata.

Argument ``-c 3,4,5`` says which columns to display and sort by. This means
plate name, sample alias, then group. This sorts first by the sequencing run.
Column 2, library name, is omitted as this is larglely redundant.

Argument ``-x 1`` (default, so not needed) indicates the filename stem can be
found in column 1, run accession.

Other files
-----------

Provided files ``BioMockStds.known.tsv``, ``BioMock.known.tsv``, and
``SynMock.known.tsv`` list the expected 25 species, 22 species, and 12
synthetic controls expected in the mock samples. Folder ``expected/`` is
created linking accession names to the appropriate species for assessing the
classifier performance.

Sub-folder under ``intermediate*/ITS2/`` are used for intermediate files,
in general there is a sub-folder for each primer-pair.
