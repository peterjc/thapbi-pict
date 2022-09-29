.. _pest_insects_sample_data:

Marker data
===========

Either clone the THAPBI PICT source code repository, or decompress the latest
source code release (``.tar.gz`` file). You should find it contains a
directory ``examples/pest_insects/`` which is for this example.

Shell scripts ``setup.sh`` and ``run.sh`` should reproduce the analysis
discussed.

FASTQ data
----------

File ``PRJNA716058.tsv`` was download from the ENA and includes the FASTQ
checksums, URLs, and sample metadata.

Script ``setup.sh`` will download the raw FASTQ files for Batovska *et al.*
(2021) from https://www.ebi.ac.uk/ena/data/view/PRJNA716058

It will download 60 raw FASTQ files (30 pairs), taking 7.9 GB on disk.

If you have the ``md5sum`` tool installed (standard on Linux), verify the FASTQ
files downloaded correctly:

.. code:: console

    $ cd raw_data/
    $ md5sum -c MD5SUM.txt
    ...
    $ cd ../

There is no need to decompress the files.

Amplicon primers & reference sequences
--------------------------------------

Three separate markers used here, as shown in the paper's Supplementary Table
S2, together with the shared Illumina adaptors used.

The authors provide their reference species level sequences as a compressed
FASTA file ``merged_arthropoda_rdp_species.fa.gz`` on the GitHub repository
for the paper: https://github.com/alexpiper/HemipteraMetabarcodingMS

The worked example applies the three primer-pairs to this FASTA file to make
an amplicon specific FASTA file for each marker.

Metadata
--------

File ``metadata.tsv`` is based on the ENA metadata and the paper text. It has
four columns:

1. run_accession, assigned by the public archive, e.g. "SRR14022295"
2. sample_alias, e.g. "100-Pool-1" or "Trap-1"
3. source, e.g. one of the mock communities like "Pool 1", or "Trap"
4. individuals, e.g. "0100" (with leading zero for sorting) or "-" for traps.

When calling THAPBI PICT, the meta data commands are given as follows:

.. code:: console

    $ thapbi_pict ... -t metadata.tsv -x 1 -c 3,4,2

Argument ``-t metadata.tsv`` says to use this file for the metadata.

Argument ``-c 3,4,2`` says which columns to display and sort by. This means
by source (i.e. which mock community, or environmental traps), then number of
individuals in the mock, and finally the human readable sample alias.
The purpose here is to group the samples logically (sorting on sample_alias
would not work), and suitable for group colouring.

Argument ``-x 1`` (default, so not needed) indicates the filename stem can be
found in column 1, run accession.

Other files
-----------

Files ``mock_community_1.known.tsv``, ..., ``mock_community_5.known.tsv`` list
the expected species in the five different mock community pools. The setup
script will create symlinks using the sample names under sub-folder
``expected/`` pointing at the relevant community known file. This is for
automatically assessing the classifier performance.

Sub-folders under ``intermediate/`` are used for intermediate files, a folder
for each primer-pair.
