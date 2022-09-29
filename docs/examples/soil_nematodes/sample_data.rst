.. _soil_nematodes_sample_data:

Marker data
===========

Either clone the THAPBI PICT source code repository, or decompress the latest
source code release (``.tar.gz`` file). You should find it contains a
directory ``examples/soil_nematodes/`` which is for this example.

Shell scripts ``setup.sh`` and ``run.sh`` should reproduce the analysis
discussed.

FASTQ data
----------

File ``PRJEB27581.tsv`` was download from the ENA and includes the FASTQ
checksums, URLs, and sample metadata.

Script ``setup.sh`` will download the raw FASTQ files for Ahmed *et al.* (2019)
from https://www.ebi.ac.uk/ena/data/view/PRJEB27581

It will download 32 raw FASTQ files (16 pairs), taking 12GB on disk.

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

There were four separate markers used here, as shown in the paper's Table 2
together with the shared Illumina adaptors used.

The authors do not provide copies of their reference sequence databases with
the paper. Instead, files ``NF1-18Sr2b.fasta``, ``SSUF04-SSUR22.fasta``,
``D3Af-D3Br.fasta`` and ``JB3-JB5GED.fasta`` were based on the accessions
listed in the paper and close matches in the NCBI found with BLASTN against
the NT database. Note many of the species names have been reduced to just
"Genus sp." in line with the mock community entries, and all the fungal
entries are listed as just "Fungi".

Metadata
--------

File ``metadata.tsv`` is based on the ENA metadata and the paper text. It has
four columns:

1. run_accession, assigned by the public archive, e.g. "ERR2678656"
2. read_count, the number of paired reads in the raw FASTQ files.
3. sample, one of "MC1", "MC2", "MC3" for the mock communities, or "Blank"
4. marker, one of "NF1-18Sr2b", "SSUF04-SSUR22", "D3Af-D3Br" or "JB3-JB5GED"

When calling THAPBI PICT, the meta data commands are given as follows:

.. code:: console

    $ thapbi_pict ... -t metadata.tsv -x 1 -c 4,3

Argument ``-t metadata.tsv`` says to use this file for the metadata.

Argument ``-c 4,3`` says which columns to display and sort by. This means
sample and then marker. The purpose here is to group the samples logically
(sorting on accession would not work), and suitable for group colouring.

Argument ``-x 1`` (default, so not needed) indicates the filename stem can be
found in column 1, run accession.

Other files
-----------

The provided ``negative_control.known.tsv`` and ``mock_community.known.tsv``
files lists the expected species in the negative controls (none) and the mock
community samples (the same 23 species). Sub-folders under ``expected/`` are
created for each primer-pair, linking each accession name to either file as
appropriate for assessing the classifier performance.

Sub-folders under ``intermediate/`` are used for intermediate files, a folder
for each primer-pair.
