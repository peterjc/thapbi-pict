.. _fungal_mock_sample_data:

Marker data
===========

Either clone the THAPBI PICT source code repository, or decompress the
latest source code release (``.tar.gz`` file). You should find it contains
a directory ``examples/fungal_mock/`` which is for this example.

Shell scripts ``setup.sh`` and ``run.sh`` should reproduce the analysis
discussed.

FASTQ data
----------

File ``PRJNA377530.tsv`` was download from the ENA and includes the FASTQ
checksums, URLs, and sample metadata.

Script ``setup.sh`` will download the raw FASTQ files for Bakker (2018) from
https://www.ebi.ac.uk/ena/data/view/PRJNA377530

It will download 122 raw FASTQ files (61 pairs), taking 346MB on disk.

If you have the ``md5sum`` tool installed (standard on Linux), verify the FASTQ
files downloaded correctly:

.. code:: console

    $ cd raw_data/AL1/
    $ md5sum -c MD5SUM.txt
    ...
    $ cd ../../

.. code:: console

    $ cd raw_data/AL1/
    $ md5sum -c MD5SUM.txt
    ...
    $ cd ../../

There is no need to decompress the files.

Amplicon primers & reference sequences
--------------------------------------

Amplicon library one (AL1) amplified a small region of ITS1 using primer pair
BITS/B58S3 (``ACCTGCGGARGGATC`` and ``GAGATCCRTTGYTRAAAGTT``), as shown in the
paper's supplementary Table S4.

Amplicon library two (AL2) amplified a larger region of ITS1 using primer pair
ITS1f/ITS2 (``CTTGGTCATTTAGAGGAAGTAA`` and ``GCTGCGTTCTTCATCGATGC``), which
includes the first library’s target region entirely. Similar yields as per
supplementary Table S4 vs S5.

Additionally, amplicon library two (AL2) amplified ITS2 using primer pair
ITS3‐KYO2 with ITS4‐KYO3 (``GATGAAGAACGYAGYRAA`` and ``CTBTTVCCKCTTCACTCG``),
with lower yields as per supplementary Table S5 vs S6.

This means we need to run THAPBI PICT three times (for each primer pair used
in each amplicon library). In fact the example runs it four times, as we can
also try the BITS/B58S3 primers on the second amplicon library, because they
amplify a subregion of what the ITS1f/ITS2 pair amplify. See the primer
discussion on the similar Redekar *et al.* (2019) worked example.

Files ``ITS1.fasta`` and ``ITS2.fasta`` were extracted from supplementary
materials appendix S2, with the species name alone added to the FASTA titles
(for input to ``thapbi_pict curated-import`` with primer trimming).

Metadata
--------

The amplicon specific files ``metadata_AL1.tsv`` and ``metadata_AL2.tsv`` are
based on the metadata downloaded from the ENA, with some reformatting. The
split into amplicon one and two was based on supplementary Tables S4, S5 and
S6 (for the mock community samples) and reading the paper (for placing the
negative controls).

They have seven columns:

1. Accession, assigned by the public archive, e.g. "SRR5314337"
2. MiSeq-name, author's filename stem, e.g. "FMockE.HC_S190"
3. Condition, based on original name without replicate suffix, e.g. "MockE_HC"
4. Replicate, numeric, e.g. "1"
5. Sample-type, either "fungal mock community" or "negative control"
6. Group, e.g. "even" or "staggered A"
7. Protocol, e.g. "high PCR cycle number" or "standard workflow"

When calling THAPBI PICT, the meta data commands are given as follows:

.. code:: console

    $ thapbi_pict ... -t metadata_AL1.tsv -c 5,6,7,3,4,2 -x 1 -g 6
    $ thapbi_pict ... -t metadata_AL2.tsv -c 5,6,7,3,4,2 -x 1 -g 6

Argument ``-t metadata.tsv`` says to use this file for the metadata.

Argument ``-c 5,6,7,3,4,2`` says which columns to display and sort by. This
means Sample-type, Group, Protocol, Condition, Replicate, MiSeq Name. The
purpose here is to group the samples logically (sorting on accession or MiSeq
Name would not work), and suitable for group colouring.

Argument ``-x 1`` (default, so not needed) indicates the filename stem can be
found in column 1, Accession. We might have downloaded the files and used the
author original names, in which case ``-x 2`` ought to work.

Argument ``-g 6`` means assign colour bands using column 6, Group. This is
used in the Excel reports.

Other files
-----------

The provided ``negative_control.known.tsv`` and ``mock_community.known.tsv``
files lists the expected species in the negative controls (none) and the mock
community samples (the same 19 species, although not always in equal ratios).

Sub-folders under ``intermediate/`` are used for intermediate files, a folder
for each amplicon library (AL1 and AL2) and primer-pair combination.
