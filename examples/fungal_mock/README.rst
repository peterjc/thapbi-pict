.. _fungal_mock_sample_data:

Introduction
============

Data source
-----------

This example is based on the two amplicon sequencing libraries from this paper:

    Bakker (2018) A fungal mock community control for amplicon sequencing
    experiments.
    https://doi.org/10.1111/1755-0998.12760
    https://www.ebi.ac.uk/ena/data/view/PRJNA377530

Both amplicon sequencing libraries were sequenced on an Illumina MiSeq
(multiplexed with other unrelated samples), giving 61 FASTQ paired files.

Amplicon library one amplified a small region of ITS1 using the BITS/B58S3
primer pair (``ACCTGCGGARGGATC`` and ``GAGATCCRTTGYTRAAAGTT``), as shown in
the paper's supplementary Table S4.

Amplicon library two amplified a larger region of ITS1 using the ITS1f/ITS2
primer pair (``CTTGGTCATTTAGAGGAAGTAA`` and ``GCTGCGTTCTTCATCGATGC``), which
includes the first library’s target region entirely. Similar yields as per
supplementary Table S4 vs S5.

Additionally, amplicon library two amplified ITS2 using the ITS3‐KYO2 and
ITS4‐KYO3 primer pair (``GATGAAGAACGYAGYRAA`` and ``CTBTTVCCKCTTCACTCG``),
with lower yields as per supplementary Table S5 vs S6.

This means we need to run THAPBI PICT three times (for each primer pair used
in each amplicon library). In fact the example runs it four times, as we can
also try the BITS/B58S3 primers on the second amplicon library, because they
amplify a subregion of what the ITS1f/ITS2 pair amplify. See the primer
discussion on the similar Redekar *et al.* (2019) worked example.

Provided files
--------------

Either clone the THAPBI PICT source code repository, or decompress the
latest source code release (``.tar.gz`` file). You should find it contains
a directory ``examples/fungal_mock/`` which is for this example.

File ``PRJNA377530.txt`` was download from the ENA and includes the FASTQ
checksums, URLs, and sample metadata.

Files ``ITS1.fasta`` and ``ITS2.fasta`` were extracted from supplementary
materials appendix S2, with the species name alone added to the FASTA titles
(for use with ``thapbi_pict curated-import``).

The provided ``negative_control.known.tsv`` and ``mock_community.known.tsv``
files lists the expected species in the negative controls (none) and the mock
community samples (the same 19 species, although not always in equal ratios).

The two folders ``amp_lib_one/`` and ``amp_lib_two/`` will be used for the
two separate amplicon libraries discussed in the paper. Each comes with a
``metadata.tsv`` file based on the metadata downloaded from the ENA, with
some reformatting. The split into amplicon one and two was based on
supplementary Tables S4, S5 and S6 (for the mock community samples) and
reading the paper (for placing the negative controls). They each have a
``raw_data/`` subdirectory containing a file named ``MD5SUM.txt`` which
can be used to validate the FASTQ files.

Shell scripts ``setup.py`` and ``run.sh`` should reproduce the analysis
discussed in the THAPBI PICT documentation.

Setup
-----

We assume you have aquired the THAPBI PICT source code, and have your command
line terminal open in the ``examples/fungal_mock/`` folder. First we run the
``setup.py`` script:

.. code:: console

   $ ./setup.py

This will download the raw gzip compressed FASTQ files from the ENA (122 files,
61 pairs, under 400MB in total), and setup appropriate per-sample symlinks to
the expected output in the ``expected/`` sub-directories for use with
classifier assessment.

If you have the ``md5sum`` tool installed (standard on Linux), verify the FASTQ
files downloaded correctly:

.. code:: console

   $ cd amp_lib_one/raw_data/
   $ md5sum -c MD5SUM.txt
   $ cd ../../

.. code:: console

   $ cd amp_lib_two/raw_data/
   $ md5sum -c MD5SUM.txt
   $ cd ../../

There is no need to decompress the files.

Running the pipeline
--------------------

Next, you can run the ``run.py`` script which will call THAPBI PICT multiple
times. For each of the two primer settings (a small fragment of ITS1 on
``amp_lib_one/``, and a larger fragment of ITS1 and an ITS2 marker on
``amp_lib_two/``), it will make a simple database using the provided
``ITS1.fasta`` or ``ITS2.fasta`` file. It will then call the pipeline using
the ``identity``, ``onebp`` and ``blast`` classifiers.

Additionally for each primer set, it will generate an additional edit-graph at
a higher minimum abundance threshold (to omit PCR noise etc).

Metadata
--------

The amplicon specific files ``metadata.tsv`` have seven columns:

1. Accession, assigned by the public archive, e.g. "SRR5314337"
2. MiSeq-name, author's filename stem, e.g. "FMockE.HC_S190"
3. Condition, based on original name without replicate suffix, e.g. "MockE_HC"
4. Replicate, numeric, e.g. "1"
5. Sample-type, either "fungal mock community" or "negative control"
6. Group, e.g. "even" or "staggered A"
7. Protocol, e.g. "high PCR cycle number" or "standard workflow"

When calling THAPBI PICT, the meta data commands are given as follows:

.. code:: console

    $ thapbi_pict ... -t metadata.tsv -c 5,6,7,3,4,2 -x 1 -g 6

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
