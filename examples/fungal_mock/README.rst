.. _fungal_mock_sample_data:

Introduction
============

Sample data
-----------

Based on the two amplicon sequencing libraries from this paper:

    Bakker (2018) A fungal mock community control for amplicon sequencing experiments
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


Setup
-----

Either clone the THAPBI PICT source code repository, on decompress the
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

Script ``setup.py`` will download the raw gzip compressed FASTQ files from
the ENA (122 files, 61 pairs, under 400MB in total), and setup appropriate
per-sample symlinks to the expected output in the ``expected/``
sub-directories for use with classifier assessment.


Setup
-----

We assume you have aquired the THAPBI PICT source code, and have your command
line terminal open in the ``examples/fungal_mock/`` folder. First we run the
``setup.py`` script:

.. code:: console

   $ ./setup.py

This will download the raw gzip compressed FASTQ files from the ENA (124 files,
61 pairs, under 400MB in total), and setup appropriate per-sample symlinks to
the expected output in the ``expected/`` sub-directories for use with classifier
assessment.

If you have the ``md5sum`` tool installed (standard on Linux), verify the FASTQ
files downloaded correctly:

.. code:: console

   $ cd amp_lib_one/raw_data/
   $ md5sum -c MD5SUM.txt
   $ cd ../../

.. code:: console

   $ cd	amp_lib_two/raw_data/
   $ md5sum -c MD5SUM.txt
   $ cd	../../

Running the pipeline
--------------------

Next, you can run the ``run.py`` script which will call THAPBI PICT multiple times.
For each of the primer settings (a small fragment of ITS1 on ``amp_lib_one/``,
a larger fragment of ITS1 and an ITS2 marker on ``amp_lib_two/``), it will make a
simple database using the provided ``ITS1.fasta`` or ``ITS2.fasta`` file. It will
then call the pipeline using the ``identity``, ``onebp`` and ``blast`` classifiers.
