.. _great_lakes_sample_data:

Marker data
===========

Either clone the THAPBI PICT source code repository, or decompress the
latest source code release (``.tar.gz`` file). You should find it contains
a directory ``examples/great_lakes/`` which is for this example.

Shell scripts ``setup.sh`` and ``run.sh`` should reproduce the analysis
discussed.

Subdirectories ``MOL16S/`` and ``SPH16S/`` are used for the different
amplicons (with different primer settings and reference databases).

FASTQ data
----------

File ``PRJNA379165.tsv`` was download from the ENA and includes the FASTQ
checksums, URLs, and sample metadata. Derived file ``metadata.tsv`` contains
report-ready metadata about the samples (see below).

Script ``setup.sh`` will download the raw FASTQ files for Klymus *et al.*
(2017) from https://www.ebi.ac.uk/ena/data/view/PRJNA379165

It will download 36 raw FASTQ files (18 pairs), taking 1.8GB on disk.

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

The ``MOL16S`` amplicon targeted a short fragment of the mtDNA 16S RNA gene
using degenerate primer pair MOL16S_F/MOL16S_R (``RRWRGACRAGAAGACCCT`` and
``ARTCCAACATCGAGGT``).

The ``SPH16S`` amplicon targeted sphaeriid mussel species where it amplified
an overlapping slightly downstream region of the mtDNA 16S RNA gene using
non-degenerate primers SPH16S_F/SPH16S_R (``TAGGGGAAGGTATGAATGGTTTG`` and
``ACATCGAGGTCGCAACC``).

This means we need to run THAPBI PICT twice (once for each primer pair,
against a different marker database each time).

Metadata
--------

The provided file ``metadata.tsv`` is based on metadata in the ENA, split into
separate columns for reporting. It has five columns:

1. Run accession, e.g. "SRR5534972"
2. Library name, e.g. "SC3PRO2"
3. Sample title, e.g. "Mock Community 2 MOL16S with Fish Block Primer"
4. Marker, "MOL16S" or "SPH16S"
5. Group, "Mock Community", "Aquarium", "River" or "Control"

When calling THAPBI PICT, the meta data commands are given as follows:

.. code:: console

    $ thapbi_pict ... -t metadata.tsv -x 1 -c 4,5,3,2

Argument ``-t metadata.tsv`` says to use this file for the metadata.

Argument ``-c 4,5,3,2`` says which columns to display and sort by. This means
Marker, Group, Sample Title, Library name. This splits up the samples first by
the expected marker, and then the group.

Argument ``-x 1`` the filename stems can be found in that column one.

Other files
-----------

Files ``MOL16S.fasta`` and ``SPH16S.fasta`` are for building reference
databases. These were constructed from the accessions in the paper listed in
Table 1, Table 8, Supplementary Table 1, Supplementary Table 3, and some
additional accessions for the mock community. The sequences were primer
trimmed using ``cutadapt`` (requiring both the left and right primer to be
present), and the description given cut to just species level (discarding
strain or isolate information).
