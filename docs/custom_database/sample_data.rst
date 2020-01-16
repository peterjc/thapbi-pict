.. _custom_database_sample_data:

Introduction
============

Sample data
-----------

This example is based on the following paper from another research group:

* Redekar *et al* (2019) Diversity of *Phytophthora*, *Pythium*, and
  *Phytopythium* species in recycled irrigation water in a container nursery.
  https://doi.org/10.1094/PBIOMES-10-18-0043-R

The THAPBI PICT default primer values reflect the Phyto-Threats project's
choice of primers where *Phytophthora* are the main target, but here the
target is wider covering *Phytophthora*, *Pythium*, and *Phytopythium*.
They therefore used a slightly different set of PCR primers:

    *PCR was performed using the primers internal transcribed spacer
    (ITS)6 (5′* ``GAAGGTGAAGTCGTAACAAGG`` *3′ *) that binds the small subunit
    rDNA upstream of the ITS1 region and ITS7 (5′* ``AGCGTTCTTCATCGATGTGC`` *3′)
    that binds within the 5.8S rDNA (Cooke et al. 2000)*

The left primer (ITS6) matches the THAPBI PICT default, but their right primer
(ITS7) matches about 60bp further downstream in *Phytophthora*. This means we
can use THAPBI PICT default settings and get meaningful but blinkered results
(for the subset of the data which our narrower primer set would have amplified,
using a *Phytophthora* centric database).

In order to classify beyound *Phytophthora*, we need to build a THABPI PICT
database including *Pythium* and *Phytopythium*. Redekar *et al* (2019)
Supplementary Table 3 provides a list of 1454 unique accessions and the
species they assigned to it (not always the same as that listed on the NCBI
record, as those annotations can change). Looking at those sequences, bar
a handful they extend though the right primer. However, only about 50 have
the left primer sequence included (depending how stringent you are), and
the rest are also missing the next 32bp.

The ITS6 primer is situated within a highly conserved region, and the next
32bp is highly conserved, usually ``TTTCCGTAGGTGAACCTGCGGAAGGATCATTA``.
Unfortunately, the majority of published *Oomycetes* ITS1 sequences omit
this. For the curated *Phytophthora* in the THAPBI PICT default database,
we have inserted the expected sequence - and have yet to find a counter
example. However, Redekar *et al* took the other obvious choice, and
remove it from their reads:

    *trimming extra bases from read1: an additional 32 bases from the 5′ end
    of read1, which mapped to 18S segment, were trimmed as the oomycete ITS
    reference database does not include the 18S segment;*

We can do something similar in THAPBI PICT by treating this typically
conserved 32bp region as part of the left primer - requiring it be present
(while allowing some ambiguity) and removing it - leaving a shorter fragment
which can be matched to a database built of those 1454 accessions.

Raw data
--------

The raw data is 384 paired FASTQ files (samples to ``SAMN08012674`` to
``SAMN08013057``, runs ``SRR6303585`` to ``SRR6303968``),
available from the NCBI Sequence Read Archive under bioproject
`PRJNA417859 <https://www.ncbi.nlm.nih.gov/bioproject/PRJNA417859/>`_,
or equivalently from the European Nucleotide Archive under project
`PRJNA417859 <https://www.ebi.ac.uk/ena/data/view/PRJNA417859>`_. Downloaded
we have 768 gzip-compressed FASTQ files (384 pairs), taking about 5GB on disk.

The ENA makes it easy to download the project metadata as a table, a copy of
which we provide as file ``PRJNA417859.txt``. With a little scripting this
can be reformatted extracting the relevant :ref:`sample metadata <metadata>`
for use with THAPBI PICT.

Setup
-----

We assume you have your command line terminal open in a new empty folder
dedicated to this analysis. Start by fetching files ``PRJNA417859.txt``
(originally downloaded from the ENA) and ``metadata.tsv`` (generated from it,
ready for use in reporting). These are including with the THAPBI PICT source
code in the folder ``tests/recycled_water/`` but can be downloaded:

.. code:: console

    $ wget https://github.com/peterjc/thapbi-pict/raw/master/tests/recycled_water/PRJNA417859.txt
    $ wget https://github.com/peterjc/thapbi-pict/raw/master/tests/recycled_water/metadata.tsv

Next, we will download the gzip-compressed FASTQ files into a sub-directory
named ``raw_data/``. You may find the ENA bulk download helper application
easier, but the following should also work in principle:

.. code:: console

    $ mkdir raw_data
    $ cd raw_data/
    $ for URL in `cut -f 6 PRJNA417859.txt | sed "s/;/ /g"`; do wget "ftp://$URL"; done
    $ cd ..

If that worked you should have 768 files named ``SRR6303585_1.fastq.gz`` and
``SRR6303585_2.fastq.gz`` (the pair for ``SRR6303585`` aka ``OSU484``) though
to ``SRR6303968_1.fastq.gz`` and ``SRR6303968_2.fastq.gz`` (the pair for
``SRR6303968`` aka ``OSU476``):

.. code:: console

    $ ls -1 raw_data/SRR*.fastq.gz | wc -l
    768

At this point it is worth checking there were no partial downloads or data
corruption by validating the MD5 checksums given by the ENA:

.. code:: console

    $ cd raw_data/
    $ wget https://github.com/peterjc/thapbi-pict/raw/master/tests/recycled_water/MD5SUM.txt
    $ md5sum -c MD5SUM.txt
    $ cd ..

There is no need to decompress the files.
