.. _custom_database_sample_data:

Introduction
============

Sample data
-----------

This example is based on the following paper from another research group:

* Redekar et al (2019) Diversity of *Phytophthora*, *Pythium*, and
  *Phytopythium* species in recycled irrigation water in a container nursery.
  https://doi.org/10.1094/PBIOMES-10-18-0043-R

Importantly, they used a different but overlapping set of PCR primers. Their
left primer matches the THAPBI PICT default, but their right primer would
match about 60bp downstream in *Phytophthora*. This means we can use the
default settings and get meaningful but blinkered results (for the subset of
the data which the narrower primer set would have amplified).

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
