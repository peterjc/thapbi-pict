.. _sample_data:

Introduction
============

Sample data
-----------

This example is based on the following paper from another research group:

* Redekar et al (2019) Diversity of *Phytophthora*, *Pythium*, and
  *Phytopythium* species in recycled irrigation water in a container nursery.
  https://doi.org/10.1094/PBIOMES-10-18-0043-R

Importantly, they used a different but overlapping set of PCR primers. Their
left primer matched the THAPBI PICT default, but thier right primer would
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
dedicated to this analysis. Start by making three sub-folders as follows:

.. code:: console

   $ mkdir raw_data/ intermediate/ summary/

We will use the file ``PRJNA417859.txt`` (originally downloaded from the ENA)
to automated downloading the FASTA files, which together are about 5GB:

.. code:: console

    $ echo TODO
    $ ls -1 raw_data/*.fastq.gz | wc -l
    768

We will also need ``metadata.tsv`` (extracted from ``PRJNA417859.txt``)::

.. code:: console

    $ echo TODO

