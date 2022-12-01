.. _custom_database_defaults:

Pipeline with defaults
======================

Running thapbi-pict pipeline
----------------------------

First, we will run the THAPBI PICT pipeline command with largely default
settings (including the default database and primers), other than including
the metadata about the water samples. Note that this dataset has no blanks or
negative controls, so we must trust the default minimum abundance threshold.

The key values which we will be changing later are the primers and database.

Assuming you have the FASTQ files in ``raw_data/``, run the pipeline command
as follows, and you should get the listed output report files:

.. code:: console

    $ mkdir -p intermediate_defaults/ summary/
    $ thapbi_pict pipeline \
      -i raw_data/ -o summary/recycled-water-defaults \
      -s intermediate_defaults/ \
      -t metadata.tsv -x 7 -c 1,2,3,4,5,6
    ...
    onebp classifier assigned species/genus to 431 of 794 unique sequences from 1 files
    Wrote summary/recycled-water-defaults.ITS1.samples.onebp.*
    Wrote summary/recycled-water-defaults.ITS1.reads.onebp.*
    ...
    $ ls -1 summary/recycled-water-defaults.*
    summary/recycled-water-defaults.ITS1.all_reads.fasta
    summary/recycled-water-defaults.ITS1.all_reads.onebp.tsv
    summary/recycled-water-defaults.ITS1.reads.onebp.tsv
    summary/recycled-water-defaults.ITS1.reads.onebp.xlsx
    summary/recycled-water-defaults.ITS1.samples.onebp.tsv
    summary/recycled-water-defaults.ITS1.samples.onebp.xlsx
    summary/recycled-water-defaults.ITS1.tally.tsv

Here we used ``-r`` (or ``--report``) to specify a different stem for the
report filenames. The :ref:`sample metadata options <metadata>` were described
earlier -- this is perhaps an idealised example in that ``metadata.tsv`` was
created so that we add the first six columns the table (sorted in that order),
where ``-x 7`` means index to the accession (filename prefix) in column seven.

Notice the output reported a taxonomic assignment for 431 of 794 unique
sequences - that's 54%, but considerably higher if we consider the reads.

Results
-------

We will compare and contrast the following four samples with the second run
using different primers and a custom database. These were deliberately picked
from the less diverse samples for clarity.

Here we pick out the four samples at the command line with ``grep``, you
can also look at the ``recycled-water-defaults.ITS1.samples.onebp.xlsx``
file in Excel:

.. code:: console

    $ cut -f 6,7,8 summary/recycled-water-defaults.ITS1.samples.onebp.tsv \
      | grep -E "(SRR6303586|SRR6303586|SRR6303588|SRR6303596|SRR6303948)"
    OSU482       SRR6303588  Phytophthora chlamydospora, Phytophthora x stagnum(*), Unknown
    OSU483       SRR6303586  Phytophthora chlamydospora, Phytophthora x stagnum(*)
    OSU536.s203  SRR6303948  Phytophthora ramorum
    OSU121       SRR6303596  Phytopythium (unknown species)

Three of these four have *Phytophthora* (and one with an unknown), while
the fourth has *Phytopythium*. However, this is discarding all the reads
which do not match the default *Phytophthora* centric primers.
