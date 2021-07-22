Classifying sequences
=====================

Running thapbi-pict-classify
----------------------------

.. tip:

   If you don't have the FASTQ files, just the FASTA files, start from here.

The second stage of the pipeline is to merge all the sample specific FASTA
files into one non-redundant FASTA file, and then classify all the unique
sequences in it. These steps can be run separately:

.. code:: console

    $ thapbi_pict fasta-nr -h
    ...
    $ thapbi_pict classify -h
    ...

There are a number of options here, but for the purpose of this worked example
we will stick with the defaults and tell it to look for FASTA files in the
``intermediate/`` directory.

.. code:: console

    $ thapbi_pict fasta-nr -i intermediate/ -o summary/thapbi-pict.ITS1.all_reads.fasta
    ...
    $ thapbi_pict classify -i summary/thapbi-pict.ITS1.all_reads.fasta
    ...

Here we have not set the output folder with ``-o`` or ``--output``, which
means the classify step will default to writing the TSV output files next to
the input FASTA file. There should now be a new FASTA and TSV file:

.. code:: console

    $ ls -1 summary/thapbi-pict.ITS1.all_reads.*
    summary/thapbi-pict.ITS1.all_reads.fasta
    summary/thapbi-pict.ITS1.all_reads.onebp.tsv

Intermediate TSV files
----------------------

For each input FASTA file ``<name>.fasta`` a plain text tab separated variable
(TSV) file is generated named ``<name>.<method>.tsv`` where the default method
is ``onebp`` (which looks for perfect matches or up to one base pair
different). The first line is a header comment line (starting with ``#``)
labelling the columns, which are:

* Unique sequence name in ``<checksum>_<abundance>`` format.
* NCBI taxid of any predictions (semi-colon separated, same order as species)
* Genus-species of any predictions (semi-colon separated, alphabetical)
* Text note field in verbose mode (arbitrary debug text from the tool)

These files are not really intended for human use, but are readable:

.. code:: console

    $ head summary/thapbi-pict.ITS1.all_reads.onebp.tsv
    <SEE TABLE BELOW>

Viewing it like this is not ideal, although there are command line tools which
help. You could open the file in R, Excel, etc:

======================================= ============= ===============================================
#sequence-name                          taxid         genus-species
======================================= ============= ===============================================
2e4f0ed53888ed39a2aee6d6d8e02206_163094 221518        Phytophthora pseudosyringae
d9bc3879fdab3b4184c04bfbb5cf6afb_83653  631361        Phytophthora austrocedri
32159de6cbb6df37d084e31c37c30e7b_28976  67594         Phytophthora syringae
ed15fefb7a3655147115fc28a8d6d671_28786  78237         Phytophthora gonapodyides
972db44c016a166de86a2bacab3f4226_28515  53983;2056922 Phytophthora cambivora;Phytophthora x cambivora
c1a720b2005f101a9858107545726123_20400  78237         Phytophthora gonapodyides
96e0e2f0475bd1617a4b05e778bb04c9_17392  78237         Phytophthora gonapodyides
f27df8e8755049e831b1ea4521ad6eb3_16369  2496075;29920 Phytophthora aleatoria;Phytophthora cactorum
21d6308d89d74b8ed493d73a2cb4adb5_16169  53983         Phytophthora cambivora
======================================= ============= ===============================================

The first entry says the most abundance sequence with MD5 checksum
``2e4f0ed53888ed39a2aee6d6d8e02206`` was seen in a total of 163094 reads, and
was classified as *Phytophthora pseudosyringae* (NCBI taxid 221518). Another
common sequence has been matched to two closely related species *Phytophthora
cambivora* (NCBI taxid 53983) and *Phytophthora x cambivora* (NCBI taxid
2056922).

If you are familiar with the command line search tool ``grep`` and the regular
expression syntax, you should find the format of these intermediate TSV files
lends itself to some simple searches. For example, you could see which samples
had matches to *Phytophthora rubi* using ``grep`` as follows:

.. code:: console

    $ grep "Phytophthora rubi" summary/thapbi-pict.ITS1.all_reads.onebp.tsv
    d8613e80b8803b13f7ea5d097f8fe46f_899  129364  Phytophthora rubi
    $ grep d8613e80b8803b13f7ea5d097f8fe46f intermediate/ITS1/*.fasta
    intermediate/ITS1/DNA10MIX_bycopynumber.fasta:>d8613e80b8803b13f7ea5d097f8fe46f_279
    intermediate/ITS1/DNA10MIX_diluted25x.fasta:>d8613e80b8803b13f7ea5d097f8fe46f_349
    intermediate/ITS1/DNA10MIX_undiluted.fasta:>d8613e80b8803b13f7ea5d097f8fe46f_271

The summary reports would also answer this particular question, but this kind
of search can be useful for exploring specific questions.
