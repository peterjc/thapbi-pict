Classifying sequences
=====================

Running thapbi-pict classify
----------------------------

.. tip:

   If you don't have the FASTQ files, just the FASTA files, start from here.

The second stage of the pipeline is to merge all the sample specific FASTA
files into one non-redundant sequence vs sample TSV file, ready to classify all
the unique sequences in it. These steps can be run separately:

.. code:: console

    $ thapbi_pict sample-tally -h
    ...
    $ thapbi_pict classify -h
    ...

There are a number of options here, but for the purpose of this worked example
we will stick with the defaults and tell it to look for FASTA files in the
``intermediate/`` directory.

.. code:: console

    $ thapbi_pict sample-tally -i intermediate/ITS1/*.fasta -o summary/thapbi-pict.ITS1.tally.tsv
    ...
    $ thapbi_pict classify -i summary/thapbi-pict.ITS1.tally.tsv
    ...

Here we have not set the output folder with ``-o`` or ``--output``, which
means the classify step will default to writing the classifier TSV output file
next to the input tally TSV file. There should now be two new files:

.. code:: console

    $ ls -1 summary/thapbi-pict.ITS1.*.tsv
    summary/thapbi-pict.ITS1.onebp.tsv
    summary/thapbi-pict.ITS1.tally.tsv

If you have the ``biom-format`` Python library installed, adding ``--biom`` to
the command line will result in a ``summary/thapbi-pict.ITS1.onebp.biom`` file
as well, equivalent to the data in ``summary/thapbi-pict.ITS1.tally.tsv`` but
potentially more useful for export to other analysis tools.

Abundance thresholds
--------------------

As you might gather from reading the sample-tally command line help, there
are a few settings for read abundance thresholds.

For the minimum read absolute abundance threshold, we have ``-a`` or
``--abundance`` (default 100 suitable for MiSeq, try 1000 for HiSeq), and
``-n`` or ``--negctrls`` for specifying negative controls (default none).

Related to this, there is a relative of fractional abundance threshold ```-f``
or ``--abundance-fraction``, and ``-y`` or ``--synctrls`` to specify and
synthetic control samples. See :ref:`negative_controls` which discusses this.

If any negative controls are specified, those paired FASTQ files are processed
*first*. If any of these contained ITS1 sequences above the specified minimum
absolute abundance threshold (default 100), that higher number is used as
the minimum abundance threshold for the non-control samples. For example, say
one control had several ITS1 sequences with a maximum abundance of 124, and
another control had a maximum ITS1 abundance of 217, while the remaining
controls had no ITS1 sequence above the default level. In that case, the tool
would take maximum 217 as the abundance threshold for the non-control samples.

.. WARNING::

    Setting the abundance threshold low (say under 50 with MiSeq samples)
    risks background contamination coming through into the results. Do not do
    this without strong justification (e.g. look at suitable controls over
    multiple plates from your own laboratory procedure).

.. WARNING::

    Setting the abundance threshold *very* low (under 10) has the additional
    problem that the number of unique sequences accepted will increase many
    times over. This will *dramatically* slow down the rest of the analysis.
    This is only advised for investigating single samples.

For the woody host data, each plate had a negative control sample which should
contain no ITS1 sequences. We can specify the negative controls with ``-n`` or
``--negctrls`` by entering the four FASTQ filenames in full, but since they
have a common prefix we can use a simple wildcard.

For this sample data, happily neither of the negative controls have any ITS1
present above the default threshold, so this would have no effect.

For the THAPBI Phyto-Threats project we now run each 96-well PCR plate with
multiple negative controls. Rather than a simple blank, these include a known
mixture of synthetic sequences of the same length, same nucleotide
composition, and also same di-nucleotide composition as real *Phytophthora*
ITS1. This means we might have say 90 biological samples which should contain
ITS1 but not the synthetics controls, and 6 negative controls which should
contain synthetic controls but not ITS1.

We therefore run ``thapbi_pict sample-tally`` separately for each plate,
where any ITS1 contamination in the synthetic controls is used to set a plate
specific minimum abundance. This means we cannot run ``thapbi_pict pipeline``
on multiple plates at once (although we could run it on each plate, we
generally want to produce reports over multiple plates).


Intermediate TSV files
----------------------

For each input tally TSV file ``<name>.tally.tsv`` another plain text TSV file
is generated named ``<name>.<method>.tsv`` where the default method is
``onebp`` (which looks for perfect matches or up to one base pair different).
These are both sequence versus sample observation tables of counts, but with
sample metadata in header lines (starting with ``#``) and additional columns
for the amplicon marker sequence, and for the classifier output also the NCBI
taxid(s), and genus-species of any classification(s).

These files are not really intended for human use, but are readable. Here we
skip ten lines of sample metadata at the start, and all the sample-specific
counts in columns 2 to 123, and the sequence in column 124, showing just the
first and final two columns:

.. code:: console

    $ tail -n +10 summary/thapbi-pict.ITS1.onebp.tsv | head | cut -f 1,125,126
    <SEE TABLE BELOW>

Viewing it like this is not ideal, although there are command line tools which
help. You could instead open the file in R, Excel, etc:

============================================ ===================== ================================================================
#Marker/MD5_abundance                        taxid                 genus-species
============================================ ===================== ================================================================
ITS1/2e4f0ed53888ed39a2aee6d6d8e02206_163094 221518                Phytophthora pseudosyringae
ITS1/d9bc3879fdab3b4184c04bfbb5cf6afb_83653  631361                Phytophthora austrocedri
ITS1/32159de6cbb6df37d084e31c37c30e7b_28976  67594                 Phytophthora syringae
ITS1/ed15fefb7a3655147115fc28a8d6d671_28786  78237                 Phytophthora gonapodyides
ITS1/972db44c016a166de86a2bacab3f4226_28515  53983                 Phytophthora x cambivora
ITS1/c1a720b2005f101a9858107545726123_20400  78237                 Phytophthora gonapodyides
ITS1/96e0e2f0475bd1617a4b05e778bb04c9_17392  78237                 Phytophthora gonapodyides
ITS1/f27df8e8755049e831b1ea4521ad6eb3_16369  2496075;2897317;29920 Phytophthora aleatoria;Phytophthora alpina;Phytophthora cactorum
ITS1/21d6308d89d74b8ed493d73a2cb4adb5_16169  53983                 Phytophthora x cambivora
============================================ ===================== ================================================================

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

    $ grep "Phytophthora rubi" summary/thapbi-pict.ITS1.onebp.tsv | cut -f 1,125,126
    ITS1/d8613e80b8803b13f7ea5d097f8fe46f_899  129364  Phytophthora rubi
    $ grep d8613e80b8803b13f7ea5d097f8fe46f intermediate/ITS1/*.fasta
    intermediate/ITS1/DNA10MIX_bycopynumber.fasta:>d8613e80b8803b13f7ea5d097f8fe46f_279
    intermediate/ITS1/DNA10MIX_diluted25x.fasta:>d8613e80b8803b13f7ea5d097f8fe46f_349
    intermediate/ITS1/DNA10MIX_undiluted.fasta:>d8613e80b8803b13f7ea5d097f8fe46f_271

The summary reports would also answer this particular question, but this kind
of search can be useful for exploring specific questions.
