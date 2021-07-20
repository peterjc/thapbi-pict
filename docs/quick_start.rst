.. _quick_start:

Quick Start
===========

Here we describe a simplified use of the THAPBI PICT tool to assess a single
Illumina MiSeq sequencing run. The input data is a set of paired FASTQ files
(one for each sample), perhaps barcoded samples from a 96-well plate.

.. image:: images/pipeline.svg
   :alt: Flowchart summarising THAPBI PICT pipeline, from raw paired FASTQ files to reports.

In this illustrative flow chart of the default pipeline, the input paired
FASTQ files are green, the intermediate per-sample FASTA and TSV files are
yellow, and the output reports are in orange. The individual steps of the
pipeline are dark blue boxes, and the ITS1 database is a pale blue cylinder.

We will now describe how to run the ``thapbi_pict pipeline`` command, which
will process the samples, make classifications, and summary reports.

.. code:: console

    $ thapbi_pict pipeline -h
    ...

Setup
-----

We assume you have a new folder dedicated to this analysis, with a sub folder
``raw_data/`` which contains the demultiplexed paired FASTQ files which are
named like ``<sample_name>_R1.fastq.gz`` and ``<sample_name>_R2.fastq.gz``
as provided by your sequencing centre. The tool understands a few widely used
naming patterns. We recommend that you *do* *not* decompress the FASTQ files
(as ``<sample_name>_R1.fastq`` and ``<sample_name>_R2.fastq``), leaving them
gzip compressed is preferable for disk space.

.. code:: console

    $ cd /path/to/my/project/
    $ ls raw_data/*.fastq.gz
    ...

We will make two additional sub-folders, ``intermediate/`` (for the per-sample
prepared FASTA files and classifier prediction TSV files), and ``summary/``
for the folder level reports.

.. code:: console

    $ mkdir intermediate/ summary/

Running
-------

With that done, we run the ``thapbi_pict pipeline`` command, which for a
single 96 sample Illumina MiSeq run should take a minute or so:

.. code:: console

    $ thapbi_pict pipeline -i raw_data/ -s intermediate/ -o summary/thapbi-pict
    ...
    Wrote summary/thapbi-pict.ITS1.samples.onebp.*
    Wrote summary/thapbi-pict.ITS1.reads.onebp.*
    All done!

This is robust to being interrupted and restarted (as long as you are not
changing settings), and will reuse intermediate files:

.. code:: console

    $ thapbi_pict pipeline  -i raw_data/ -s intermediate/ -o summary/thapbi-pict
    ...
    Skipped 120 previously prepared ITS1 samples
    ...
    Wrote summary/thapbi-pict.ITS1.samples.onebp.*
    Wrote summary/thapbi-pict.ITS1.reads.onebp.*
    All done!

All being well, this will produce a set of report files, with names starting
with the prefix ``summary/thapbi-pict.*`` given as follows:

.. code:: console

    $ ls -1 summary/thapbi-pict.*
    summary/thapbi-pict.ITS1.all_reads.fasta
    summary/thapbi-pict.ITS1.all_reads.onebp.tsv
    summary/thapbi-pict.ITS1.reads.onebp.tsv
    summary/thapbi-pict.ITS1.reads.onebp.xlsx
    summary/thapbi-pict.ITS1.samples.onebp.tsv
    summary/thapbi-pict.ITS1.samples.onebp.txt
    summary/thapbi-pict.ITS1.samples.onebp.xlsx

.. WARNING::

    This minimal example omits a key consideration - telling the tool which
    samples are negative controls, and/or manually setting the minimum read
    abundance.

Intermediate FASTA files
------------------------

The first stage of the pipeline can be run separately as the
``thapbi_pict prepare`` command. Here each pair of FASTQ files named something
like ``<sample_name>_R1.fastq.gz`` and ``<sample_name>_R2.fastq.gz`` is
processed to give a much smaller FASTA format file ``<sample_name>.fasta``
containing all the unique sequences from that sample which resemble an ITS1
sequence (or a synthetic control).

In these FASTA files, each sequence is named as ``<checksum>_<abundance>``
where the `MD5 checksum <https://en.wikipedia.org/wiki/MD5>`_ of the
sequence and is used as a unique shorthand - a 32 character string of the
digits ``0`` to ``9`` and lower cases letters ``a`` to ``f`` inclusive.
These MD5 checksums are used later in the pipeline, including in some reports.

Unusually the intermediate FASTA files start with a header made of multiple
lines starting with ``#``, which record information about the sample for use
in reporting. This includes how many raw reads the FASTQ files had, how many
were left after quality trimming, pair merging, primer trimming and finally
the abundance threshold. Many tools will accept these files as FASTA without
complaint, but some tools require the header be removed.

The second stage of the pipeline can be run separately as the ``thapbi_pict
fasta-nr`` command. This produces a pooled non-redundant FASTA file with all
the observed marker sequences in it (and the total read abundance). This is
file ``summary/thapbi-pict.ITS1.all_reads.fasta`` in the above example.

Intermediate TSV files
----------------------

The third stage of the pipeline can be run separately as the ``thapbi_pict
classify`` command. Here each species predictions are made for each sequence
in the prepared non-redundant FASTA file, generating a tab separated variable
(TSV) file where the first column is the sequence name in
``<checksum>_<abundance>`` format. This is file
``summary/thapbi-pict.ITS1.all_reads.onebp.tsv`` in the above example.

Sample Reports
--------------

The first set of reports from the pipeline or ``thapbi_pict summary`` command
are the sample reports:

* Human readable file ``XXX.samples.onebp.txt`` (plain text).
* Plain table ``XXX.samples.onebp.tsv`` (tab separated variables, TSV) which
  can be opened in R, Excel, or similar.
* Visually formatted table ``thapbi-pict.samples.onebp.xlsx`` (Microsoft Excel
  format), with the same content but with colors etc applied.

These aim to give a summary of the species identified within each sample. The
human readable text report deliberately does not include read counts as the
method is only semi-quantitative - as long as it passed the minimum read
abundance, any unique sequence is included.

The tables have one row for each sample. The main columns give total read
counts, those not matched to anything ("Unknown"), reads matched at species
level (with ambiguous combinations listed explicitly), and reads matched only
to genus level.

In the Excel version, conditional formatting is used to highlight the non-zero
counts with a red background.

Read Reports
------------

The other set of reports from the pipeline or ``thapbi_pict summary`` command
are more detailed being at the level of the unique sequences or reads:

* Plain table ``XXX.reads.onebp.tsv`` (tab separated variables, TSV)
  which can be opened in R, Excel, or similar.
* Visually formatted table ``XXX.reads.onebp.xlsx`` (Microsoft Excel
  format), with the same content but with colors etc applied.

This read report has a row for each unique sequence. The first columns are the
unique sequence MD5 checksum, any species prediction, the sequence itself, the
number of samples it was detected in above the threshold, and the total number
of times this was seen (in samples where it was above the threshold). Then
the main columns (one per sample) list the abundance of each unique sequence
in that sample (if above the threshold).

In the Excel version, conditional formatting is used to highlight the non-zero
counts with a red background.

Edit Graph
----------

While not run by the pipeline, there is a separate ``thapbi_pict edit-graph``
command, where the default output is:

* Edit-distance graph ``XXX.edit-graph.xgmml`` (XGMML, eXtensible
  Graph Markup and Modeling Language) which we recommend opening in `Cytoscape
  <https://cytoscape.org/>`_.

Note that ``thapbi_pict edit-graph`` supports other node-and-edge graph file
formats, and can produce a static PDF image as well using `GraphViz
<http://graphviz.org/>`_ and other dependencies, or a distance matrix.

Next Steps
----------

This minimal example omits a key consideration which is telling the tool which
of the samples are your negative controls and/or manually setting the minimum
read abundance.

Also, interpreting the main reports is much easier if you can provide suitably
formatted :ref:`metadata <metadata>`. Happily, you can quickly re-run the
pipeline and it will re-use any already generated intermediate files.

.. image:: images/pipeline-meta.svg
   :alt: Flowchart summarising THAPBI PICT pipeline, from raw paired FASTQ files to reports, using metadata.

The :ref:`first worked example <woody_hosts>` covers these issues, with
excerpts of the expected output.
