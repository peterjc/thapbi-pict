Worked Example
==============

The *Quick Start* described a simplified use of the THAPBI PICT tool to
assess a single Illumina MiSeq sequencing run using the ``thapbi_pict
pipeline`` command. Here we will run over the same process using real data,
the individual commands within the default pipeline - and include metadata
for the reporting. We will close with the equivalent all in one pipeline
command.

Sample Data
-----------

This example is based on the following paper from earlier in the THAPBI
project, where the original analysis used the precursor pipeline ``metapy``:

* Riddell et al (2019) Metabarcoding reveals a high diversity of woody
  host-associated Phytophthora spp. in soils at public gardens and amenity
  woodlands in Britain. https://doi.org/10.7717/peerj.6931

The raw data is from two Illumina MiSeq runs (a whole 96-well plate from 2016,
and about half the samples a second 96-well plate sequenced in 2017, the rest
of the plate being samples from a separate ITS1 study).

The raw FASTQ files are too large to include with the THAPBI PICT source code,
so to follow the complete example you must download 244 ``*.fastq.gz`` files
separately (122 pairs, a little over 200MB in total).

There are multiple replicates from each of 14 sample sites, with FASTQ files
``Site_<N>_sample_<X>_R1.fastq.gz`` and ``Site_<N>_sample_<X>_R2.fastq.gz``
(plus controls with a different pattern), which as the first step of the
typical THAPBO PICT workflow (``thapbi_pict prepare-reads``) are transformed
in FASTA files named ``Site_<N>_sample_<X>.fasta`` (etc). We provide these
FASTA files as a compressed file with the THAPBI PICT source code, so after
decompression they can be used to follow the rest of a typical analysis.

We also provide metadata for the samples for use in the reports.

Setup
-----

We assume you have your command line terminal open in a new empty folder
dedicated to this analysis. Start by making three sub-folders as follows:

.. code:: console

   $ mkdir raw_data/ intermediate/ summary/

We will need the ``site_metadata.tsv`` (included with the THAPBI PICT source
code under ``tests/woody_hosts/``, or easily download). This is a table of
metadata (based on table S1 in the paper), with one row for each of the 14
samples plus controls, with a cross reference to the 122 sequenced FASTQ
filename stems. Download or copy this to your project folder:

.. code:: console

    $ wget https://github.com/peterjc/thapbi-pict/raw/master/tests/woody_hosts/site_metadata.tsv

If you have downloaded the 244 paired FASTQ files, put them in the raw data
sub-folder as ``raw_data/*.fastq.gz``, but this is only needed for the very
first step which can be skipped if you instead decompress the provided 122
FASTA file (included in the THAPBI PICT source code, but can easily be
re-downloaded):

.. code:: console

   $ cd intermediate/
   $ wget https://github.com/peterjc/thapbi-pict/raw/master/tests/woody_hosts/woody_hosts_fasta.tar.bz2
   $ tar -jxvf /path/to/downloads/woody_hosts_fasta.tar.bz2
   $ cd ..

Note that four of the FASTA files are empty, ``Site_13_sample_7.fasta`` and
``Site_9_sample_4-3.fasta`` (nothing above the minimum threshold), and both
negative controls (good).

thapbi-pict prepare-reads
=========================

Calling ``thapbi-pict prepare-reads`` is the first action done by the top
level ``thapbi_pict pipeline`` command, as illustrated here:

.. image:: images/pipeline.svg
   :alt: Flowchart summarising THAPBI PICT pipeline, from raw paired FASTQ files to reports.

In this illustrative flow chart of the default pipeline, the input paired
FASTQ files are green, the intermediate per-sample FASTA and TSV files are
yellow, and the output reports are in orange. The individual steps of the
pipeline are dark blue boxes, and the ITS1 database is a pale blue cylinder.

.. code:: console

    $ thapbi_pict prepare-reads -h
    ...

Assuming you have the FASTQ files in ``raw_data/`` as described above:

.. code:: console

    $ thapbi_pict prepare-reads -i raw_data/ -o intermediate/
    ...

You should then find 122 small FASTQ files in the ``intermediate/`` folder
(or you can get these from the compressed file as described above). Note
this is robust to being interupted and restarted (e.g. a job might time
out on the cluster).

.. WARNING::
 
    So far this example omits a key consideration - telling the tool which
    samples are negative controls, and/or manually setting the minimum read
    abundance. See below.

Intermediate FASTA files
------------------------

What the prepare command does can be briefly summarised as follows:

* Quality trim the FASTQ reads (pairs where either read becomes too short are
  discarded).
* Merge the overlapping paired FASTQ reads into single sequences (pairs which
  do not overlap discarded, for example from unexpectedly long fragements, or
  not enough left after quality trimming).
* Primer trim (reads without both primers are discarded).
* Convert into a non-redundant FASTA file, with the sequence name recording
  the abundance, discarding sequences of low abundance.
* Filter with Hidden Markov Models (HMMs) of ITS1 and our four synthetic
  controls (non-matching sequences are discarded).

For each input ``<sample_name>_R1.fastq.gz`` and ``<sample_name>_R2.fastq.gz``
FASTQ pair we get a single much smaller FASTA file ``<sample_name>.fasta``,
the contents of which reflects those last two stages.

.. WARNING::

   The intermediate FASTA files can legitimately be empty when no sequences
   passed the thresholds. This can happen when a PCR failed, and is expected
   to happen on blank negative controls.

The sequence entries in the FASTA file are named ``<checksum>_<abundance>``.
Here ``<checksum>`` is the `MD5 checksum <https://en.wikipedia.org/wiki/MD5>`_
of the sequence, and this is used as a unique shorthand. It is a 32 character
string of the digits ``0`` to ``9`` and lower cases letters ``a`` to ``f``
inclusive. These MD5 checksums are used later in the pipeline, including in
reports. The ``<abundance>`` is just an integer, the number of paired reads
which after processing had this unique sequence.

The description entry in the FASTA file is the name of the HMM it matched,
allowing us to distinguish the biological ITS1 sequences from the synthetic
controls.

Finally, the sequence in the FASTA file is written as a single line in upper
case. With standard line wrapping at 60 or 80 characters, the ITS1 sequences
would need a few lines each. However, they are still short enough that having
them on line line without line breaks is no hardship - and it is extremely
helpful for simple tasks like using ``grep`` to look for a particualr sequence
at the command line.

Abundance thresholds
--------------------

As you might gather from reading the command line help, there are two settings
to do with the minimum read abundance threshold, ``-a`` or ``--abundance``
(default 100), and ``-n`` or ``--negctrls`` for specifying negative controls.

.. WARNING::

   By default ``thapbi_pict prepare-reads`` and	``thapbi_pict pipeline`` will
   reuse existing intermediate FASTA files, so you must	explicitly delete any
   old FASTA files before the new abundance threshold will have any effect.

For example, to	lower the threshold from the default to	50, you	could use:

.. code:: console

    $ rm -rf intermediate/*.fasta
    $ thapbi_pict prepare-reads -i raw_data/ -o intermediate/ -a 50
    ...

.. WARNING::

    Setting the abundance threhold low (say under 50) risks letting background
    contamination through into the results. Do not do this without strong
    justification (e.g. look at suitable controls over multiple plates from
    your own laboratory procedure).

.. WARNING::

    Setting the abundance threshold very low (under 10) has the additional
    problem that the number of unique sequences accepted will increase many
    times over. This will *dramatically* slow down the rest of the analysis.
    This is only advised for investigating single samples.

For the woody host data, each plate had a negative control sample which should
contain no ITS1 sequences (and at the default threshold happily none are
found). We can specify the negative controls with ``-n`` or ``--negctrls`` by
entering their filenames in full, or with an appropriate wild card:

.. code:: console

    $ thapbi_pict prepare-reads -i raw_data/ -o intermediate/ -n raw_data/NEGATIVE*
    ...

In this case happily neither of the negative controls have any ITS1 present
above the default threshold, so this would have no effect.

For the THAPBI project we now run each 96-well PCR plate with multiple
negative controls. Rather than a simple blank, these include a known mixture
of synthetic sequences of the same length, nucelotide composition, and also
di-nucleotide composition as real *Phytophthora* ITS1. This means we might
have say 90 biological samples which should contain ITS1 but not the
synthetics controls, and 9 negative controls which should contain synthetic
controls but not ITS1. We then run ``thapbi_pict prepare-reads`` for each
plate, where any ITS1 contamination in the synthetic controls is used to set
a plate specific minimum abundance.

thapbi-pict classify
--------------------

.. tip:

   If you don't have the FASTQ files, just the FASTA files, start from here.

The second stage of the pipeline can be run separately as the
``thapbi_pict classify`` command.

Intermediate TSV files
----------------------

For each FASTA file a tab separated variable (TSV) file is generated where
the first column is the sequence name in ``<checksum>_<abundance>`` format.

Metadata
--------

The *Quick Start* introduced the typical pipeline taking paired FASTQ files
though to reports, and mentioned the idea of enhancing the reports with
sample metadata.

.. image:: images/pipeline-meta.svg
   :alt: Flowchart summarising THAPBI PICT pipeline, from raw paired FASTQ files to reports, using metadata.

In the following we will show the reports with and without metadata.

Sample Reports
--------------

Two of the output reports from the pipeline can also be generated by the
``thapbi_pict sample-report`` sub-command:

* Human readable file ``thapbi-pict.samples.txt`` (plain text).
* Computer readable file ``thapbi-pict.samples.tsv`` (tab separated
  variables, TSV) which can be openend in R, Excel, or similar.

These aim to give a summary of the species identified within each sample. The
human readable report deliberately does not include read counts as the method
is only semi-quantative - as long as it passed the minimum read abundance,
any unique sequence is included.

The computer readable file is intended to facilitate downstream analysis.

Read Reports
------------

The next two output reports from the pipeline can also be generated by the
``thapbi_pict read-summary`` sub-command:

* Plain table ``thapbi-pict.reads.tsv`` (tab separated variables, TSV) which
  can be opened in R, Excel, or similar.
* Visually formatted table ``thapbi-pict.reads.xlsx`` (Microsoft Excel
  format), with the same content but with colors etc applied.

This read report has a row for each unique sequences. The first columns are
the unique sequence MD5 checksum, any species prediction, the sequence itself,
the number of samples it detected in above the threshold, and the total number
of times this was seen (in samples where it was above the threshold). Then
the main columns (one per sample) list the abundance of each unique sequence
in that sample (if above the threshold).

In the Excel version, conditional formatting is used to highlight the non-zero
counts with a red background.

Edit Graph
----------

The final output report from the pipeline can also be generated by the
``thapbi_pict edit-graph`` sub-command:

* Edit-distance graph ``thapbi-pict.edit-graph.xgmml`` (XGMML, eXtensible
  Graph Markup and Modeling Language) which we recommend opening in `Cytoscape
  <https://cytoscape.org/>`_.

Note that ``thapbi_pict edit-graph`` supports other node-and-edge graph file
formats, and can produce a static PDF image as well using `GraphViz
<http://graphviz.org/>`_ and other dependencies.
