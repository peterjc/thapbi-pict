Worked Example
==============

The *Quick Start* described a simplified use of the THAPBI PICT tool to
assess a single Illumina MiSeq sequencing run using the ``thapbi_pict
pipeline`` command, as a flowchart:

.. image:: images/pipeline.svg
   :alt: Flowchart summarising THAPBI PICT pipeline, from raw paired FASTQ files to reports.

Here we will run over the same process using real data, calling the individual
commands within the default pipeline - and include metadata for the reporting.
We will close with the equivalent all-in-one pipeline command.

.. image:: images/pipeline-meta.svg
   :alt: Flowchart summarising THAPBI PICT pipeline, from raw paired FASTQ files to reports, using metadata.

In these illustrative flow charts of the default pipeline, the input paired
FASTQ files (and metadata) are green, the intermediate per-sample FASTA and
TSV files are yellow, and the output reports are in orange. The individual
steps of the pipeline are dark blue boxes, and the ITS1 database is a pale
blue cylinder.

Sample Data
-----------

This example is based on the following paper from earlier in the THAPBI
Phyto-Threats project, where the original analysis used the precursor pipeline
``metapy``:

* Riddell et al (2019) Metabarcoding reveals a high diversity of woody
  host-associated Phytophthora spp. in soils at public gardens and amenity
  woodlands in Britain. https://doi.org/10.7717/peerj.6931

The raw data is from two Illumina MiSeq runs, a whole 96-well plate from 2016,
and about half the samples from a second 96-well plate sequenced in 2017
(where the rest of the plate was samples from a separate ITS1 study). There
are multiple replicates from each of 14 sample sites, plus controls.

The raw FASTQ files are too large to include with the THAPBI PICT source code,
so to follow the complete example you must download 244 ``*.fastq.gz`` files
separately (122 pairs, a little over 200MB in total).

The first step of a typical THAPBI PICT workflow is to transform the paired
FASTQ files into much smaller FASTA files. We provide those FASTA files
compressed with the THAPBI PICT source code, so they can be used to follow the
rest of a typical analysis. We also provide metadata for the samples for use
in the reports.

Setup
-----

We assume you have your command line terminal open in a new empty folder
dedicated to this analysis. Start by making three sub-folders as follows:

.. code:: console

   $ mkdir raw_data/ intermediate/ summary/

We will need file ``site_metadata.tsv`` (included with the THAPBI PICT source
code as ``tests/woody_hosts/site_metadata.tsv``) which can be downloaded:

.. code:: console

    $ wget https://github.com/peterjc/thapbi-pict/raw/master/tests/woody_hosts/site_metadata.tsv

The FASTQ files are only needed for the very first step of the worked example.
If you have downloaded the 244 paired FASTQ files, put them in the raw data
sub-folder as ``raw_data/*.fastq.gz``.

If you don't have the FASTQ files, you need get the pre-prepared 122 FASTA
files into your intermediate data sub-folder as ``intermediate/*.fasta``.
These are provided as a small compressed file included in the THAPBI PICT
source code ``tests/woody_hosts/woody_hosts_fasta.tar.bz2``, or can easily be
downloaded:

.. code:: console

   $ wget https://github.com/peterjc/thapbi-pict/raw/master/tests/woody_hosts/woody_hosts_fasta.tar.bz2
   $ tar -jxvf woody_hosts_fasta.tar.bz2 -C intermediate/

Note that four of the FASTA files are empty, ``Site_13_sample_7.fasta`` and
``Site_9_sample_4-3.fasta`` (nothing above the minimum threshold), and both
negative controls (good).

Running thapbi-pict prepare-reads
---------------------------------

Calling ``thapbi-pict prepare-reads`` is the first action done by the top
level ``thapbi_pict pipeline`` command.

.. code:: console

    $ thapbi_pict prepare-reads -h
    ...

Assuming you have the FASTQ files in ``raw_data/`` as described above:

.. code:: console

    $ thapbi_pict prepare-reads -i raw_data/ -o intermediate/
    ...

For each input FASTQ file pair ``raw_data/<sample_name>_R1.fastq.gz`` and
``raw_data/<sample_name>_R2.fastq.gz`` you should get a small FASTA file
``intermediate/<sample_name>.fasta``. In this cases there are multiple
replicates from each of 14 sample sites where the file name stem is
``Site_<N>_sample_<X>``, plus the controls.

.. code:: console

    $ ls -1 intermediate/*.fasta | wc -l
    122

You should find 122 small FASTQ files in the ``intermediate/`` folder (or you
can get these from the compressed file as described above). Note this is
robust to being interupted and restarted (e.g. a job might time out on the
cluster).

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
  do not overlap are discarded, for example from unexpectedly long fragements,
  or not enough left after quality trimming).
* Primer trim (reads without both primers are discarded).
* Convert into a non-redundant FASTA file, with the sequence name recording
  the abundance (discarding sequences of low abundance).
* Filter with Hidden Markov Models (HMMs) of ITS1 and our four synthetic
  controls (non-matching sequences are discarded).

For each input ``<sample_name>_R1.fastq.gz`` and ``<sample_name>_R2.fastq.gz``
FASTQ pair we get a single *much* smaller FASTA file ``<sample_name>.fasta``.

.. WARNING::

   The intermediate FASTA files can legitimately be empty when no sequences
   passed the thresholds. This can happen when a PCR failed, and is expected
   to happen on blank negative controls.

The sequence entries in the FASTA file are named ``<checksum>_<abundance>``.
Here ``<checksum>`` is the `MD5 checksum <https://en.wikipedia.org/wiki/MD5>`_
of the sequence, and this is used as a unique shorthand. It is a 32 character
string of the digits ``0`` to ``9`` and lower cases letters ``a`` to ``f``
inclusive, like ``a559aa4d00a28f11b83012e762391259``. These MD5 checksums are
used later in the pipeline, including in reports. The ``<abundance>`` is just
an integer, the number of paired reads which after processing had this unique
sequence.

The description entry in the FASTA file is currently just the name of the HMM
it matched, allowing us to distinguish the biological ITS1 sequences from the
synthetic controls.

Finally, the sequence in the FASTA file is written as a single line in upper
case. With standard FASTA line wrapping at 60 or 80 characters, the ITS1
sequences would need a few lines each. However, they are still short enough
that having them on one line without line breaks is no hardship - and it is
*extremely* helpful for simple tasks like using ``grep`` to look for a
particular sequence fragment at the command line.

For example,

.. code:: console

    $ cat intermediate/Site_1_sample_1.fasta
    >a559aa4d00a28f11b83012e762391259_2303 phytophthora_its1
    CCACACCTAAAAAACTTTCCACGTGAACTGTATCGAACAACTAGTTGGGGGTCTTGTTTGGCGTGCGGCTGCTTCGGTAGCTGCTGCTAGGCGAGCCCTATCACGGCGAGCGTTTGGACTTCGGTCTGAGCTAGTAGCTATTTTTTAAACCCATTCTTTAATACTGATTATACT
    >140ccd03a87b423a1f06521f08131464_724 phytophthora_its1
    CCACACCTAAAAAAACTTTCCACGTGAACCGTATCAACCCCTATAATTTGGGGGCTTGCTCGGCGGCGTGTGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCGGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAGCTAGTAGCTATCAATTTTAAACCCTTTCTTAAATACTGAACATACT
    >868e1ad838c7ec587dfd05b9dd4556ec_339 phytophthora_its1
    CCACACCTAAAAAAAACTTTCCACGTGAACCGTATCAACCCCTATAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCGGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAGCTAGTAGCTATCAATTTTAAACCCTTTCTTAAATACTGAACATACT
    >742f1f7a934f2df075be6f2eea756fc9_210 phytophthora_its1
    CCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAACCGTTAGTTGGGGGCTTCTGTTCGGCTGGCTTCGGCTGGCTGGGCGGCGGCTCTATCATGGCGAGCGCTTGAGCCTTCGGGTCTGAGCTAGTAGCCCACTTTTTAAACCCATTCCTAAATACTGAATATACT
    >7f27d3a8f7150e0ee7ad64073e6da6b5_193 phytophthora_its1
    CCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAACCCTTAGTTGGGGGCTTCTGTTCGGCTGGCTTCGGCTGGCTGGGCGGCGGCTCTATCATGGCGAGCGCTTGAGCCTTCGGGTCTGAGCTAGTAGCCCACTTTTTAAACCCATTCCTAAATACTGAATATACT
    >eaf42569c8b95c8bf4f9bf1b65a96ce4_183 phytophthora_its1
    CCACACCTAAAAAACTTTCCACGTGAACCGTATCAACCCACTTAGTTGGGGGCTAGTCCCGGCGGCTGGCTGTCGATGTCAAAGTTGACGGCTGCTGCTGTGTGTCGGGCCCTATCATGGCGAGCGTTTGGGTCCCTCTCGGGGGAACTGAGCCAGTAGCCCTTATTTTTTAAACCCATTCTTGAATACTGAATATACT
    >ffb8fbb83fa26a101c2fddf2af13cf95_167 phytophthora_its1
    CCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAATCCTTTTATTGGGGGCTTCTGTCTGGTCTGGCTTCGGCTGGTCTGGGTGGCGGCTCTATCATGGTGACCGCTCTGGGCTTCGGCTTGGAGTTAGTAGCCCACTTTTTAAACCCATTCTTAATTACTGAACATACT
    >af3654932ad7a06c5f4af3c738706c76_114 phytophthora_its1
    CCACACCTAAAAAAACTTTCCACGTGAACCGTATCAACCCCTATAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCGGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAGCTAGTAGCTATCAATTTTAAACCCTTTCTTAAATACTGAACATACT

We see this sample had eight unique sequences accepted, all matched the ITS1
HMM (happily none match the synthetic controls). The most common had MD5
checksum ``a559aa4d00a28f11b83012e762391259`` and was seen in 2303 reads.

You could easily find out which other samples had this unique sequence using
the command line search tool ``grep`` as follows:

.. code:: console

    $ grep a559aa4d00a28f11b83012e762391259 intermediate/*.fasta
    ...

You can also answer this example question from the read report produced later.

Abundance thresholds
--------------------

As you might gather from reading the command line help, there are two settings
to do with the minimum read abundance threshold, ``-a`` or ``--abundance``
(default 100), and ``-n`` or ``--negctrls`` for specifying negative controls
(default none).

If any negative controls are specified, those paired FASTQ files are processed
*first*, using the specified minimum abundance (default 100). If any of these
contained ITS1 sequences above the threshold, that higher number is used as
the minimum abundance threshold for the non-control samples. For example, say
one control had several ITS1 sequences with a maximum abundance of 124, and
another control had a maximum ITS1 abundance of 217, while the remaining
controls had no ITS1 sequence above the default level. In that case, the tool
would take maximum 217 as the abundance threshold for the non-control samples.

For example, to lower the threshold from the default to 50, you could use:

.. code:: console

    $ rm -rf intermediate/*.fasta
    $ thapbi_pict prepare-reads -i raw_data/ -o intermediate/ -a 50
    ...

.. WARNING::

   By default ``thapbi_pict prepare-reads`` and ``thapbi_pict pipeline`` will
   reuse existing intermediate FASTA files, so you must explicitly delete any
   old FASTA files before the new abundance threshold will have any effect.

.. WARNING::

    Setting the abundance threhold low (say under 50) risks letting background
    contamination through into the results. Do not do this without strong
    justification (e.g. look at suitable controls over multiple plates from
    your own laboratory procedure).

.. WARNING::

    Setting the abundance threshold *very* low (under 10) has the additional
    problem that the number of unique sequences accepted will increase many
    times over. This will *dramatically* slow down the rest of the analysis.
    This is only advised for investigating single samples.

For the woody host data, each plate had a negative control sample which should
contain no ITS1 sequences. We can specify the negative controls with ``-n`` or
``--negctrls`` by entering the four FASTQ filenames in full, but since they
have a common prefix we can use a simple wild card:

.. code:: console

    $ thapbi_pict prepare-reads -i raw_data/ -o intermediate/ -n raw_data/NEGATIVE*.fastq.gz
    ...

For this sample data, happily neither of the negative controls have any ITS1
present above the default threshold, so this would have no effect.

For the THAPBI project we now run each 96-well PCR plate with multiple
negative controls. Rather than a simple blank, these include a known mixture
of synthetic sequences of the same length, same nucelotide composition, and
also same di-nucleotide composition as real *Phytophthora* ITS1. This means we
might have say 90 biological samples which should contain ITS1 but not the
synthetics controls, and 6 negative controls which should contain synthetic
controls but not ITS1. We then run ``thapbi_pict prepare-reads`` separately
for each plate, where any ITS1 contamination in the synthetic controls is
used to set a plate specific minimum abundance. This means we cannot run
``thapbi_pict pipeline`` on multiple plates at once (although we could run it
on each plate, we generally want to produce reports over multiple plates).

Running thapbi-pict classify
----------------------------

.. tip:

   If you don't have the FASTQ files, just the FASTA files, start from here.

The second stage of the pipeline can be run separately as the ``thapbi_pict
classify`` command:

.. code:: console

    $ thapbi_pict classify -h
    ...

There are a number of options here, but for the purpose of this worked example
we will stick with the defaults and tell it to look for FASTA files in the
``intermediate/`` directory.

.. code:: console

    $ thapbi_pict classify -i intermediate/
    ...

Here we have not set the output folder with ``-o`` or ``--output``, which
means the tool will default to writing the TSV output files next to each
input FASTA file. There should now be 122 TSV files, one for each FASTA:

.. code:: console

    $ ls -1 intermediate/*.tsv | wc -l
    122

Intermediate TSV files
----------------------

For each FASTA file named ``<sample_name>.fasta`` a plain text tab separated
variable (TSV) file is generated named ``<sample_name>.<method>.tsv`` where
the default method is ``onebp`` (this looks for perfect matches or up to one
base pair different). The first line is a header comment line (starting with
``#``) labelling the columns, which are:

* Unique sequence name in ``<checksum>_<abundance>`` format.
* NCBI taxid of any predictions (semi-colon separated, as order as species)
* Genus-species of any predictions (semi-colon separated, alphabetical)
* Text note field (arbitrary debug text from the tool)

These files are not really intended for human use, but are readable:

.. code:: console

    $ cat intermediate/Site_1_sample_1.onebp.tsv
    ...

Viewing it like this is not ideal, although there are command line tools which
help. You could open the file in R, Excel, etc. Slightly abridged, we have:

========================================= ============= ================================================= ====
#sequence-name                            taxid         genus-species:...                                 note
========================================= ============= ================================================= ====
``a559aa4d00a28f11b83012e762391259_2303`` 221518        *Phytophthora pseudosyringae*                     ...
``140ccd03a87b423a1f06521f08131464_724``  78237         *Phytophthora gonapodyides*                       ...
``868e1ad838c7ec587dfd05b9dd4556ec_339``  78237         *Phytophthora gonapodyides*                       ...
``742f1f7a934f2df075be6f2eea756fc9_210``  164328        *Phytophthora ramorum*                            ...
``7f27d3a8f7150e0ee7ad64073e6da6b5_193``  164328        *Phytophthora ramorum*                            ...
``eaf42569c8b95c8bf4f9bf1b65a96ce4_183``  53983;2056922 *Phytophthora cambivora;Phytophthora x cambivora* ...
``ffb8fbb83fa26a101c2fddf2af13cf95_167``  631361        *Phytophthora austrocedri*                        ...
``af3654932ad7a06c5f4af3c738706c76_114``  78237         *Phytophthora gonapodyides*                       ...
========================================= ============= ================================================= ====

This says most of the unique sequences here have been assigned a single unique
*Phytophthora* species, except for ``eaf42569c8b95c8bf4f9bf1b65a96ce4`` (found
in 183 reads for this sample) which has matched *Phytophthora cambivora* (NCBI
taxid 53983) and close relative *Phytophthora x cambivora* (NCBI taxid
2056922).

If you are familiar with the command line search tool ``grep`` and the regular
expression syntax, you should find the format of these intermediate TSV files
lends itself to some simple searches. For example, you could see which samples
had matches to *Phytophthora rubi* using ``grep`` twice as follows (exclude
header lines, then find species):

.. code:: console

    $ grep -v "^#" intermediate/*.tsv | grep "Phytophthora rubi"
    intermediate/DNA10MIX_bycopynumber.onebp.tsv:2ba87367bdbb87cc37521bed773ffa37_285  129364  Phytophthora rubi  Unique taxonomy match
    intermediate/DNA10MIX_diluted25x.onebp.tsv:2ba87367bdbb87cc37521bed773ffa37_363    129364  Phytophthora rubi  Unique taxonomy match
    intermediate/DNA10MIX_undiluted.onebp.tsv:2ba87367bdbb87cc37521bed773ffa37_274     129364  Phytophthora rubi  Unique taxonomy match

The summary reports would also answer this paricular question, but this kind
of search can be useful for exploring specific questions.

Metadata
--------

The *Quick Start* introduced the typical pipeline taking paired FASTQ files
though to reports, and mentioned the idea of enhancing the reports with
sample metadata.

.. image:: images/pipeline-meta.svg
   :alt: Flowchart summarising THAPBI PICT pipeline, from raw paired FASTQ files to reports, using metadata.

In the following we will show the reports with and without metadata.
File ``site_metadata.tsv`` is a table of metadata (based on table S1 in the
paper), in plain text tab separated variable format (TSV). It has one row for
each of the 14 samples plus controls, with a new column cross referencing the
122 sequenced FASTQ filename stems.

This metadata file is perhaps unusual in that it has a header of comment lines
(starting ``#``) which some tools like R and THAPBI PICT need to be told to
ignore explicitly. Quoting from that header::

    # Lines 1-19, human readable header text
    # Lines 20, colum headers
    # Lines 21 onwards, data for 14 field sites and 3 controls

When calling THAPBI PICT we need to include ``-f 20`` or ``--metafields 20``
indicating the column headers are on line 20 (and lines 1 to 19 should be
ignored). Typically the column header will be on line one, so this is not
needed.

As to the columns, again quoting from that header::

    # Column 1 (A), Anonymised site number with leading zero added (01 to 14), or control name
    # Column 2 (B), Approximate altitude at centre
    # Column 3 (C), underlying soil type
    # Columns 4 to 15 (D to O): Tree/shrub broad taxonomic grouping and health status (H, healthy; D, symptoms/stump/dead)
    # Column 16 (P): Semi-colon separated list of Illumina MiSeq sample names

We will be using ``-c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15`` (or ``--metacols``)
meaning show columns 1 to 15 inclusive in the reports (in that order).

Finally, we will use ``-x 16`` or ``--metaindex 16`` to indicate column 16
contains cross references to the sequenced sample filename stems (semi-colon
separated).

This cross referencing idea is key to getting the best results from attaching
metadata to your sequenced samples. Here is an abridged representation of the
table, showing column one (site or control name), column two (altitute), and
finally column 16 which has the filename stems of the sequence data belonging
to this row of the table (semi-colon separated list).

======== ======== === ========================================================
#Site    Altitude ... MiSeq Sample(s)
======== ======== === ========================================================
01             30 ... Site_1_sample_1;Site_1_sample_2;Site_1_sample_3;Site_1_sample_4;Site_1_sample_5;Site_1_sample_6;Site_1_sample_7;Site_1_sample_8;Site_1_sample_9-2;Site_1_sample_10
02             55 ... Site_2_sample_1;Site_2_sample_2;Site_2_sample_3;Site_2_sample_4;Site_2_sample_5;Site_2_sample_6;Site_2_sample_7;Site_2_sample_8;Site_2_sample_9;Site_2_sample_10
03             45 ... Site_3_sample_1;Site_3_sample_2;Site_3_sample_4;Site_3_sample_7;Site_3_sample_8;Site_3_sample_9
04             20 ... Site_4_sample_1;Site_4_sample_2;Site_4_sample_3;Site_4_sample_3-2;Site_4_sample_4;Site_4_sample_5;Site_4_sample_6;Site_4_sample_8;Site_4_sample_9;Site_4_sample_10
05            100 ... Site_5_sample_1;Site_5_sample_2;Site_5_sample_4;Site_5_sample_5;Site_5_sample_6;Site_5_sample_8;Site_5_sample_9
06              5 ... Site_6_sample_1;Site_6_sample_2-2;Site_6_sample_3-1;Site_6_sample_4;Site_6_sample_5-3;Site_6_sample_6;Site_6_sample_7-1;Site_6_sample_8-2;Site_6_sample_9;Site_6_sample_10
07            105 ... Site_7_sample_1;Site_7_sample_2;Site_7_sample_3;Site_7_sample_5;Site_7_sample_6;Site_7_sample_7;Site_7_sample_8;Site_7_sample_9;Site_7_sample_10
08             45 ... Site_8_sample_1;Site_8_sample_2;Site_8_sample_3;Site_8_sample_4;Site_8_sample_5-2;Site_8_sample_6;Site_8_sample_7;Site_8_sample_7-2;Site_8_sample_8;Site_8_sample_9
09             15 ... Site_9_sample_1;Site_9_sample_4-3;Site_9_sample_6;Site_9_sample_7;Site_9_sample_8;Site_9_sample_9;Site_9_sample_10
10             30 ... Site_10_sample_7;Site_10_sample_8
11             80 ... Site_11_sample_1;Site_11_sample_2;Site_11_sample_3;Site_11_sample_4;Site_11_sample_5;Site_11_sample_6;Site_11_sample_7;Site_11_sample_8;Site_11_sample_9;Site_11_sample_10
12             30 ... Site_12_sample_1;Site_12_sample_2;Site_12_sample_3-3;Site_12_sample_4;Site_12_sample_5-3;Site_12_sample_6;Site_12_sample_8;Site_12_sample_9;Site_12_sample_10
13            300 ... Site_13_sample_1;Site_13_sample_2;Site_13_sample_4;Site_13_sample_5;Site_13_sample_6;Site_13_sample_7;Site_13_sample_8;Site_13_sample_9;Site_13_sample_10
14             30 ... Site_14_sample_1-2;Site_14_sample_2;Site_14_sample_3;Site_14_sample_4;Site_14_sample_5;Site_14_sample_6;Site_14_sample_10
DNA10MIX          ... DNA10MIX_diluted25x;DNA10MIX_undiluted;DNA10MIX_bycopynumber
DNA16MIX          ... DNA16MIX
NEGATIVE          ... NEGATIVE_firstplate;NEGATIVE_secondplate
======== ======== === ========================================================

Also note that in column one we have listed the numerical site names with
leading zeros giving ``01`` to ``14`` to ensure they sort as expected.

Running thapbi_pict sample-report
---------------------------------

The first output reports from the pipeline can be generated separately by the
``thapbi_pict sample-summary`` command:

.. code:: console

    $ thapbi_pict sample-summary -h
    ...

To mimic what the pipeline command would do, run the following, where we have
to give names to the computer readable TSV (``-o`` or ``--output``) and human
readable TXT outputs (``-r`` or ``--report``):

.. code:: console

    $ thapbi_pict sample-summary -i intermediate/ \
      -o summary/thapbi-pict.samples.tsv -r summary/thapbi-pict.samples.txt
    ...

Note the trailing slash ``\`` at the end of the first line indicates the
command continues on the next line. You can actually type this at the standard
Linux command prompt (or include it in a copy and paste), or just enter this
as one very long command.

We will look at the output in a moment, along side the equivalent reports
generated with metadata (see discussion above about column and row numbers):

.. code:: console

    $ time thapbi_pict sample-summary -i intermediate/ -o \
      summary/with-metadata.samples.tsv -r summary/with-metadata.samples.txt \
      -t site_metadata.tsv -c 1,2,3,4,5,6,7,8,9,10,11,12,13,15 -x 16 -f 20
    ...

Without metadata the results are sorted by sample name (and thus the first entry
here should be ``DNA10MIX_bycopynumber``) but with a metadata table given, that
order is used (with any metadata-less entries at the end).

Because the computer readable TSV file does not include the metadata directly,
the two versions differ only in the line order:

.. code:: console

    $ diff <(sort summary/thapbi-pict.samples.tsv) <(sort summary/with-metadata.samples.tsv)
    $ diff summary/thapbi-pict.samples.tsv summary/with-metadata.samples.tsv
    ...

The human readable text file is more interesting, and is discussed next.

Sample Report
-------------

Here we will discuss the high level human readable summary report from
``thapbi_pict sample-report``, produced as plain text.

Without metadata, the samples are sorted by filename alphabetically. In this
example that means we get DNA controls, Negative controls, Site 1, 10, 11, 12,
13, 14, 2, ..., 9. This is unfortunate, can can be resolved by providing some
minimal metadata - or better yet, using leading zeros in all sample names.

Pulling out the ``Site_1_sample_1`` example, we see:

    Site_1_sample_1

     - Phytophthora austrocedri
     - Phytophthora cambivora (uncertain/ambiguous)
     - Phytophthora gonapodyides
     - Phytophthora pseudosyringae
     - Phytophthora ramorum
     - Phytophthora x cambivora (uncertain/ambiguous)

    Site_1_sample_10

    ...

Note we get ``Site_1_sample_1`` then ``Site_1_sample_10`` and then
``Site_1_sample_2`` etc. The sort order problem strikes again!

As to the meaning of this list, those are the species identified - in some
cases with a cavaet. The file starts with a tiny explanation:

    NOTE: Species listed with (uncertain/ambiguous) in brackets are where
    sequences matched multiple species equally well. For example,
    *Phytophthora andina*, *P. infestans*, and *P. ipomoeae*, share an
    identical marker.

In this case, as you may recall from when we looked at the classifier output
file ``intermediate/Site_1_sample_1.onebp.tsv``, one of the sequences matched
both *Phytophthora cambivora* and *Phytophthora x cambivora*.

In comparison, with metadata, all the samples matched to a metadata row get
grouped with a shared metadata header::

    Site: 01
    Altitude (m): 30
    Underlying soil type: Brown earth, moderately well drained
    H/Cupressaceae: 0
    D/Cupressaceae: 1
    H/Other conifers: 0
    D/Other conifers: 1
    H/Ericaceae: 0
    D/Ericaceae: 4
    H/Fagaceae or Nothofagaceae: 2
    D/Fagaceae or Nothofagaceae: 1
    H/Other angiosperms: 0
    D/Other angiosperms: 1
    D/Other: 0

    Sequencing sample: Site_1_sample_1

     - Phytophthora austrocedri
     - Phytophthora cambivora (uncertain/ambiguous)
     - Phytophthora gonapodyides
     - Phytophthora pseudosyringae
     - Phytophthora ramorum
     - Phytophthora x cambivora (uncertain/ambiguous)

    Sequencing sample: Site_1_sample_2

     - Unknown
     - Phytophthora austrocedri
     - Phytophthora cambivora (uncertain/ambiguous)
     - Phytophthora gonapodyides
     - Phytophthora pseudosyringae
     - Phytophthora ramorum
     - Phytophthora x cambivora (uncertain/ambiguous)

    ...

If any of the requested metadata is missing, then it does not get shown. For
example, this applies to the DNA mixes and the negative controls.

Note that for ``Site_1_sample_2``, at least one unique sequence was not given
a species or even genus level classification. This likely reflects a gap in
the default database, and/or the default method being too strict.

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
