.. _summary_reports:

Summary reports
===============

Running thapbi_pict summary
---------------------------

The reports from the pipeline can be generated separately by the ``thapbi_pict
summary`` command:

.. code:: console

    $ thapbi_pict summary -h
    ...

To mimic what the pipeline command would do, run the following:

.. code:: console

    $ thapbi_pict summary -i intermediate/ \
      summary/thapbi-pict.ITS1.all_reads.onebp.tsv \
      -o summary/thapbi-pict.ITS1
    ...

Note the trailing slash ``\`` at the end of the first line indicates the
command continues on the next line. You can actually type this at the standard
Linux command prompt (or include it in a copy and paste), or just enter this
as one very long command.

We will look at the output in a moment, along side the equivalent reports
generated with :ref:`metadata <metadata>` (see linked discussion about column
numbers):

.. code:: console

    $ thapbi_pict summary -i intermediate/ \
      summary/thapbi-pict.ITS1.all_reads.onebp.tsv \
      -o summary/with-metadata.ITS1 \
      -t metadata.tsv -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -x 16
    ...

Both the read report and sample report are tables, produced as both
computer-friendly plain text tab-separated variable (TSV), and human-friendly
Excel (with colors and conditional formatting).

Read Report
-----------

The heart of the read report is a large table, of unique sequences (ASVs rows)
versus sequenced samples (columns), with read abundance counts. There are
additional columns with sequence information, and when :ref:`metadata` is
present, extra rows at the start with sample information.

This read report has a row for each unique sequence. The first columns are the
marker name (here always "ITS1"), the unique sequence MD5 checksum, any
species prediction, the sequence itself, the number of samples it was detected
in above the threshold, the maximum number of reads with this sequence in any
one sample, and the total number of reads (from samples where it was above the
threshold). Then the main columns (one per sample) list the abundance of each
unique sequence in that sample (if above the threshold).

In the Excel version, conditional formatting is used to highlight the non-zero
counts with a red background. Furthermore, with metadata it will attempt to
assign repeated bands of background color to groups (pink, orange, yellow,
green, blue). In this example, each sample site gets a new color:

.. image:: https://user-images.githubusercontent.com/63959/60735578-ebdcf200-9f4b-11e9-8856-1ab66bd1245b.png
   :alt: Screenshot of Excel showing ``summary/with-metadata.samples.onebp.xlsx`` file.

Typical sample naming schemes will result in replicates as neighbouring
columns - meaning you should see very similar patterns of red (non-zero).
Certainly in this dataset scanning horizontally we do see some sequences
clearly show presence/absence patterns consistent with the samples.

The default row sorting will result in a dominant sequence being followed by
any close variants assigned to the same species. Many of these rows will
represent PCR artefacts found in just one or two samples. This contributes
to the "halo" effect seen in the :ref:`edit_graph` representation, discussed
next.

Sample Report
-------------

The heart of the sample report is a table of samples (rows) versus species
predictions (columns), with read abundance counts. There are additional
columns with sample read counts, and when :ref:`metadata` is present, extra
columns at the start with sample information.

Here is a screenshot of the ``summary/with-metadata.ITS1.samples.onebp.xlsx``
file opened in Excel:

.. image:: https://user-images.githubusercontent.com/63959/76231207-cf046700-621c-11ea-9f3a-cdb0cf539483.png
   :alt: Excel screenshot showing with-metadata.ITS1.samples.onebp.xlsx

The metadata is in the first columns, then the sequence filename stem, a text
summary of the species predictions, some inferred sequence count data, and the
one column for each unique species or ambiguous species combinations.

Using the metadata each site has one or more rows in the same background
color (pink, orange, yellow, green, blue, repeated), with one row for each
time it was sequenced (the per-site sampling).

The values are total read counts for that row/column, with conditional
formatting applied so non-zero entries have a bright red background.

For example, the final rows are the two DNA mixture controls (blue and pink)
and the negative controls (orange). These have almost no metadata, and the
negative controls read counts are all zero.

The plain text table ``with-metadata.ITS1.samples.onebp.xlsx`` is the same,
but without the colors and formatting. The files generated without metadata
(``thapbi-pict.ITS1.samples.onebp.xlsx`` etc) lack the extra columns and the
background color bands.

The files without metadata start with the FASTQ filename stem as the inferred
sample name in column 1:

.. code:: console

    $ cut -f 1 summary/thapbi-pict.ITS1.samples.onebp.tsv | head
    #Sequencing sample
    DNA10MIX_bycopynumber
    DNA10MIX_diluted25x
    DNA10MIX_undiluted
    DNA15MIX
    NEGATIVE_firstplate
    NEGATIVE_secondplate
    Site_10_sample_7
    Site_10_sample_8
    Site_11_sample_1

In contrast, the 15 extra metadata columns are inserted before this, and are
used to sort the samples:

.. code:: console

    $ cut -f 1,16 summary/with-metadata.ITS1.samples.onebp.tsv | head
    #Site  Sequencing sample
    01     Site_1_sample_1
    01     Site_1_sample_2
    01     Site_1_sample_3
    01     Site_1_sample_4
    01     Site_1_sample_5
    01     Site_1_sample_6
    01     Site_1_sample_7
    01     Site_1_sample_8
    01     Site_1_sample_9-2

Like the FASTQ filename stems, the metadata is still sorted as strings, but by
using leading zeros and ``YYYY-MM-DD`` style for any dates, you can achieve a
logical presentation.

After the sequencing sample name (the FASTQ filename stem), we have the
classification summary as a comma separated list - attempting to summarise
the later per-species columns. Species listed here with (*) are where
sequences matched multiple species equally well. For example, *Phytophthora
andina*, *P. infestans*, and *P. ipomoeae*, share an identical ITS1 marker.

The next columns are derived from the data itself, reads counts in the samples
as raw FASTQ, after read merging with Flash, primer trimming with Cutadapt,
the maximum ASV read count for non-spike-in or spike-in sequences, number of
singletons, total number of reads for the accepted ASVs (i.e. passing the
abundance threshold), and the number of unique ASVs accepted.
It may be easier to look at this in Excel, but at the command line:

.. code:: console

    $ cut -f 16,18-25 summary/with-metadata.ITS1.samples.onebp.tsv | head
    <SEE TABLE BELOW>

As a table:

================= ========= ===== ======== ============= ============ ========== ======== ======
Sequencing sample Raw FASTQ Flash Cutadapt Max non-spike Max spike-in Singletons Accepted Unique
================= ========= ===== ======== ============= ============ ========== ======== ======
Site_1_sample_1   6136      5900  5886     2269          0            692        4180     5
Site_1_sample_2   6135      5955  5947     2532          0            671        4548     6
Site_1_sample_3   6778      6484  6470     2146          0            579        5060     4
Site_1_sample_4   4145      3984  3974     1499          0            469        2852     4
Site_1_sample_5   4722      4232  4213     3130          0            433        3130     1
Site_1_sample_6   12633     12070 12034    5864          0            1217       9208     4
Site_1_sample_7   7560      7170  7141     3372          0            741        5402     5
Site_1_sample_8   6324      5956  5942     2037          0            630        4524     4
Site_1_sample_9-2 4542      4335  4331     2780          0            385        3436     2
================= ========= ===== ======== ============= ============ ========== ======== ======

Finally, we get to the main part of the sample table, one column per
classifier result, with the number of reads. Picking out some examples:

.. code:: console

    $ cut -f 16,28,40,60 summary/with-metadata.ITS1.samples.onebp.tsv | head
    <SEE TABLE BELOW>

As a table:

================= ======================== ========================= =======
Sequencing sample Phytophthora austrocedri Phytophthora gonapodyides Unknown
================= ======================== ========================= =======
Site_1_sample_1   165                      1158                      0
Site_1_sample_2   445                      718                       101
Site_1_sample_3   0                        1110                      1313
Site_1_sample_4   204                      861                       0
Site_1_sample_5   0                        3130                      0
Site_1_sample_6   0                        0                         0
Site_1_sample_7   0                        902                       161
Site_1_sample_8   0                        1863                      116
Site_1_sample_9-2 0                        0                         656
================= ======================== ========================= =======

Generally we hope to see single species predictions for each ASV, however when
there are conflicts such as equally good matches, or a reference sequence that
is shared between species, both are reported:

.. code:: console

    $ cut -f 16,29,30 summary/with-metadata.ITS1.samples.onebp.tsv | head
    <SEE TABLE BELOW>

As a table:

================= ====================== ===============================================
Sequencing sample Phytophthora cambivora Phytophthora cambivora;Phytophthora x cambivora
================= ====================== ===============================================
Site_1_sample_1   0                      182
Site_1_sample_2   0                      538
Site_1_sample_3   0                      0
Site_1_sample_4   0                      186
Site_1_sample_5   0                      0
Site_1_sample_6   0                      0
Site_1_sample_7   390                    0
Site_1_sample_8   0                      0
Site_1_sample_9-2 0                      0
================= ====================== ===============================================

In this example, while ``Site_1_sample_7`` had sequences uniquely matching
*Phytophthora cambivora*, ``Site_1_sample_1``, ``Site_1_sample_1`` and
``Site_1_sample_4`` instead had sequences which could be either *Phytophthora
cambivora* or *Phytophthora x cambivora*. These species are listed with a
``(*)`` suffix in the earlier classification summary column:

.. code:: console

    $ grep Site_1_sample_4 summary/with-metadata.ITS1.samples.onebp.tsv | cut -f 16,17
    Site_1_sample_4  Phytophthora austrocedri, Phytophthora cambivora(*), Phytophthora gonapodyides, Phytophthora pseudosyringae, Phytophthora x cambivora(*)
