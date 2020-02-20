.. _metadata:

Metadata
========

The :ref:`quick_start` introduced the typical pipeline taking paired FASTQ
files though to reports, and mentioned the idea of enhancing the reports with
sample metadata.

.. image:: ../../images/pipeline-meta.svg
   :alt: Flowchart summarising THAPBI PICT pipeline, from raw paired FASTQ files to reports, using metadata.

In the following we will show the reports with and without metadata.
File ``metadata.tsv`` is a table of metadata (based on table S1 in the
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
separated). They will be shown in this order.

This cross referencing idea is key to getting the best results from attaching
metadata to your sequenced samples. Here is an abridged representation of the
table, showing column one (site or control name), column two (altitute), and
finally column 16 which has the filename stems of the sequence data belonging
to this row of the table (semi-colon separated list).

======== ======== === ========================================================
#Site    Altitude ... MiSeq Sample(s)
======== ======== === ========================================================
01             30 ... Site_1_sample_1; Site_1_sample_2; Site_1_sample_3; Site_1_sample_4; Site_1_sample_5; Site_1_sample_6; Site_1_sample_7; Site_1_sample_8; Site_1_sample_9-2; Site_1_sample_10
02             55 ... Site_2_sample_1; Site_2_sample_2; Site_2_sample_3; Site_2_sample_4; Site_2_sample_5; Site_2_sample_6; Site_2_sample_7; Site_2_sample_8; Site_2_sample_9; Site_2_sample_10
03             45 ... Site_3_sample_1; Site_3_sample_2; Site_3_sample_4; Site_3_sample_7; Site_3_sample_8; Site_3_sample_9
04             20 ... Site_4_sample_1; Site_4_sample_2; Site_4_sample_3; Site_4_sample_3-2; Site_4_sample_4; Site_4_sample_5; Site_4_sample_6; Site_4_sample_8; Site_4_sample_9; Site_4_sample_10
05            100 ... Site_5_sample_1; Site_5_sample_2; Site_5_sample_4; Site_5_sample_5; Site_5_sample_6; Site_5_sample_8; Site_5_sample_9
06              5 ... Site_6_sample_1; Site_6_sample_2-2; Site_6_sample_3-1; Site_6_sample_4; Site_6_sample_5-3; Site_6_sample_6; Site_6_sample_7-1; Site_6_sample_8-2; Site_6_sample_9; Site_6_sample_10
07            105 ... Site_7_sample_1; Site_7_sample_2; Site_7_sample_3; Site_7_sample_5; Site_7_sample_6; Site_7_sample_7; Site_7_sample_8; Site_7_sample_9; Site_7_sample_10
08             45 ... Site_8_sample_1; Site_8_sample_2; Site_8_sample_3; Site_8_sample_4; Site_8_sample_5-2; Site_8_sample_6; Site_8_sample_7; Site_8_sample_7-2; Site_8_sample_8; Site_8_sample_9
09             15 ... Site_9_sample_1; Site_9_sample_4-3; Site_9_sample_6; Site_9_sample_7; Site_9_sample_8; Site_9_sample_9; Site_9_sample_10
10             30 ... Site_10_sample_7; Site_10_sample_8
11             80 ... Site_11_sample_1; Site_11_sample_2; Site_11_sample_3; Site_11_sample_4; Site_11_sample_5; Site_11_sample_6; Site_11_sample_7; Site_11_sample_8; Site_11_sample_9; Site_11_sample_10
12             30 ... Site_12_sample_1; Site_12_sample_2; Site_12_sample_3-3; Site_12_sample_4; Site_12_sample_5-3; Site_12_sample_6; Site_12_sample_8; Site_12_sample_9; Site_12_sample_10
13            300 ... Site_13_sample_1; Site_13_sample_2; Site_13_sample_4; Site_13_sample_5; Site_13_sample_6; Site_13_sample_7; Site_13_sample_8; Site_13_sample_9; Site_13_sample_10
14             30 ... Site_14_sample_1-2; Site_14_sample_2; Site_14_sample_3; Site_14_sample_4; Site_14_sample_5; Site_14_sample_6; Site_14_sample_10
DNA10MIX          ... DNA10MIX_undiluted; DNA10MIX_diluted25x; DNA10MIX_bycopynumber
DNA16MIX          ... DNA16MIX
NEGATIVE          ... NEGATIVE_firstplate; NEGATIVE_secondplate
======== ======== === ========================================================

Also note that in column one we have listed the numerical site names with
leading zeros giving ``01`` to ``14`` to ensure they sort as expected.
