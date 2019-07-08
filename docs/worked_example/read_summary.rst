Read level report
=================

Running thapbi_pict read-summary
--------------------------------

The second report from the pipeline can be generated separately by the
``thapbi_pict read-summary`` command:

.. code::

    $ thapbi_pict read-summary -h
    ...

To mimick the filenames the pipeline command would use, we must set the
two output filenames explicitly, with ``-o`` or ``-output`` for the TSV
format table and ``-e`` or ``--excel`` for the Excel format table:

.. code::

    $ thapbi_pict read-summary -i intermediate/ \
      -o summary/thapbi-pict.reads.tsv -e summary/thapbi-pict.reads.xlxs
    ...

The contents of the two files are the same - here we will focus on the Excel
version which adds some visual formatting to improve usability.

This also command will accept the same :ref:`metadata <metadata>` arguments as
used earlier:

.. code::

    $ thapbi_pict read-summary -i intermediate/ \
      -o summary/with-metadata.samples.tsv -e summary/with-metadata.samples.xlsx \
      -t site_metadata.tsv -c 1,2,3,4,5,6,7,8,9,10,11,12,13,15 -x 16 -f 20
    ...

This will again affect the sort order of the sequences samples (here as
columns), and includes extra header rows containing the requested metadata.

Read Report
-----------

The heart of the read report is a large table, of unique sequences (rows)
versus sequenced samples (columns), with read abundance counts. There are
additional columns with sequence information, and when :ref:`metadata` is
present, extra rows at the start with sample information.

This read report has a row for each unique sequence. The first columns are
the unique sequence MD5 checksum, any species prediction, the sequence itself,
the number of samples it was detected in above the threshold, and the total
number of reads (from samples where it was above the threshold). Then the
main columns (one per sample) list the abundance of each unique sequence in
that sample (if above the threshold).

In the Excel version, conditional formatting is used to highlight the non-zero
counts with a red background. Furthermore, with metadata it will attempt to
assign repeated bands of background color to groups (pink, orange, yellow,
green, blue). In this example, each sample site gets a new color:

.. image:: https://user-images.githubusercontent.com/63959/60735578-ebdcf200-9f4b-11e9-8856-1ab66bd1245b.png
   :alt: Screenshot of Excel showing ``summary/with-metadata.samples.xlsx`` file.

Typical sample naming schemes will result in replicates as neighbouring
columns - meaning you should see very similar patterns of red (non-zero).
Certainly in this dataset scanning horizontally we do see some sequences
clearly show presence/absence patterns consistent with the samples.

The default row sorting will result in a dominant sequence being followed by
any close variants assigned to the same species. Many of these rows will
represent PCR artefacts found in just one or two samples. This contributes
to the "halo" effect seen in the :ref:`edit_graph` representation, discusssed
next.
