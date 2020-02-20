.. _sample_summary:

Sample level report
===================

Running thapbi_pict sample-summary
----------------------------------

The first output reports from the pipeline can be generated separately by the
``thapbi_pict sample-summary`` command:

.. code:: console

    $ thapbi_pict sample-summary -h
    ...

To mimic what the pipeline command would do, run the following, where we have
to give names to the computer readable TSV (``-o`` or ``--output``), equivalent
Excel (``-e`` or ``--excel``), and human readable plain text TXT outputs
(``-r`` or ``--report``):

.. code:: console

    $ thapbi_pict sample-summary -i intermediate/ \
      -o summary/thapbi-pict.samples.onebp.tsv \
      -e summary/thapbi-pict.samples.onebp.xlxs \
      -r summary/thapbi-pict.samples.onebp.txt
    ...

Note the trailing slash ``\`` at the end of the first line indicates the
command continues on the next line. You can actually type this at the standard
Linux command prompt (or include it in a copy and paste), or just enter this
as one very long command.

We will look at the output in a moment, along side the equivalent reports
generated with :ref:`metadata <metadata>` (see linked discussion about column
and row numbers):

.. code:: console

    $ time thapbi_pict sample-summary -i intermediate/ \
      -o summary/with-metadata.samples.onebp.tsv \
      -e summary/with-metadata.samples.onebp.xlsx \
      -r summary/with-metadata.samples.onebp.txt \
      -t metadata.tsv -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -x 16 -f 20
    ...

This time the even longer command has been shown split over five lines.

The computer readable TSV and equivalent Excel file will include the metadata
as additional leading columns, and also differ in the line order.

We next focus on the changes in the human readable text file.

Sample Report
-------------

Here we will discuss the high level human readable summary report from
``thapbi_pict sample-summary``, produced as plain text.

Without :ref:`metadata <metadata>`, the samples are sorted by filename
alphabetically. In this example that means we get DNA controls, negative
controls, Site 1, 10, 11, 12, 13, 14, 2, ..., 9. This is unfortunate, so if it
is too late to change your sequence sample naming scheme (e.g. use leading
zeros, and ``YYYY-MM-DD`` style for any dates), you can at least use nicely
sorting names in your metadata.

When a metadata table is given, the rows are sorted by the displayed columns
(in the order requested), with any sequenced files without metadata entries
shown at the very end. Thus we get site ``01`` to ``15``, ``DNA10MIX``,
``DNA16MIX`` and ``NEGATIVE`` last. Within site ``01``, we get the sequenced
samples in the order given in column 16, i.e. ``Site_1_sample_1``,
``Site_1_sample_2``, ..., ``Site_1_sample_10`` as desired.

Pulling out the ``Site_1_sample_1`` example, without metadata in file
``summary/thapbi-pict.samples.onebp.txt`` we see:

.. code:: text

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
``Site_1_sample_2`` etc. The alphabetic sort order problem again.

As to the meaning of this list, those are the species identified - in some
cases with a caveat. The file starts with a tiny explanation:

.. code:: text

    NOTE: Species listed with (uncertain/ambiguous) in brackets are where
    sequences matched multiple species equally well. For example,
    *Phytophthora andina*, *P. infestans*, and *P. ipomoeae*, share an
    identical marker.

In this case, as you may recall from when we looked at the classifier output
file ``intermediate/Site_1_sample_1.onebp.tsv``, one of the sequences matched
both *Phytophthora cambivora* and *Phytophthora x cambivora*.

In comparison, with metadata in file ``summary/with-metadata.samples.onebp.txt``,
all the samples matched to a metadata row get grouped with a shared metadata
header:

.. code:: text

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

Note that for ``Site_1_sample_2``, at least one unique sequence was not given
a species or even genus level classification, thus the ``Unknown`` entry. This
likely reflects a gap in the default database, and/or the default method being
too strict.

If any of the requested metadata is missing (i.e. a blank entry in the table
for a requested field), then it does not get shown. For example, this applies
to the DNA mixes and the negative controls.
