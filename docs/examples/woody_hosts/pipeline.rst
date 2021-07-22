Pipeline with metadata
======================

Running thapbi-pict pipeline
----------------------------

Having run all the steps of the typical pipeline individually, we now return
to the top level ``thapbi_pict pipeline`` command:

.. code:: console

    $ thapbi_pict pipeline -h
    ...

Assuming you have the FASTQ files in ``raw_data/``, we could run the pipeline
command as follows, and should get multiple output report files:

.. code:: console

    $ thapbi_pict pipeline -i raw_data/ -s intermediate/ \
      -o summary/thapbi-pict
    ...
    $ ls -1 summary/thapbi-pict.*
    summary/thapbi-pict.ITS1.all_reads.fasta
    summary/thapbi-pict.ITS1.all_reads.onebp.tsv
    summary/thapbi-pict.ITS1.reads.onebp.tsv
    summary/thapbi-pict.ITS1.reads.onebp.xlsx
    summary/thapbi-pict.ITS1.samples.onebp.tsv
    summary/thapbi-pict.ITS1.samples.onebp.txt
    summary/thapbi-pict.ITS1.samples.onebp.xlsx
    summary/thapbi-pict.edit-graph.onebp.pdf
    summary/thapbi-pict.edit-graph.onebp.xgmml

As described for the :ref:`prepare-reads step <prepare_reads>` we should also
specify which of the samples are negative controls, which may be used to
increase the plate level minimum abundance threshold:

.. code:: console

    $ thapbi_pict pipeline -i raw_data/ -s intermediate/ \
      -o summary/thapbi-pict -n raw_data/NEGATIVE*.fastq.gz
    ...

And, as described for the :ref:`summary reports <summary_reports>`, we can
provide metadata:

.. code:: console

    $ thapbi_pict pipeline -i raw_data/ -s intermediate/ \
      -o summary/with-metadata -n raw_data/NEGATIVE*.fastq.gz \
      -t metadata.tsv -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -x 16 -f 20
    ...

Finally, as we will review next, we can ask the pipeline to assess the results
against any expected sample species classifications:

.. code:: console

    $ thapbi_pict pipeline -i raw_data/ expected/ -s intermediate/ \
      -o summary/with-metadata -n raw_data/NEGATIVE*.fastq.gz \
      -t metadata.tsv -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -x 16 -f 20
    ...
    $ ls -1 summary/with-metadata.*
    summary/with-metadata.ITS1.all_reads.fasta
    summary/with-metadata.ITS1.all_reads.onebp.tsv
    summary/with-metadata.ITS1.assess.confusion.onebp.tsv
    summary/with-metadata.ITS1.assess.onebp.tsv
    summary/with-metadata.ITS1.assess.tally.onebp.tsv
    summary/with-metadata.ITS1.reads.onebp.tsv
    summary/with-metadata.ITS1.reads.onebp.xlsx
    summary/with-metadata.ITS1.samples.onebp.tsv
    summary/with-metadata.ITS1.samples.onebp.txt
    summary/with-metadata.ITS1.samples.onebp.xlsx

Here we also used ``-r`` (or ``--report``) to specify a different stem
for the report filenames.

Conclusions
-----------

For the THAPBI Phyto-Threats project our datasets span multiple plates, but we
want to set plate-specific minimum abundance thresholds. That is taken care of
as long as each plate is in its own directory. For example, you might have
``raw_data/plate_NNN/*.fastq.gz`` and run the pipeline with ``-i raw_data/``).

However, while you could run the pipeline command on all the data in one go,
with access to a computer cluster it will likely be faster to run at least the
(slowest)  ``prepare-reads`` stage on separate cluster nodes (e.g. one cluster
job for each plate).
