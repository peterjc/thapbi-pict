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
command as follows, and should get five output report files:

.. code:: console

    $ thapbi_pict pipeline -i raw_data/ -s intermediate/ -o summary/
    ...
    $ ls -1 summary/thapbi-pict.*
    thapbi-pict.reads.tsv
    thapbi-pict.reads.xlsx
    thapbi-pict.samples.tsv
    thapbi-pict.samples.txt
    thapbi-pict.edit-graph.xgmml

As described for the :ref:`prepare-reads step <prepare_reads>` we should also
specify which of the samples are negative controls, which may be used to
increase the plate level minimum abundance threshold:

.. code:: console

    $ thapbi_pict pipeline -i raw_data/ -s intermediate/ -o summary/ \
      -n raw_data/NEGATIVE*.fastq.gz
    ...

And, as described for the :ref:`sample-summary <sample_summary>` and
:ref:`read-summary <read_summary>` steps, we can provide metadata:

.. code:: console

    $ thapbi_pict pipeline -i raw_data/ -s intermediate/ -o summary/ \
      -n raw_data/NEGATIVE*.fastq.gz -r with-metadata \
      -t site_metadata.tsv -c 1,2,3,4,5,6,7,8,9,10,11,12,13,15 -x 16 -f 20
    ...
    $ ls -1 summary/with-metadata.*
    with-metadata.reads.tsv
    with-metadata.reads.xlsx
    with-metadata.samples.tsv
    with-metadata.samples.txt
    with-metadata.edit-graph.xgmml

Here we also used ``-r`` (or ``--report``) to specify a different stem
for the report filenames.


