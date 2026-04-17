.. _custom_database_pipeline:

Pipeline with custom database
=============================

Running thapbi-pict pipeline
----------------------------

Having created the dual-marker database ``pooled.sqlite``, running the pipeline
is quite straightforward - but I have modified a lot of the default settings:

.. code:: console

    $ mkdir -p intermediate/ summary/
    $ thapbi_pict pipeline -d pooled.sqlite --synthetic '' \
        -m 1s3g --denoise unoise-l --unoise_alpha 6 \
        -i raw_data/ expected/ \
        --merged-cache tmp_merged/ \
        -s intermediate/ -o summary/ \
        -t metadata.tsv -x 1 -c 3,7,4,6
    ...
    Running 1s3g classifier on summary/rps10.tally.tsv, 186 sequences
    1s3g classifier assigned species/genus to 62 of 186 unique sequences from 1 files
    ...

In addition to setting the database ``-d pooled.sqlite`` we also use
``--synthetic ''`` tell the pipeline this dataset lacks any synthetic controls.
Minimal denoising has been enabled with ``--denoise unoise-l --unoise_alpha 6``
to remove the long tail of unclassified unique amplicon sequence variants.

The classifier method was set explicitly with ``-m 1s3g`` to allow fuzzier
matching than the more conservative default reflecting the relatively sparse
rps10 database coverage.

Results
-------

Pending.
