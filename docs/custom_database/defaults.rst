Pipeline with defaults
======================

Running thapbi-pict pipeline
----------------------------

First, we will run the THAPBI PICT pipeline command with largely default
settings, other than including the metadata about the water samples. Note
that this dataset has no blanks or negative controls, so must trust to the
default minimum abundance threshold.

The key values which we will be changing later are the primers and database.

Assuming you have the FASTQ files in ``raw_data/``, run the pipeline command
as follows, and you should get five output report files:

.. code:: console

    $ mkdir intermediate_defaults/ summary_defaults/
    $ thapbi_pict pipeline \
      -i raw_data/ -o summary_defaults/ -s intermediate_defaults/ \
      -r recycled-water-defaults -t metadata.tsv -c 1,2,3,4,5,6,7
    ...
    $ ls -1 intermediate_defaults/SRR*.fasta | wc -l
    384
    $ ls -1 intermediate_defaults/SRR*.onebp.tsv | wc -l
    384
    $ ls -1 summary/summary_defaults/thapbi-pict.*
    recycled-water-defaults.reads.tsv
    recycled-water-defaults.reads.xlsx
    recycled-water-defaults.samples.tsv
    recycled-water-defaults.samples.txt
    recycled-water-defaults.edit-graph.xgmml

Here we used ``-r`` (or ``--report``) to specify a different stem for the
report filenames. The metadata options were described earlier in -- this is
perhaps an idealised example in that ``metadata.tsv`` was created so that
we could just all seven columns of the table, and the sample name (filename
prefix) is in the first column.


Results
-------

Will pick a couple of samples to compare and contrast with the second run...
