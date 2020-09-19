.. _woody_hosts:

Environmental *Phytophthora* ITS1
=================================

The :ref:`quick_start` described a simplified use of the THAPBI PICT tool to
assess a single Illumina MiSeq sequencing run using the ``thapbi_pict
pipeline`` command, as a flowchart:

.. image:: ../../images/pipeline-meta.svg
   :alt: Flowchart summarising THAPBI PICT pipeline, from raw paired FASTQ files to reports, using metadata.

Here we will run over the same process using real *Phytophthora* ITS1 data,
calling the individual commands within the default pipeline - and include
metadata for reporting. We then run the equivalent all-in-one pipeline command.

Finally, since the sample data includes some positive controls, we can look at
assessing the classifier performance.

.. toctree::
   :maxdepth: 1

   sample_data
   prepare
   classify
   metadata
   summary
   edit_graph
   pipeline
   assess
