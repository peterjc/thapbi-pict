Worked Example
==============

The *Quick Start* described a simplified use of the THAPBI PICT tool to
assess a single Illumina MiSeq sequencing run using the ``thapbi_pict
pipeline`` command, as a flowchart:

.. image:: ../images/pipeline.svg
   :alt: Flowchart summarising THAPBI PICT pipeline, from raw paired FASTQ files to reports.

Here we will run over the same process using real data, calling the individual
commands within the default pipeline - and include metadata for the reporting.
We will close with the equivalent all-in-one pipeline command.

.. image:: ../images/pipeline-meta.svg
   :alt: Flowchart summarising THAPBI PICT pipeline, from raw paired FASTQ files to reports, using metadata.

In these illustrative flow charts of the default pipeline, the input paired
FASTQ files (and metadata) are green, the intermediate per-sample FASTA and
TSV files are yellow, and the output reports are in orange. The individual
steps of the pipeline are dark blue boxes, and the ITS1 database is a pale
blue cylinder.

.. toctree::
   :maxdepth: 1
   :caption: Worked Example:

   sample_data
   prepare
   classify
   sample_summary
   read_summary
   edit_graph
