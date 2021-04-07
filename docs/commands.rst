Command Line
============

THAPBI PICT is a command line tool, meaning you must open your command line
terminal window and key in instructions to use the tool. The documentation
examples use the ``$`` (dollar sign) to indicate the prompt, followed by text
to be entered. For example, this should run the tool with no instructions:

.. code:: console

    $ thapbi_pict
    ...

Rather than literally printing dot dot dot, the tool should print out some
terse help, listing various sub-command names, and an example of how to get
more help.

For example, ``-v`` (minus sign, lower case letter v) or ``--version`` (minus,
minus, version in lower case) can be added to find out the version of the tool
installed:

.. code:: console

    $ thapbi_pict -v
    THAPBI PICT v0.8.1

THAPBI PICT follows the sub-command style popularised in bioinformatics by
``samtools`` (also used in the version control tool ``git``). This means most
of the instructions take the form ``thapbi_pict sub-command ...``, where the
dots indicate some additional options.

The main sub-commands are to do with classifying sequence files and reporting
the results, and these are described in the first :ref:`worked example
<woody_hosts>`:

* ``prepare`` - turn paired FASTQ input files for each sample, giving
  de-duplicated FASTA files
* ``classify`` - produce genus/species level predictions as
  tab-separated-variable TSV files
* ``summary`` - summarise a set of predictions by sample (with human readable
  report), and by unique sequence and sample (both with Excel reports)
* ``edit-graph`` - draw the unique sequences as nodes on a graph, connected by
  edit-distance
* ``pipeline`` - run all of the above in sequence
* ``assess`` - compare classifier output to known positive controls

There are further sub-commands to do with making or inspecting an SQLite3
format barcode marker sequence database, most of which are covered in the
second :ref:`worked example, with a custom database <recycled_water>`:

* ``dump`` - export a DB as TSV or FASTA format
* ``load-tax`` - import a copy of the NCBI taxonomy
* ``ncbi-import`` - import a FASTA file using the NCBI style naming
* ``curated-import`` - import a FASTA file where the descriptions are just the
  species names
* ``curated-seq`` - label prepared reads with known species assignment (single
  isolate positive controls)
* ``conflicts`` - report on genus or species level conflicts in the database

And some other miscellaneous commands:

* ``ena-submit`` - write a TSV table of your paired FASTQ files for use with
  the ENA interactive submission system.

Start with reading the help for any command using ``-h`` or ``--help`` as
follows:

.. code:: console

    $ thapbi_pict pipeline -h
    ...

Most of the commands have required arguments, and if you omit a required
argument it will stop with an error:

.. code:: console

    $ thapbi_pict pipeline
    ...
    thapbi_pict pipeline: error: the following arguments are required: -i/--input, -o/--output
