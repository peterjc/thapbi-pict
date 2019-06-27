Command Line
============

THAPBI PICT is a command line tool, meaning you must open your command line
terminal window and key in instructions to use the tool. The documentation
example use the ``$`` (dollar sign) to indicate the prompt, followed by text
to be entered. For example, this should run the tool with no instructions:

.. code:: bash

    $ thapbi_pict
    ...

Rather than literally printing dot dot dot, the tool should print out some
terse help, listing various sub-command names, and an example of how to get
more help.

For example, ``-v`` (minus sign, lower case letter v) or ``--version`` (minus,
minus, version in lower case) can be added to find out the version of the tool
installed:

.. code:: bash

    $ thapbi_pict -v
    THAPBI PICT v0.3.1   

THAPBI PICT follows the sub-command style popularised in bioinformatics by
``samtools`` (also used in the version control tool ``git``). This means most
of the instructions take the form ``thapbi_pict sub-command ...``, where the
dots indicate some additional options. Start with reading the help for any
command using ``-h`` or ``--help`` as follows:

.. code:: bash

    $ thapbi_pict pipeline -h
    ...

The main sub-commands are to do with classifying sequence files and reporting
the results:

* ``prepare`` - turn paired FASTQ input files for each sample, giving de-duplicated FASTA files
* ``classify`` - produce genus/species level predictions as tab-separated-variable TSV files
* ``sample-summary`` - summarise a set of predictions by sample (with human readable report)
* ``read-summary`` - summarise a set of predictions by unique sequence and sample (with Excel report)
* ``edit-graph`` - draw the unique sequences as nodes on a graph, connected by edit-distance
* ``pipeline`` - run all of the above in sequence

There are further sub-commands to do with making or inspecting an SQLite3
format ITS1 database:

* ``dump`` - export a DB as TSV or FASTA format
* ``load-tax`` - import a copy of the NCBI taxonomy
* ``ncbi-import`` - import a FASTA file using the NCBI style naming
* ``import-seq`` - import prepared reads with known species assignment (single isolate positive controls)
* ``legacy-import`` - import a FASTA file using the style of our legacy curated ITS1 dataset
