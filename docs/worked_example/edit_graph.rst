Running thapbi_pict edit-graph
------------------------------

The final step of the pipeline command can be run alone as ``thapbi_pict
edit-graph``:

.. code:: console

    $ thapbi_pict edit-graph -h
    ...

This command does not use the intermediate TSV files or metadata, just the
intermediate FASTA files and the ITS1 database.

To mimic the pipeline output, we must set the output filename explicitly
with ``-o`` or ``--output``:

.. code:: console

    $ thapbi_pict edit-graph -i intermediate/ -o summary/thapbi-pict.edit-graph.xgmml
    ...

This will generate an XGMML (eXtensible Graph Markup and Modeling Language)
file by default, but you can also request other formats including PDF
(which requires additional dependencies including GraphViz):

.. code:: console

    $ thapbi_pict edit-graph -i intermediate/ -o summary/thapbi-pict.edit-graph.pdf -f pdf
    ...

.. WARNING:

    With larger datasets, the edit graph easily the slowest of the report
    commands, and the PDF output even more so.

Edit Graph
----------

You should be able to open the PDF file easily, and while it is interesting
it is read only and non-interactive. This is where the XGMML output shines.
You will need to install the free open source tool  `Cytoscape
<https://cytoscape.org/>`_ to use this..
