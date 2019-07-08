.. _edit_graph:

Edit Graph
==========

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

Viewing the PDF
---------------

You should be able to open the PDF file easily, and get something like this:

.. image:: https://user-images.githubusercontent.com/63959/60818082-d7d80100-a194-11e9-8002-7d855c0c0bd1.png

In this representation, graph nodes (circles) represent unique sequences, and
graph edges (straight lines) represent the edit distances between them.

In the PDF output, nodes are coloured by genus (red for *Phytophthora*), and
labelled if in the database at pecies level. Otherwise they are grey.

The edges are solid for a one base pair edit distance, dashed for a two base
pair edit distance, and dotted for a three base pair edit distance. All grey.

Viewing the XGMML
-----------------

You should be able to open the PDF file easily, and while it is interesting
it is read only and non-interactive. This is where the XGMML output shines.
You will need to install the free open source tool  `Cytoscape
<https://cytoscape.org/>`_ to use this.

Open Cytoscape, and from the top level menu select ``File``, ``Import``,
``Network from file...``, then select ``summary/thapbi-pict.edit-graph.xgmml``
(the XGMML file created above).

You should get something like this, were initially all the nodes are drawn
on top of each other:

.. image:: https://user-images.githubusercontent.com/63959/60818958-8af52a00-a196-11e9-9d89-27f5d027f1e7.png

From the top level menu select "Layout", "Perfuse Force Directed Layout",
"Edit-distance-weight", and you should then see something prettier - if
you zoom in you should see something like this:

.. image:: https://user-images.githubusercontent.com/63959/60818965-8d578400-a196-11e9-87ff-467da6fa0353.png

This time you can interact with the graph, moving nodes about with the mouse,
try different layouts, view and search the attributes of the nodes and edges.

Here the nodes are labelled with the species if they were in the database
at species level, or otherwise as the start of the MD5 checkum in curly
brackets (so that they sort nicely). The default node colors as as in the
PDF output, likewise the edge styles.

The node attributes include the full MD5 (so you can lookup the full sequence
or classification results for any node of interest), sample count, total read
abundance (both numbers shown in the ``thapbi_pict read-summary`` report),
genus (allowing you to do your own color scheme), and species if known.

The edge attributes include ``Edit-distance`` (values ``1``, ``2``, ``3``
for number of base pairs difference between sequences) and matching
``Edit-distance-weight`` (values ``3``, ``2``, ``1`` used earlier for the
layout where we prioritise the small edit distance edges).

Halo effect
-----------

In this final screenshot we have zoomed in and selected all 11 nodes in the
connected component centered on *P. pseudosyringae* (Cytoscape highlights
selected nodes in yellow):

.. image:: https://user-images.githubusercontent.com/63959/60819693-f4296d00-a197-11e9-9605-9eb573666f37.png

The node table view is automatically filtered to show just these nodes, and we
can see that all the grey nodes appeared in only one sample each - with the
*P. pseudosyringae* entry in the database in 66 samples, while the one base
away *P. ilics* sequence was in 6 samples.

This kind of grey-node halo around highly abundance sequences is more common
when plotting larger datasets. It is consistent with PCR artefacts occuring in
just one (or two) samples giving rise to (almost) unique sequences based on
the template sequence.
