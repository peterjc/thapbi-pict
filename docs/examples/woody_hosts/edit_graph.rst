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

    $ thapbi_pict edit-graph -i intermediate/ -o summary/thapbi-pict.edit-graph.onebp.xgmml
    ...

This will generate an XGMML (eXtensible Graph Markup and Modeling Language)
file by default, but you can also request other formats including PDF
(which requires additional dependencies including GraphViz):

.. code:: console

    $ thapbi_pict edit-graph -i intermediate/ -o summary/thapbi-pict.edit-graph.onebp.pdf -f pdf
    ...

.. WARNING:

    With larger datasets, the edit graph is easily the slowest of the report
    commands, and the PDF output even more so.

Nodes and edges
---------------

In this context, we are talking about a graph in the mathematical sense of
nodes connected by edges. Our nodes are unique sequences (which we can again
label by the MD5 checksum), and the edges are how similar two sequences are.
Specially, we are using the Levenshtein edit distance. This means an edit
distance of one could be a single base substitution, insertion or deletion.

The tool starts by compiling a list of all the unique sequences in your
samples (i.e. all the rows in the ``thapbi_pict read-summary`` report), and
optionally all the unique sequences in the database. It then computes the
edit distance between them all (this can get slow).

We build the network graph by adding edges for edits of up to three base pairs
(by default). This gives small connected components or sub-graphs which are
roughly at the species level.

Redundant edges are dropped, for example if *A* is one edit away from *B*,
and *B* is one edit away from *C*, there is need to draw the two edit line
from *A* to *C*.

We draw the nodes as circles, scaled by the number of samples that unique
sequence appeared in. If that exact sequence is in the database, is it colored
according to genus, defaulting to grey.

=========== ========== ===================
Color       RGB value  Meaning
=========== ========== ===================
Red         ``FF0000``  *Phytophthora*
Lime        ``00FF00``  *Peronospora*
Blue        ``0000FF``  *Hyaloperonospora*
Yellow      ``FFFF00``  *Bremia*
Cyan        ``00FFFF``  *Pseudoperonospora*
Magenta     ``FF00FF``  *Plasmopara*
Maroon      ``800000``  *Nothophytophthora*
Olive       ``808000``  *Peronosclerospora*
Green       ``008000``  *Perofascia*
Purple      ``800080``  *Paraperonospora*
Teal        ``008080``  *Protobremia*
Dark red    ``8B0000``  Other known genus
Dark orange ``FF8C00``  Conflicting genus
Orange      ``FFA500``  Synthetic sequence
Grey        ``808080``  Not in the database
=========== ========== ===================

The edges are all grey, solid for a one base pair edit distance, dashed for a
two base pair edit distance, and dotted for a three base pair edit distance.

Viewing the PDF
---------------

You should be able to open the PDF file easily, and get something like this -
lots of red circles for *Phytophthora*, some grey circles for sequences not
in the database, and plenty of grey straight line edges between them.

.. image:: https://user-images.githubusercontent.com/63959/94339236-3551ab00-fff0-11ea-96cf-d30ef04fbaf3.png

In the PDF (and XGMML) output, nodes are coloured by genus (red for
*Phytophthora*), but only labelled if in the database at species level.

The edges are solid for a one base pair edit distance, dashed for a two base
pair edit distance, and dotted for a three base pair edit distance. All grey.

Viewing the XGMML
-----------------

You should be able to open the PDF file easily, and while it is interesting
it is read only and non-interactive. This is where the XGMML output shines.
You will need to install the free open source tool  `Cytoscape
<https://cytoscape.org/>`_ to use this.

Open Cytoscape, and from the top level menu select ``File``, ``Import``,
``Network from file...``, then select
``summary/thapbi-pict.edit-graph.onebp.xgmml`` (the XGMML file created above).

You should get something like this, where initially all the nodes are drawn
on top of each other:

.. image:: https://user-images.githubusercontent.com/63959/94338640-d2f6ab80-ffeb-11ea-9256-760b6f0dfe19.png

From the top level menu select "Layout", "Perfuse Force Directed Layout",
"Edit-distance-weight", and you should then see something prettier - if
you zoom in you should see something like this:

.. image:: https://user-images.githubusercontent.com/63959/94338441-1819de00-ffea-11ea-82b8-ef1b83dd03d9.png

This time you can interact with the graph, moving nodes about with the mouse,
try different layouts, view and search the attributes of the nodes and edges.

Here the nodes are labelled with the species if they were in the database
at species level, or otherwise as the start of the MD5 checksum in curly
brackets (so that they sort nicely). The default node colors are as in the
PDF output, likewise the grey edge styles.

The node attributes include the full MD5 (so you can lookup the full sequence
or classification results for any node of interest), sample count, total read
abundance (both numbers shown in the ``thapbi_pict read-summary`` report),
genus (allowing you to do your own color scheme), and species if known.

The edge attributes include ``Edit-distance`` (values ``1``, ``2``, ``3``
for number of base pairs difference between sequences) and matching
``Edit-distance-weight`` (values ``3``, ``2``, ``1`` used earlier for the
layout where we prioritise the small edit distance edges).

.. _halo_effect:

Halo effect
-----------

In this final screenshot we have zoomed in and selected all 11 nodes in the
connected component centered on *P. pseudosyringae* (Cytoscape highlights
selected nodes in yellow):

.. image:: https://user-images.githubusercontent.com/63959/94338444-1b14ce80-ffea-11ea-8cde-cc4971ba5853.png

The node table view is automatically filtered to show just these nodes, and we
can see that all the grey nodes appeared in only one sample each - with the
*P. pseudosyringae* entry in the database in 66 samples, while the one base
away *P. ilics* sequence was in 6 samples.

This kind of grey-node halo around highly abundance sequences is more common
when plotting larger datasets. It is consistent with PCR artefacts occurring
in just one (or two) samples giving rise to (almost) unique sequences based on
the template sequence.
