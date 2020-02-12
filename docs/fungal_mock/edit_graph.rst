Edit Graphs
===========

The sequence edit-graph is very useful for understanding what came off the
sequencer - although you may need to play with the thresholds to find a sweet
spot for hiding the noise. If you have loaded the XGMML network file into
Cytoscape you can do this interactively by selecting nodes based on the
``Max-sample-abundance`` attribute and hiding or removing them. This will
not update the ``Sample-count`` and node sizes. For that you can re-run
``thapbi_pict edit-graph`` with the higher sample level minimum abundance
setting (``-a`` or ``--abundance``). You do not need to regenerate the
intermediate per-sample FASTA files unless you want to use a lower threshold.

The following figures are from the example script ``run.sh`` which called
``thapbi_pict edit-graph`` with ``-a 75``, meaning a unique sequence had to be
in a sample from at least 75 reads to be considered. Using a lower value gives
a much noiser picture (see the halo effect discussed earlier).

Additionally this used ``-s`` (or ``--showdb``) to force including all of the
database sequences (dark red nodes), as some did not appear in the samples
(shown as the smallest dark red dots, typically the bottom row of the image).

Amplicon library one - ITS1
---------------------------

Starting with amplicon library one, where the BITS/B58S3 primers were used for
a short fragment of ITS1.

.. image:: ../images/amp_lib_one.BITS_B58S3.edit-graph.a75.svg
   :target: https://raw.githubusercontent.com/peterjc/thapbi-pict/master/docs/images/amp_lib_one.BITS_B58S3.edit-graph.a75.svg?sanitize=true
   :alt: Sequence edit-graph for amplicon library one using BITS/B58S3 primers for ITS1.

This is from file ``amp_lib_one.BITS_B58S3.edit-graph.a75.xgmml`` created by
``run.sh``.

The dark red nodes represent sequences in the database - given how this
database was constructed to match the mock community, we would hope to see all
the database entries represented in the samples. Some are missing at this
abundance threshold (four bottom row entries *Saccharomyces cerevisiae*,
*Ustilago maydis*, *Rhizomucor miehei* and *Chytriomyces hyalinus*, plus one
of the multiple *Rhizomucor irregularis* entries shown middle right).

The large red nodes are the well represented community members, starting with
*Naganishia albida* shown top left, which has four different 1bp variants some
of which appear in many samples - you can see the sample counts in you load
the XGMML file for this graph in CytoScape. These are common enough to suggest
they could be alternative versions of the ITS1 region in the genomes of these
community members?

The large grey nodes not connected to a red node represent unwanted reads,
likely contaminations discussed later.

In general each species is reprented by a single connected component. The
exceptions are *Rhizomucor irregularis* (multiple distantly related entries)
and the *Fusarium*. The expected sequence for *Fusarium verticillioides* was
not seen at all, however there are a great many copies one base away from
the expected *Fusarium oxysporum* sequence (abbreviated MD5 checksum
``bb28f2``, in full ``bb28f2b57f8fddefe6e7b5d01eca8aea``). Is this perhaps
coming from the *Fusarium verticillioides* strain?

Amplicon library two - ITS1
---------------------------

First, analysed using the same BITS/B58S3 primers as for ITS1 as in amplicon
library one - the unique sequence MD5 checksums overlap with those seen in
amplicon one:
Notice the presence/absense is different:

.. image:: ../images/amp_lib_two.BITS_B58S3.edit-graph.a75.svg
   :target: https://raw.githubusercontent.com/peterjc/thapbi-pict/master/docs/images/amp_lib_two.BITS_B58S3.edit-graph.a75.svg?sanitize=true
   :alt: Sequence edit-graph for amplicon library two using BITS/B58S3 primers for ITS1 (although actually amplified with ITS1f/ITS2 primers).

This is from file ``amp_lib_two.BITS_B58S3.edit-graph.a75.xgmml`` created by
``run.sh``.

Broadly the same as from amplicon library one, but notice the
presence/absense patterns are different. Also there are more variants of the
``bb28f2`` *Fusarium*, and a pair of unexpected grey nodes 3bp apart
(``e055cb`` and ``ee5482``, which as discussed later match *glomeromycetes*,
perhaps *Rhizophagus* which was in the mock community).

Now, using the actual primer pair, ITS1f/ITS2, which give a longer ITS1
fragment. Note that the sequences are extended so the checksums are different
to those in the preceeding images, but again broadly the same picture as the
two images above:

.. image:: ../images/amp_lib_two.ITS1f_ITS2.edit-graph.a75.svg
   :target: https://raw.githubusercontent.com/peterjc/thapbi-pict/master/docs/images/amp_lib_two.ITS1f_ITS2.edit-graph.a75.svg?sanitize=true
   :alt: Sequence edit-graph for amplicon library two using ITS1f/ITS2 primers for ITS1.

This is from file ``amp_lib_two.ITS1f_ITS2.edit-graph.a75.xgmml`` created by
``run.sh``.

Amplicon library two - ITS2
---------------------------

Finally, amplicon library two using the ITS3-KYO and ITS4-KYO3 primers for ITS2.

.. image:: ../images/amp_lib_two.ITS3-KYO2_ITS4-KYO3.edit-graph.a75.svg
   :target: https://raw.githubusercontent.com/peterjc/thapbi-pict/master/docs/images/amp_lib_two.ITS3-KYO2_ITS4-KYO3.edit-graph.a75.svg?sanitize=true
   :alt: Sequence edit-graph for amplicon library two using ITS3-KYO and ITS4-KYO3 primers for ITS2.

This is from file ``amp_lib_two.ITS3-KYO2_ITS4-KYO3.edit-graph.a75.xgmml``
created by ``run.sh``.

Some more noteworthy changes to presence/absense, including much more
*Saccharomyces cerevisiae* (still drawn bottom left). Also there are no
unexpected grey nodes, and perhaps most interestingly from a species
classification point of view, now the three *Fusarium* species fall into
separate connected components.
