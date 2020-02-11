Edit Graphs
===========

The sequence edit-graph is very useful for understanding that came off the
sequencer - although you may need to play with the thresholds to find a sweet
spot for hiding the noise. If you have run the pipeline (or at least the
prepare-reads step), you can re-run ``thapbi_pict edit-graph`` with a higher
sample level minimum abundance (``-a`` or ``--abundance``).

The following figures are from the example script ``run.sh`` which called
``thapbi_pict edit-graph` with ``-a 75``, meaning a unique sequence had to be
in a sample from at least 75 reads to be considered. Using a lower value gives
a much noiser picture (see the halo effect discussed earlier).

Additionally this used ``-s`` (or ``--showdb``) to force including all of the
database sequences (dark red nodes), as some did not appear in the samples
(shown as the smallest dark red dots, typically the bottom row of the image).

Amplicon library one - ITS1
---------------------------

Starting with amplicon library one, where the BITS/B58S3 primers we used for
a short fragment of ITS1:

.. image:: ../images/amp_lib_one.BITS_B58S3.edit-graph.a75.svg
   :alt: Sequence edit-graph for amplicon library one using BITS/B58S3 primers for ITS1. Minimum abundance threshold 75.

Amplicon library two - ITS1
---------------------------

First, analysed using the same BITS/B58S3 primers as for ITS1 as in amplicon
library one - the unique sequence MD5 checksums overlap with those seen in
amplicon one. Notice the presence/absense is different:

.. image:: ../images/amp_lib_two.BITS_B58S3.edit-graph.a75.svg
   :alt: Sequence edit-graph for amplicon library two using BITS/B58S3 primers for ITS1 (although actually amplified with ITS1f/ITS2 primers). Minimum abundance threshold 75.

Now, using the *actual* primer pair, ITS1f/ITS2 which give a longer ITS1
fragment. Note that the sequences are extended so the checksums are different
to those in the preceeding images:


.. image:: ../images/amp_lib_two.ITS1f_ITS2.edit-graph.a75.svg
   :alt: Sequence edit-graph for amplicon library two using ITS1f/ITS2 primers for ITS1. Minimum abundance threshold 75.

Amplicon library two - ITS2
---------------------------

Finally, amplicon library two using the ITS3-KYO and ITS4-KYO3 primers for ITS2.

.. image:: ../images/amp_lib_two.ITS3-KYO2_ITS4-KYO3.edit-graph.a75.svg
   :alt: Sequence edit-graph for amplicon library two using ITS3-KYO and ITS4-KYO3 primers for ITS2. Minimum abundance threshold 75.
