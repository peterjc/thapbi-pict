Pooled animal and plant DNA barcodes
====================================

We have very briefly reviewed the output of each of the animal and plant
markers, noting some have no sequences at the THAPBI PICT default minimum
abundance threshold. Now we discuss the pooled results.

Sample report
-------------

Please open the ``summary/pooled.samples.onebp.xlsx`` sample report, zoomed
out you should have something like this:

.. image:: https://user-images.githubusercontent.com/63959/118682572-7c4b1580-b7f8-11eb-8973-15f1b50543f9.png
   :alt: Excel screenshot showing summary/pooled.samples.onebp.xlsx

The final column is the unknowns - and even at this zoom it is possible to see
a solid red region for the two traditional medicine samples (wide green
background bands).

Read report
-----------

To look at the unknown reads see ``summary/pooled.reads.onebp.xlsx``. Sorting
be the species prediction and zooming out should show something like this
where the top half of the rows are those sequences with a species prediction.
It is clear that the majority of the unknown sequences are from the two
traditional medicine samples (wide green bands):

.. image:: https://user-images.githubusercontent.com/63959/118682235-29715e00-b7f8-11eb-8dfb-bf18153a1ffa.png
   :alt: Excel screenshot showing pooled.reads.onebp.xlsx

Overall the replicates are reassuringly consistent - look at neighbouring
rows/columns within the colour bands in the two reports.

Pooled classifier assessment
----------------------------

The automated model assessment output in ``summary/pooled.assess.onebp.tsv``
is also worth review. Note this only looks at the experimental mixtures where
there is a ground truth (S1, S2, S4, S5, S6, S7, S9 and S10) - not the
traditional medicine samples where the true species content is unknown.

.. code:: console

    $ cut -f 1-5,9,11 summary/pooled.assess.onebp.tsv
    <SEE TABLE BELOW>

Working at the command line or using Excel should show the following:

====================== ==== === === ==== ==== ===========
#Species               TP   FP  FN  TN   F1   Ad-hoc-loss
====================== ==== === === ==== ==== ===========
OVERALL                1058 727 240 7699 0.69 0.478
Acipenser schrenckii   0    20  0   123  0.00 1.000
Aloe reynoldsii        0    114 0   29   0.00 1.000
Aloe variegata         110  0   25  8    0.90 0.185
Anguilla anguilla      3    0   3   137  0.67 0.500
Beta vulgaris          0    0   16  127  0.00 1.000
Bos taurus             139  2   0   2    0.99 0.014
Brassica juncea        0    127 0   16   0.00 1.000
Brassica napus         7    0   9   127  0.61 0.562
Brassica nigra         0    127 0   16   0.00 1.000
Brassica oleracea      128  6   0   9    0.98 0.045
Brassicaceae (misc)    0    70  0   73   0.00 1.000
Cactaceae (misc)       0    3   0   140  0.00 1.000
Carica papaya          16   0   0   127  1.00 0.000
Crocodylus niloticus   122  0   12  9    0.95 0.090
Cullen sp.             0    16  0   127  0.00 1.000
Cycas revoluta         3    0   3   137  0.67 0.500
Dendrobium sp.         131  0   3   9    0.99 0.022
Echinocactus sp.       6    0   0   137  1.00 0.000
Euphorbia sp.          3    0   3   137  0.67 0.500
Gallus gallus          6    1   0   136  0.92 0.143
Glycine max            16   0   0   127  1.00 0.000
Gossypium hirsutum     16   0   0   127  1.00 0.000
Homo sapiens           0    2   0   141  0.00 1.000
Huso dauricus          112  0   16  15   0.93 0.125
Lactuca altaica        0    66  0   77   0.00 1.000
Lactuca sativa         74   2   0   67   0.99 0.026
Lactuca serriola       0    66  0   77   0.00 1.000
Lactuca tatarica       0    39  0   104  0.00 1.000
Lactuca virosa         0    66  0   77   0.00 1.000
Meleagris gallopavo    16   0   0   127  1.00 0.000
Parapenaeopsis sp.     0    0   6   137  0.00 1.000
Pieris brassicae       6    0   0   137  1.00 0.000
Pleuronectes platessa  64   0   0   79   1.00 0.000
Solanum lycopersicum   16   0   0   127  1.00 0.000
Sus scrofa             64   0   0   79   1.00 0.000
Triticum aestivum      0    0   16  127  0.00 1.000
Zea mays               0    0   128 15   0.00 1.000
OTHER 31 SPECIES IN DB 0    0   0   4433 0.00 0.000
====================== ==== === === ==== ==== ===========

Most of the false positives (FP) are alternative genus level matches in
*Brassica* and *Lactuca* (as discussed in the paper). The two sequences we
recorded in the Mini-rbcL reference set as the family Brassicaceae are likely
also *Brassica*. The trnL-P6-loop marker had references for *Aloe reynoldsii*
but these matches are most likely from *Aloe variegata*.

A couple of the unique sequences are in the Mini-rbcL reference as the family
Cactaceae, and since they only appeared in the experimental mixes, these are
likely *Echinocactus* sp. or *Euphorbia* sp.

There are more interesting FP for *Acipenser schrenckii* (the authors found
this was accidentally included from the *Huso dauricus* caviar used), human
(*Homo sapiens*, presumed laboratory contamination), and finally *Lactuca
sativa*, cow (*Bos taurus*) and chicken (*Gallus gallus*) which the authors
traced to cross-contamination during sample preparation or DNA isolation.

Why do *Cullen* sp. show up in ``S1`` from the trnL P6 loop marker (as well
as ``S3`` which the authors found too, see their Table 8)?

If the sample database had been more inclusive there would have been many
more false positives. For example, the trnL-UAA sequence perfectly matching
AP007232.1 *Lactuca sativa* is also a perfect match for MK064549.1 *Luisia
teres*. Similarly, the Mini-rbcL sequence perfectly matching AP012989.1
*Brassica nigra* and MG872827.1 *Brassica juncea* also matches MN056359.2
*Raphanus sativus* (and more). This demonstrates the difficulties in curating
an appropriate marker database - and the content should depend in part on your
target samples.

Currently the provided references sequences (and thus classification databases
used) lack any markers for *Beta vulgaris*, *Parapenaeopsis* sp.,
*Triticum aestivum* or *Zea mays*. Most of these were present at only a few
percent dry weight, and are likely present below the default minimum abundance
threshold. This explains the false negatives.

Conclusion
----------

It appears that the THAPBI PICT default minimum abundance threshold of 100
reads is too stringent for detecting all the markers in a complex pool like
this. Including negative sequencing controls would help set an objective
lower bound.

There also appear to be marker sequences in these control samples which have
not yet been published, which would help by filling in gaps in the reference
set used for classification.

Also note we did not look at the multi-primer COI long marker, and perhaps the
default ``onebp`` classifier is not appropriate for the Mini-rbcL marker.
