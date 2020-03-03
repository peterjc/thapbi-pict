Universal plant DNA barcodes and mini-barcodes
==============================================

As in the animal primers, for rbcL the paper used two targets, a long barcode
and a shorter mini-barcode. The same names have been used in the ``run.sh``
script provided, the output of which is referred to below.

matK
----

The paper described two sets of primers for matK, although only one was used
for the MiSeq sequencing. This gave no sequences at the default abundance
threshold.

rbcL - long target
------------------

This gave no sequences at the default abundance threshold.

Mini-rbcL - short target
------------------------

This was by far and above the most diverse in terms of unique sequences
recovered:

.. code:: console

    $ grep -c ">" summary/Mini-rbcL.all.fasta
    215

We see expected plant species like *Lactuca sativa*, *Brassica oleracea*,
*Aloe variegata* and *Dendrobium sp.* - exactly how they are classified
depends critically on how the database is built.

The traditional medicine samples have multiple unknown sequences likely of
plant origin.

The edit-graph is the most complicated of those in this dataset - not
simply in terms of the number of nodes. This marker needs more careful
review before using THAPBI PICT's default ``onebp`` classifier.

trnl-UAA
--------

Not very diverse, only nine unique sequences recovered:

.. code:: console

    $ grep -c ">" summary/trnl-UAA.all.fasta
    9

We see lots of *Brassica*, the difficulties with *Brassica oleracea* vs
*Brassica napus* (and the genus in general) are discussed in the paper too.

trnL-P6-loop
------------

This gave no sequences at the default abundance threshold.

ITS2
----

Quite diverse, with over fifty unique sequences recovered:

.. code:: console

    $ grep -c ">" summary/ITS2.all.fasta
    9

Finds all the *Brassica* and *Echinocactus* sp., most of the *Euphorbia* sp.
Unexpected matches to *Lactuca* sp.
