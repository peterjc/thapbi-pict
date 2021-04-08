Universal plant DNA barcodes and mini-barcodes
==============================================

As in the animal primers, for rbcL the paper used two targets, a long barcode
and a shorter mini-barcode. The same names have been used in the ``run.sh``
script provided, the output of which is referred to below.

matK
----

The paper described two sets of primers for matK, although only one was used
for the MiSeq sequencing. This gave no sequences at the default abundance
threshold, dropping to 50 showed three uniques sequences in three files, and
even dropping to 10 only gave results from ``EM_2``, ``EM_14`` and ``S8``.

NCBI BLAST of these sequence gave no perfect matches, but suggested
*Sanguisorba* sp. was present, noted in the original paper for ``S8`` which
is one of the traditional medicine samples.

rbcL - long target
------------------

Using our default abundance threshold and the author's minimum length of 140bp,
we got no sequences at all. Allowing a minimum length of 100 (our default)
gave the following sequence and a one SNP variant, all from ``S3``::

    >3ec67342f519461a0ad40fef436b1b1d
    GACTGCGGGGTTCAAAGCTGGTGTTAAAGATTATAGATTGACGTATTATACTCCTGAATTGGGGTTATCCGCTAAGAATT
    ACGGTAGAGCAGTTTATGAATGTCTT

The best NCBI BLAST matches are *Astragalus*, but with a break point. The
authors of the original paper report finding *Astragalus danicus* in ``S3``.

Mini-rbcL - short target
------------------------

This was by far and above the most diverse in terms of unique sequences
recovered:

.. code:: console

    $ grep -c ">" summary/Mini-rbcL.all_reads.fasta
    279

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

Not very diverse, only eight unique sequences recovered:

.. code:: console

    $ grep -c ">" summary/trnL-UAA.all_reads.fasta
    8

We see lots of *Brassica*, the difficulties with *Brassica oleracea* vs
*Brassica napus* (and the genus in general) are discussed in the paper too.

trnL-P6-loop
------------

This gave no sequences at the default abundance threshold, and disabling
the abundance threshold no unique trimmed sequence was seen more than once.
This is strange, given the authors report finding *Lactuca sativa* and
*Cycas revoluta* from this primer. It is however easily explained, quoting
the paper:

   We implemented a minimum DNA barcode length of 200 nt, except for DNA
   barcodes with a basic length shorter than 200 nt, in which case the
   minimum expected DNA barcode length is set to 100 nt for ITS2, 140 nt
   for mini-rbcL, and 10 nt for the trnL (P6 loop) marker.

We should have changed the THAPBI PICT minimum length from 100 (our default)
to 10 as well.

ITS2
----

Quite diverse, with over fifty unique sequences recovered:

.. code:: console

    $ grep -c ">" summary/ITS2.all_reads.fasta
    59

Finds all the *Brassica* and *Echinocactus* sp., most of the *Euphorbia* sp.

We do see unexpected matches to *Lactuca* sp. where *Lactuca sativa* was in
the experimental mixture. The dominant sequence present is just one base pair
away from a published sequence from that species (KM210323.1), but perfectly
matches published sequences from *Lactuca altaica*, *L. serriola* and
*L. virosa* - and that is what was in the sample database. If you open the
associated edit-graph file (``ITS2.edit-graph.onebp.xgmml``) in Cytoscape,
you can see this quite clearly.
