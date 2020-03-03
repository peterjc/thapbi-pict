.. _worked_example:

Endangered Species Mixes 16S etc
================================

This is the most complicated of the examples considered, where most of the
samples are "Experimental mixtures" of multiple plants and animals (plus two
traditional medicine mixtures where the exact content is unknown), which have
all be sequenced with about a dozen different primer pairs for multiple
metabarcoding markers including 16S, COI, cyt-b, matK, rbcL, trnL and ITS2.

This example pushes THAPBI PICT beyond its current design goals. It requires
setting up multiple databases, and running the tool multiple times - and then
any cross-barcode reporting or summation is left to the user. This includes
the built in sample assessment, where applying the same expected species
to all the barcodes is overly simplistic (e.g. the 16S markers should never
report any plant species).

.. toctree::
   :maxdepth: 1

   sample_data
   animals
   plants
