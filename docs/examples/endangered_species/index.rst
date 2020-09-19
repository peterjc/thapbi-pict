.. _endangered_species:

Endangered Species Mixes 16S etc
================================

This is the most complicated of the examples considered, where most of the
samples are "Experimental mixtures" of multiple plants and animals (plus two
traditional medicine mixtures where the exact content is unknown), which have
all been sequenced with about a dozen different primer pairs for multiple
metabarcoding markers including 16S, COI, cyt-b, matK, rbcL, trnL and ITS2:

    Arulandhu *et al.* (2017) Development and validation of a multi-locus DNA
    metabarcoding method to identify endangered species in complex samples.
    https://doi.org/10.1093/gigascience/gix080

This example pushes THAPBI PICT beyond its original design goals. It requires
setting up multiple databases (all of which ought to be properly curated),
and running the tool multiple times (where potentially different thresholds
might be needed). In order to provide cross-barcode reporting, the ``run.sh``
script ends by pooling all the per-marker intermediates by sample, and
generating combined reports.

.. toctree::
   :maxdepth: 1

   sample_data
   animals
   plants
   pooled
