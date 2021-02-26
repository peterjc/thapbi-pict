.. _great_lakes:

Great Lakes Mock Community 16S
==============================

This example is based on:

    Klymus *et al.* (2017) Environmental DNA (eDNA) metabarcoding assays to
    detect invasive invertebrate species in the Great Lakes.
    https://doi.org/10.1371/journal.pone.0177643

Our main focus is 5 mock communities of 11 marine species in different ratios
(Table 2). The target amplicon copy number varies from trace level (14 reads)
to high copy number (9090 reads), making this an interesting example to
examine THAPBI PICT's minimum read abundance setting.

Two different sets primers were used targeting overlapping regions of the
mtDNA 16S RNA marker gene, named MOL16S and SPH16S, which were sequenced
separately (and not as in some of the other examples pooled together for the
Illumina sequencing). The current THAPBI PICT design requires running the
pipeline twice.

The full dataset includes aquarium and river environmental samples too, but
public sequence databases lack many of the sequences detected.

.. toctree::
   :maxdepth: 1

   sample_data
   abundance
   presence_absence
   edit_graph
