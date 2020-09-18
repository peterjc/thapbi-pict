.. _worked_examples:

Worked Examples
===============

While THAPBI PICT stands for *Phytophthora* ITS1 Classification Tool, with
appropriate primer settings and a custom database, it can be applied to other
organisms and/or barcode marker sequences.

These worked examples use public datasets from published papers, with various
markers covering oomycetes, fungi, animals and plants. The main criteria has
been mock communities with known species composition.

.. toctree::
   :hidden:
   :maxdepth: 1

   woody_hosts/index
   recycled_water/index
   fungal_mock/index
   microalgal_mock/index
   fecal_sequel/index
   endangered_species/index

* :ref:`woody_hosts` - Simple ITS1 example using the default primers and
  database. Based on a paper from earlier in the THAPBI Phyto-Threats project:

    Riddell *et al.* (2019) Metabarcoding reveals a high diversity of woody
    host-associated *Phytophthora* spp. in soils at public gardens and
    amenity woodlands in Britain. https://doi.org/10.7717/peerj.6931

* :ref:`recycled_water` - ITS1 example where defaults can be used, but
  ideally requires different primers and a custom database. Based on:

    Redekar *et al.* (2019) Diversity of *Phytophthora*, *Pythium*, and
    *Phytopythium* species in recycled irrigation water in a container
    nursery. https://doi.org/10.1094/PBIOMES-10-18-0043-R

* :ref:`fungal_mock` - ITS1 and ITS2 example requiring multiple primers and
  databases, based on:

    Bakker (2018) A fungal mock community control for amplicon sequencing
    experiments. https://doi.org/10.1111/1755-0998.12760

* :ref:`microalgal_mock` - 18S rRNA example requiring multiple primers and
  databases, based on:

    Bradley *et al.* (2016) Design and Evaluation of Illumina
    MiSeq-Compatible, 18S rRNA Gene-Specific Primers for Improved
    Characterization of Mixed Phototrophic Communities.
    https://doi.org/10.1128/AEM.01630-16

* :ref:`fecal_sequel` - COI example in bats, showing importance of the
  database content with the default classifier. Based on:

    Walker *et al.* (2019) A fecal sequel: Testing the limits of a genetic
    assay for bat species identification.
    https://doi.org/10.1371/journal.pone.0224969

* :ref:`endangered_species` - A dozen markers in animals and plants, pushing
  THAPBI PICT beyond its original design goals by pooling markers for
  sample level analysis. Based on:

    Arulandhu *et al.* (2017) Development and validation of a multi-locus DNA
    metabarcoding method to identify endangered species in complex samples.
    https://doi.org/10.1093/gigascience/gix080

For each worked example there is a different sub-folder in the THAPBI PICT
source code under ``examples/`` containing ``setup.sh`` to download the public
data, ``run.sh`` to execute the main analysis discussed, and assorted other
files like ``metadata.tsv``.
