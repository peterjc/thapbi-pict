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
   drained_ponds/index
   fungal_mock/index
   great_lakes/index
   fecal_sequel/index
   synthetic_mycobiome/index
   soil_nematodes/index
   pest_insects/index
   endangered_species/index

* :ref:`woody_hosts` - A simple example using the default primers and
  database. Based on a paper from earlier in the THAPBI Phyto-Threats project:

    Riddell *et al.* (2019) Metabarcoding reveals a high diversity of woody
    host-associated *Phytophthora* spp. in soils at public gardens and
    amenity woodlands in Britain. https://doi.org/10.7717/peerj.6931

* :ref:`recycled_water` - An example where the defaults can be used, but
  ideally requires a different primer pair and a custom database. Based on:

    Redekar *et al.* (2019) Diversity of *Phytophthora*, *Pythium*, and
    *Phytopythium* species in recycled irrigation water in a container
    nursery. https://doi.org/10.1094/PBIOMES-10-18-0043-R

* :ref:`drained_ponds` - An example with a single marker and custom database.
  Based on:

    Muri *et al.* (2020) Read counts from environmental DNA (eDNA)
    metabarcoding reflect fish abundance and biomass in drained ponds.
    https://doi.org/10.3897/mbmg.4.56959

* :ref:`fungal_mock` - An example with multiple markers (including two
  sequenced together) requiring separate primers settings and databases, based
  on:

    Bakker (2018) A fungal mock community control for amplicon sequencing
    experiments. https://doi.org/10.1111/1755-0998.12760

* :ref:`great_lakes` - An example with two mitochondrial markers (sequenced
  separately), with mock communities, where we focus on the minimum abundance
  threshold. Based on:

    Klymus *et al.* (2017) Environmental DNA (eDNA) metabarcoding assays to
    detect invasive invertebrate species in the Great Lakes.
    https://doi.org/10.1371/journal.pone.0177643

* :ref:`fecal_sequel` - A single marker example in bats, showing importance of
  the database content with the default classifier. Based on:

    Walker *et al.* (2019) A fecal sequel: Testing the limits of a genetic
    assay for bat species identification.
    https://doi.org/10.1371/journal.pone.0224969

* :ref:`synthetic_mycobiome` - A single marker example in fungi, with mock
  biological communities *and* synthetic control sequences. Based on:

    Palmer *et al.* (2018) Non-biological synthetic spike-in controls and the
    AMPtk software pipeline improve mycobiome data.
    https://doi.org/10.7717/peerj.4925

* :ref:`soil_nematodes` - Four markers (sequenced separately) in a soil nematode
  mock community. Based on:

    Ahmed *et al.* (2019) Metabarcoding of soil nematodes: the importance of
    taxonomic coverage and availability of reference sequences in choosing
    suitable marker(s)
    https://doi.org/10.3897/mbmg.3.36408

* :ref:`pest_insects` - Three markers (sequenced together) in insect mock
  communities. Based on:

    Batovska *et al.* (2021) Developing a non-destructive metabarcoding
    protocol for detection of pest insects in bulk trap catches
    https://doi.org/10.1038/s41598-021-85855-6

* :ref:`endangered_species` - A dozen markers in animals and plants (sequenced
  together). Based on:

    Arulandhu *et al.* (2017) Development and validation of a multi-locus DNA
    metabarcoding method to identify endangered species in complex samples.
    https://doi.org/10.1093/gigascience/gix080

For each worked example there is a different sub-folder in the THAPBI PICT
source code under ``examples/`` containing at least ``setup.sh`` to do one-off
setup like downloading the public data, and ``run.sh`` to execute the main
analysis discussed. There will usually be assorted other files like reference
sequences, or ``metadata.tsv``.

Running the examples will create or use subdirectories ``raw_data/`` for the
downloaded FASTQ files, ``intermediate/`` for per-sample working files, and
``summary/`` for the final output reports. Where the example includes positive
controls like mock communities, the expected species content is recorded under
``expected/`` in per-sample files.
