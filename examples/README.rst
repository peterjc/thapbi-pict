This folder contains multiple examples for use with THAPBI PICT, each with
their own ``README.rst`` file for context (part of a larger set of worked
example documentation).

You can use the ``setup_all.sh`` to download and setup all the raw data,
and likewise ``run_all.sh`` to general all the intermediate files and the
output reports.

1. ``woody_hosts/`` - simple ITS1 example using default primers and database,
   based on a paper from earlier in the THAPBI Phyto-Threats project:

   * Riddell *et al.* (2019) Metabarcoding reveals a high diversity of woody
     host-associated *Phytophthora* spp. in soils at public gardens and
     amenity woodlands in Britain. https://doi.org/10.7717/peerj.6931

2. ``recycled_water/`` - ITS1 example where defaults can be used but idealy
   requires different primers and a custom database, based on:

   * Redekar *et al.* (2019) Diversity of *Phytophthora*, *Pythium*, and
     *Phytopythium* species in recycled irrigation water in a container
     nursery. https://doi.org/10.1094/PBIOMES-10-18-0043-R

3. ``fungal_mock/`` - ITS1 and ITS2 example requiring multiple primers and
   databases, based on:

   * Bakker (2018) A fungal mock community control for amplicon sequencing
     experiments. https://doi.org/10.1111/1755-0998.12760

4. ``microalgal_mock/`` - 18S rRNA example requiring multiple primers and
   databases, based on:

   * Bradley *et al.* (2016) Design and Evaluation of Illumina MiSeq-Compatible,
     18S rRNA Gene-Specific Primers for Improved Characterization of Mixed
     Phototrophic Communities. https://doi.org/10.1128/AEM.01630-16

5. ``fecal_sequel/`` - COI example in bats, showing importance of the database
   content with the default classifier, based on:

   * Walker *et al.* (2019) A fecal sequel: Testing the limits of a genetic
     assay for bat species identification.
     https://doi.org/10.1371/journal.pone.0224969

6. ``endangered_species/`` - A dozen markers in animals and plants, pushing
   THAPBI PICT beyond its original design goals by pooling markers for
   sample level analysis. Based on:

   * Arulandhu *et al.* (2017) Development and validation of a multi-locus DNA
     metabarcoding method to identify endangered species in complex samples.
     https://doi.org/10.1093/gigascience/gix080

See https://thapbi-pict.readthedocs.io/en/latest/examples/ for discussion of
these examples.
