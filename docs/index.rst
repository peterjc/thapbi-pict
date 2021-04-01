.. THAPBI PICT documentation master file, created by
   sphinx-quickstart on Thu Jun 27 12:25:27 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

THAPBI *Phytophthora* ITS1 Classifier Tool (PICT)
=================================================

THAPBI PICT is a sequence based diagnostic/profiling tool from the UK funded
Tree Health and Plant Biosecurity Initiative (THAPBI) `Phyto-Threats project
<https://www.forestresearch.gov.uk/research/global-threats-from-phytophthora-spp/>`_,
focused on identifying *Phytophthora* species present in Illumina sequenced
environmental samples.

*Phytophthora* (from Greek meaning plant-destroyer) species are
economically important plant pathogens, important in both agriculture
and forestry. ITS1 is short for Internal Transcribed Spacer one, which
is a region of eukaryotes genomes between the 18S and 5.8S rRNA genes.
This is commonly used for molecular barcoding, where sequencing this
short region can identify species.

With appropriate primer settings and a custom database of full length markers,
THAPBI PICT can be applied to other organisms and/or barcode marker sequences
- not just *Phytophthora* ITS1. It requires overlapping paired-end Illumina
reads which can be merged to cover the *full* amplicon marker. Longer markers
or fragmented amplicons are not supported. Internally it works by tracking
amplicon sequence variants (ASVs), using MD5 checksums as identifiers.

The worked examples include oomycetes, fungi, microalgae, and bats, and cover
markers in ITS1, ITS2, 18S rRNA and COI and more. The main criteria has been
mock communities with known species composition.

The THAPBI Phyto-Threats project is supported by a grant funded jointly by the
Biotechnology and Biological Sciences Research Council (`BBSRC
<https://bbsrc.ukri.org/>`_), the Department for Environment, Food and Rural
affairs (`DEFRA <https://www.gov.uk/government/organisations/department-for-environment-food-rural-affairs>`_),
the Economic and Social Research Council (`ESRC <https://esrc.ukri.org>`_),
the `Forestry Commission <https://www.gov.uk/government/organisations/forestry-commission>`_,
the Natural Environment Research Council (`NERC <https://nerc.ukri.org>`_)
and the `Scottish Government <https://www.gov.scot/>`_, under the Tree
Health and Plant Biosecurity Initiative (THAPBI).

Key links:

* Documentation on Read The Docs: https://thapbi-pict.readthedocs.io/
* Source code repository on GitHub: https://github.com/peterjc/thapbi-pict/
* Software released on PyPI: https://pypi.org/project/thapbi-pict/
* Zenodo DOI for software: https://doi.org/10.5281/zenodo.4529395

.. toctree::
   :maxdepth: 1
   :caption: Documentation contents:

   introduction
   installation
   quick_start
   examples/index
   database
   negatives
   commands
   api/index
   release_history
   development_notes

This documentation was generated for THAPBI PICT version |version|.
