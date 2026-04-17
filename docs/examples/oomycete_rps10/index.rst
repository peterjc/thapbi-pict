Oomycete ITS1 & rps10
=====================

The first :ref:`worked example <woody_hosts>` looked at *Phytophthora* ITS1
data from woody-host trees, using the same PCR primers as the THAPBI PICT
defaults, and the default database of *Phytophthora* ITS1 data provided.

Here we re-analyse a published dataset from a different group, doing Illumina
MiSeq ITS1 amplicon sequencing as a base-line to evaluate a mitochondrial 40S
ribosomal protein S10 (rps10) gene based marker as an alternative barcode for
the wider group of oomycetes:

    Foster *et al.* (2021) A New Oomycete Metabarcoding Method Using the rps10
    Gene. https://doi.org/10.1094/PBIOMES-02-22-0009-R

The example combines our default *Phytophthora* focused ITS1 database with
their broader rps10 reference sequences to process both markers together.

.. toctree::
   :maxdepth: 1

   sample_data
   build_db
   pipeline
