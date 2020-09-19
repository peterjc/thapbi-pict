.. _recycled_water:

Environmental Oomycetes ITS1
============================

The first :ref:`worked example <woody_hosts>` looked at *Phytophthora* ITS1
data from woody-host trees, using the same PCR primers as the THAPBI PICT
defaults, and the default database of *Phytophthora* ITS1 data provided.

Here we re-analyse a published dataset from a different group, doing Illumina
MiSeq ITS1 amplicon sequencing of irrigation water samples from Oregon:

    Redekar *et al.* (2019) Diversity of *Phytophthora*, *Pythium*, and
    *Phytopythium* species in recycled irrigation water in a container nursery.
    https://doi.org/10.1094/PBIOMES-10-18-0043-R

Different PCR primers were used to cover *Pythium* and *Phytopythium* as well
as *Phytophthora*. This requires a new database of markers.

.. toctree::
   :maxdepth: 1

   sample_data
   defaults
   primers
   build_db
   examine_db
   pipeline
