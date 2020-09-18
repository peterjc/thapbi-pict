.. _recycled_water:

Environmental Oomycetes ITS1
============================

The first :ref:`worked example <woody_hosts>` applied the default analysis
pipeline to some real *Phytophthora* ITS1 data from woody-host trees, using the
same PCR primers as the THAPBI PICT defaults, and with the default database
of *Phytophthora* ITS1 data the tool provides.

Here we re-analyse a published dataset from a different group, doing Illumina
MiSeq ITS1 amplicon sequencing of irrigation water samples from Oregon. This
covered *Pythium* and *Phytopythium* as well as *Phytophthora*, and therefore
used different PCR primers. That requires a different database of markers.

.. toctree::
   :maxdepth: 1

   sample_data
   defaults
   primers
   build_db
   examine_db
   pipeline
