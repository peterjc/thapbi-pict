Sample Data for THAPBI PICT Worked Example
==========================================

This example is based on the following paper from earlier in the THAPBI
project, where the original analysis used the precursor pipeline ``metapy``:

* Riddell et al (2019) Metabarcoding reveals a high diversity of woody
  host-associated Phytophthora spp. in soils at public gardens and amenity
  woodlands in Britain. https://doi.org/10.7717/peerj.6931

The raw FASTQ files are too large to include with the THAPBI PICT source code,
so to follow the complete example you must download 244 ``*.fastq.gz`` files
separately (122 pairs, a little over 200MB in total).

There are multiple replicates from each of 14 sample sites, with FASTQ files
``Site_<N>_sample_<X>_R1.fastq.gz`` and ``Site_<N>_sample_<X>_R2.fastq.gz``
(plus controls with a different pattern), which as the first step of the
typical THAPBO PICT workflow (``thapbi_pict prepare-reads``) are transformed
in FASTA files named ``Site_<N>_sample_<X>.fasta`` (etc). We include these
FASTA files here as a compressed file, so after decompression they can be
used to follow the rest of a typical analysis:

.. code:: console

   $ tar -jxvf woody_hosts_fasta.tar.bz2

Note that four of the FASTA files are empty, ``Site_13_sample_7.fasta`` and
``Site_9_sample_4-3.fasta`` (nothing above the minimum threshold), and both
negative controls (good).
   
File ``site_metadata.tsv`` is a table of metadata (based on table S1 in the
paper), with one row for each of the 14 samples plus controls, with a cross
reference to the 122 sequenced FASTQ filename stems.

These files are used for the tutorial on https://thapbi-pict.readthedocs.io/
and also for automated testing during the tool development.
