Sample Data
-----------

This example is based on the following paper from earlier in the THAPBI
Phyto-Threats project, where the original analysis used the precursor pipeline
``metapy``:

* Riddell et al (2019) Metabarcoding reveals a high diversity of woody
  host-associated Phytophthora spp. in soils at public gardens and amenity
  woodlands in Britain. https://doi.org/10.7717/peerj.6931

The raw data is from two Illumina MiSeq runs, a whole 96-well plate from 2016,
and about half the samples from a second 96-well plate sequenced in 2017
(where the rest of the plate was samples from a separate ITS1 study). There
are multiple replicates from each of 14 sample sites, plus controls.

The raw FASTQ files are too large to include with the THAPBI PICT source code,
so to follow the complete example you must download 244 ``*.fastq.gz`` files
separately (122 pairs, a little over 200MB in total).

The first step of a typical THAPBI PICT workflow is to transform the paired
FASTQ files into much smaller FASTA files. We provide those FASTA files
compressed with the THAPBI PICT source code, so they can be used to follow the
rest of a typical analysis. We also provide metadata for the samples for use
in the reports.

Setup
-----

We assume you have your command line terminal open in a new empty folder
dedicated to this analysis. Start by making three sub-folders as follows:

.. code:: console

   $ mkdir raw_data/ intermediate/ summary/

We will need file ``site_metadata.tsv`` (included with the THAPBI PICT source
code as ``tests/woody_hosts/site_metadata.tsv``) which can be downloaded:

.. code:: console

    $ wget https://github.com/peterjc/thapbi-pict/raw/master/tests/woody_hosts/site_metadata.tsv

The FASTQ files are only needed for the very first step of the worked example.
If you have downloaded the 244 paired FASTQ files, put them in the raw data
sub-folder as ``raw_data/*.fastq.gz``.

If you don't have the FASTQ files, you need to get the pre-prepared 122 FASTA
files into your intermediate data sub-folder as ``intermediate/*.fasta``.
These are provided as a small compressed file included in the THAPBI PICT
source code ``tests/woody_hosts/woody_hosts_fasta.tar.bz2``, or can easily be
downloaded:

.. code:: console

   $ wget https://github.com/peterjc/thapbi-pict/raw/master/tests/woody_hosts/woody_hosts_fasta.tar.bz2
   $ tar -jxvf woody_hosts_fasta.tar.bz2 -C intermediate/

Note that four of the FASTA files are empty, ``Site_13_sample_7.fasta`` and
``Site_9_sample_4-3.fasta`` (nothing above the minimum threshold), and both
negative controls (good).
