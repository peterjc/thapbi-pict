.. _sample_data:

Introduction
============

Sample data
-----------

This example is based on the following paper from earlier in the THAPBI
Phyto-Threats project, where the original analysis used the precursor pipeline
``metapy``:

* Riddell et al (2019) Metabarcoding reveals a high diversity of woody
  host-associated *Phytophthora* spp. in soils at public gardens and amenity
  woodlands in Britain. https://doi.org/10.7717/peerj.6931

Importantly, they used the same PCR primers, and therefore analysis with the
tool's default settings including the provided database is appropriate.

The raw data is from two Illumina MiSeq runs, a whole 96-well plate from 2016,
and about half the samples from a second 96-well plate sequenced in 2017
(where the rest of the plate was samples from a separate ITS1 study). There
are multiple replicates from each of 14 sample sites, plus controls.

The raw FASTQ files are too large to include with the THAPBI PICT source code,
so to follow the complete example you must download 244 ``*.fastq.gz`` files
separately (122 pairs, 215MB compressed) from `DOI 10.5281/zenodo.3342957
<https://doi.org/10.5281/zenodo.3342957>`_.

The first step of a typical THAPBI PICT workflow is to transform the paired
FASTQ files into much smaller FASTA files. We provide those FASTA files
compressed with the THAPBI PICT source code, so if you skip the FASTQ files,
you can still follow the rest of a typical analysis.

We also provide :ref:`metadata for the samples <metadata>` for use in the
reports.

Setup
-----

We assume you have your command line terminal open in a new empty folder
dedicated to this analysis. Start by making three sub-folders as follows:

.. code:: console

   $ mkdir raw_data/ intermediate/ summary/ expected/

We will need file ``metadata.tsv`` (included with the THAPBI PICT source
code as ``tests/woody_hosts/metadata.tsv``), containing :ref:`metadata
about the sites <metadata>`. This file can be downloaded:

.. code:: console

    $ wget https://github.com/peterjc/thapbi-pict/raw/master/tests/woody_hosts/metadata.tsv

We will also need four files ``expected/*.known.tsv`` which describe
the species expected in the positive control sequences. These are included
with the source under ``tests/woody_hosts/``, or can be downloaded:

.. code:: console

    $ cd expected/
    $ wget https://github.com/peterjc/thapbi-pict/raw/master/examples/woody_hosts/expected/DNA15MIX.known.tsv
    $ wget https://github.com/peterjc/thapbi-pict/raw/master/examples/woody_hosts/expected/DNA10MIX_bycopynumber.known.tsv
    $ wget https://github.com/peterjc/thapbi-pict/raw/master/examples/woody_hosts/expected/DNA10MIX_diluted25x.known.tsv
    $ wget https://github.com/peterjc/thapbi-pict/raw/master/examples/woody_hosts/expected/DNA10MIX_undiluted.known.tsv
    $ cd ..

The FASTQ files are only needed for the very first step of the worked example,
or for running the combined pipeline command. The download is about 215MB:

.. code:: console

    $ wget https://zenodo.org/record/3342957/files/raw_data.tar.gz
    $ tar -zxvf raw_data.tar.gz
    $ ls -1 raw_data/*.fastq.gz | wc -l
    244

Ideally confirm the MD5 checksum is ``3435aa5b567897ab23e6aac46fb822a9`` by
running ``md5sum raw_data.tar.gz`` (on Linux) or ``md5 raw_data.tar.gz`` (on
macOS). This is just to double check there were no problems with the download.

If you skip the raw data, instead you must download the pre-prepared 122 FASTA
files into your intermediate data sub-folder as ``intermediate/*.fasta``.
These are provided as a small compressed file included in the THAPBI PICT
source code ``tests/woody_hosts/intermediate.tar.bz2``, or can easily be
downloaded:

.. code:: console

   $ wget https://github.com/peterjc/thapbi-pict/raw/master/tests/woody_hosts/intermediate.tar.bz2
   $ tar -jxvf intermediate.tar.bz2
   $ ls -1 intermediate/*.fasta | wc -l
   122

Note that four of the FASTA files are empty, ``Site_13_sample_7.fasta`` and
``Site_9_sample_4-3.fasta`` (nothing above the minimum threshold), and both
negative controls (good).
