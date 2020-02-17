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

Provided files
--------------

Either clone the THAPBI PICT source code repository, or decompress the
latest source code release (``.tar.gz`` file). You should find it contains
a directory ``examples/woody_hosts/`` which is for this example.

File ``metadata.tsv`` contains metadata about the sampling sites (see below).

Compressed archive file ``intermediate.tar.bz2`` contains the results of
the first step in the analysis, in case you have problems downloading the
raw FASTQ files (see below).

Subdirectory ``expected/`` contains four plain text tab-separated files,
describing the expected species in some mock community positive controls:

* ``DNA15MIX.known.tsv``
* ``DNA10MIX_bycopynumber.known.tsv``
* ``DNA10MIX_diluted25x.known.tsv``
* ``DNA10MIX_undiluted.known.tsv``

Shell scripts ``setup.py`` and ``run.sh`` should reproduce the analysis
discussed in the THAPBI PICT documentation.

Setup
-----

We assume you have aquired the THAPBI PICT source code, and have your command
line terminal open in the ``examples/woody_hosts/`` folder. First we run the
``setup.py`` script:

.. code:: console

   $ ./setup.py

This will download the raw gzip compressed FASTQ files from Zenodo (244 files,
122 pairs, a little over 200MB in total) into the ``raw_data/`` subdirectory.

If you skip the raw data, instead you must decompress the pre-prepared 122
FASTA files into your ``intermediate/`` subdirectory:
.. code:: console

   $ tar -jxvf intermediate.tar.bz2
   $ ls -1 intermediate/*.fasta | wc -l
   122

Note that four of these FASTA files are empty, ``Site_13_sample_7.fasta`` and
``Site_9_sample_4-3.fasta`` (nothing above the minimum threshold), and both
negative controls (good).

Running the pipeline
--------------------

The documentation goes through running each step of the analysis gradually,
before finally calling pipeline command to do it all together. We provide
script ``run.sh`` to do the final run-though automatically (first without
any metadata, then again with it), but encourage you to follow along the
individual steps first.

If you skip the raw FASTQ files, you won't be able to run the entire pipeline
or the prepare-reads steps.
