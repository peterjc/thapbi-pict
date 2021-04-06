.. _sample_data:

Marker data
===========

Either clone the THAPBI PICT source code repository, or decompress the
latest source code release (``.tar.gz`` file). You should find it contains
a directory ``examples/woody_hosts/`` which is for this example.

Shell scripts ``setup.sh`` and ``run.sh`` should reproduce the analysis
discussed.

The documentation goes through running each step of the analysis gradually,
before finally calling pipeline command to do it all together. We provide
script ``run.sh`` to do the final run-though automatically (first without
any metadata, then again with it), but encourage you to follow along the
individual steps first.

If you skip the raw FASTQ files, you won't be able to run the entire pipeline
or the prepare-reads steps.

FASTQ data
----------

The raw data is from two Illumina MiSeq runs, a whole 96-well plate from 2016,
and about half the samples from a second 96-well plate sequenced in 2017
(where the rest of the plate was samples from a separate ITS1 study). There
are multiple replicates from each of 14 sample sites, plus controls.
The raw FASTQ files are too large to include with the THAPBI PICT source code.

Script ``setup.sh`` will download the raw FASTQ files for Riddell *et al.*
(2019) from https://doi.org/10.5281/zenodo.3342957

It will download 244 raw FASTQ files (122 pairs), about 215MB on disk

The first step of a typical THAPBI PICT workflow is to transform the paired
FASTQ files into much smaller FASTA files. We provide those FASTA files
compressed with the THAPBI PICT source code, so if you skip the FASTQ files,
you can still follow the rest of a typical analysis.

If you skip the raw data, instead you must decompress the pre-prepared 122
FASTA files into your ``intermediate/`` subdirectory:

.. code:: console

    $ tar -jxvkf intermediate.tar.bz2
    ...
    $ ls -1 intermediate/*.fasta | wc -l
    122

Note that four of these FASTA files are empty, ``Site_13_sample_7.fasta`` and
``Site_9_sample_4-3.fasta`` (nothing above the minimum threshold), and both
negative controls (good).

We also provide :ref:`metadata for the samples <metadata>` for use in the
reports.

Amplicon primers & reference sequences
--------------------------------------

The ITS1 primers used here match the THAPBI PICT defaults, so the default
database can also be used.

Metadata
--------

The provided file ``metadata.tsv`` is an expanded version of Supplementary
Table 1 from the original paper, adding a column for the Illumina MiSeq sample
names, rows for the controls, and a comment block at the start.

The 16 columns are as follows, where 4 to 15 are in pairs for tree/shrub broad
taxonomic grouping and health status (H, healthy; D, symptoms/stump/dead):

1. Anonymised site number (with leading zero, "01" to "14"), or control name
2. Approximate altitude at centre
3. Underlying soil type
4. Healthy Cupressaceae
5. Diseased Cupressaceae
6. Healthy other conifers
7. Diseased other conifers
8. Healthy Ericaceae
9. Diseased Ericaceae
10. Healthy Fagaceae or Nothofagaceae
11. Diseased Fagaceae or Nothofagaceae
12. Healthy other angiosperms
13. Diseased other angiosperms
14. Healthy other
15. Diseased other
16. MiSeq Sample(s) (semi-colon separated list)

Lines 1 to 19 are human readable header text, line 20 is the column headers.
Lines 21 onwards are data for 14 field sites and 3 controls.

When calling THAPBI PICT, the meta data commands are given as follows:

.. code:: console

    $ thapbi_pict ... -t metadata.tsv -x 16 -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -f 20

These settings are described in detail later. This example is important in
that column 16 contains multiple entries where a site had multiple sequenced
samples (replicates). It is unusual in having comment lines before the column
header line which must be specified.

Other files
-----------

Compressed archive file ``intermediate.tar.bz2`` contains the results of
the first step in the analysis, in case you have problems downloading the
raw FASTQ files (see below).

Subdirectory ``expected/`` contains four plain text tab-separated files,
describing the expected species in some mock community positive controls:

* ``DNA15MIX.known.tsv``
* ``DNA10MIX_bycopynumber.known.tsv``
* ``DNA10MIX_diluted25x.known.tsv``
* ``DNA10MIX_undiluted.known.tsv``
