.. _custom_database_sample_data:

Marker data
===========

Either clone the THAPBI PICT source code repository, or decompress the
latest source code release (``.tar.gz`` file). You should find it contains
a directory ``examples/oomycete_rps10/`` which is for this example.

Shell scripts ``setup.sh`` and ``run.sh`` should reproduce the analysis
discussed.

The documentation goes through running each step of the analysis gradually,
including building a combined database, before calling pipeline command.
We provide script ``run.sh`` to do the final run-though automatically, but
encourage you to follow along the individual steps first.

FASTQ data
----------

File ``PRJNA699663.tsv`` was download from the ENA and includes the FASTQ
checksums, URLs, and sample metadata. With a little scripting to extract the
relevant :ref:`sample metadata <metadata>` for use with THAPBI PICT this was
reformatted as ``metadata.tsv`` (see below).

Script ``setup.sh`` will download the raw FASTQ files for Foster *et al.*
(2021) from https://www.ebi.ac.uk/ena/data/view/PRJNA699663 - you could also
use https://www.ncbi.nlm.nih.gov/bioproject/PRJNA699663/ to get this.

It will download 64 raw FASTQ files (32 pairs), taking about 2.3GB on disk

If you have the ``md5sum`` tool installed (standard on Linux; we suggest
``conda install coreutils`` to install this on macOS), verify the FASTQ files
downloaded correctly:

.. code:: console

    $ cd raw_data/
    $ md5sum -c MD5SUM.txt
    ...
    $ cd ..

There is no need to decompress the files.

Amplicon primers & reference sequences
--------------------------------------

This example looks at samples amplified with the same ITS1 primers as the
THAPBI PICT default database, but also rps10 primers. The rps10 assay is a
multiplex PCR reaction comprising two rps10 forward primers and seven rps10
reverse primers that differ slightly in sequence but anneal to the same
position. We follow the authors in using IUPAC codes to approximate the
rps10 marker primers, with `GTTGGTTAGAGYARAAGACT` for the left and
`ATRYYTAGAAAGAYTYGAACT` for the right (reverse complement
`AGTTCRARTCTTTCTARRYAT`).

In order to classify the rps10 sequences, we need to build a THABPI PICT
database of full-length primer-trimmed references. Happily we can use the
authors' own reference sequences from OomyceteDB
https://oomycetedb.cgrb.oregonstate.edu
(after reformatting the FASTA headers to a pattern THAPBI PICT recognises).

Metadata
--------

The provided file ``metadata.tsv`` has seven columns:

1. run_accession, eg "SRR13658667"
2. tax_id, eg "410658"
3. scientific_name, eg "soil metagenome"
4. library_name, eg "C5"
5. experiment_title, eg "Illumina MiSeq sequencing; ITS1 sequences of agricultural soil"
6. Marker, "ITS1" or "rps10"
7. Source, eg  "Agricultural soil"

When calling THAPBI PICT, the meta data commands are given as follows:

.. code:: console

    $ thapbi_pict ... -t metadata.tsv -x 1 -c 3,7,4,6

Argument ``-t metadata.tsv`` says to use this file for the metadata.

The ``-x 1`` argument indicates the filename stem can be found in column 1,
Accession.

Argument ``-c 3,7,4,6`` says which columns to display and sort by (do
not include the indexed column again).

We have not given a ``-g`` argument to assign colour bands in the Excel
reports, so it will default to the first column in ``-c``, "scientific name"
meaning we get four coloured bands for "aquatic metagenome", "plant
metagenome", "soil metagenome", and "synthetic metagenome".

Other files
-----------

The setup script will create symlinks using the sample names under sub-folder
``expected/`` pointing at the relevant mock-community known file. This is for
automatically assessing the classifier performance.

Sub-folders under ``intermediate/`` are used for intermediate files, a folder
for each primer-pair. Subfolder ``summary/`` is used for the generated reports.
