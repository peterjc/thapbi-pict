.. _microalgal_mock_sample_data:

Marker data
===========

Either clone the THAPBI PICT source code repository, or decompress the
latest source code release (``.tar.gz`` file). You should find it contains
a directory ``examples/microalgal_mock/`` which is for this example.

Shell scripts ``setup.sh`` and ``run.sh`` should reproduce the analysis
discussed.

Subdirectories ``V4/`` and ``V8V9/`` are used for the different amplicons
(with different primer settings and reference databases).

FASTQ data
----------

File ``PRJNA314977.tsv`` was download from the ENA and includes the FASTQ
checksums, URLs, and sample metadata. File ``metadata.tsv`` contains metadata
about the samples (see below).

Script ``setup.sh`` will download the pre-trimmed FASTQ files for Bradley
*et al.* (2016) from https://www.ebi.ac.uk/ena/data/view/PRJNA314977

It will download 248 raw FASTQ files (124 pairs), taking 1.4GB on disk.

If you have the ``md5sum`` tool installed (standard on Linux), verify the FASTQ
files downloaded correctly:

.. code:: console

    $ cd raw_data/V4/
    $ md5sum -c MD5SUM.txt
    ...
    $ cd ../../

.. code:: console

    $ cd raw_data/V8V9/
    $ md5sum -c MD5SUM.txt
    ...
    $ cd ../../

There is no need to decompress the files.

Amplicon primers & reference sequences
--------------------------------------

Referring to Table 1 in the paper, for the V4 target region they used
Reuk454FWD1/V4r primer pair (``CCAGCASCYGCGGTAATTCC`` and
``ACTTTCGTTCTTGAT``), while for the V8-V9 target they used the V8f/1510r pair
(``ATAACAGGTCTGTGATGCCCT`` and ``CCTTCYGCAGGTTCACCTAC``). However, the FASTQ
files provided have already been primer trimmed.

This means we need to run THAPBI PICT twice (once for each primer pair,
against a different marker database each time).

Metadata
--------

The provided file ``metadata.tsv`` is based on Supplementary Material Table S2,
although consistently using just one row per sample (representing three replicates),
and with the location and mock community ratio split into separate columns. This
has been combined with the MiSeq accessions from ``PRJNA314977.tsv``. It has seven
columns:

1. Environment, "Freshwater", "Wastewater", "Marine", "Control" or "Mock community"
2. Sample, e.g. "3F" or "MC4" with letters from environment, without replicate suffix
3. Description, e.g. "marsh" or "ocean"
4. Location, e.g. "Bray's Bayou, TX"
5. Freshwater:Marine, ratio for mock communities, e.g. "1:100"
6. V4 MiSeq, semi-colon separated list of three SRR accessions.
7. V8V9 MiSeq, semi-colon separated list of three SRR accessions.

Note that rather than providing separate files for the V4 and V8V9 sequencing,
they have been given separate cross-reference columns. Where a replicate is
missing, "Undetermined" is listed rather than an SRR accession.

When calling THAPBI PICT, the meta data commands are given as follows for V4:

.. code:: console

    $ thapbi_pict ... -t metadata.tsv -c 1,2,3,4,5 -x 6

And as follows for V8V9:

.. code:: console

    $ thapbi_pict ... -t metadata.tsv -c 1,2,3,4,5 -x 7

Argument ``-t metadata.tsv`` says to use this file for the metadata.

Argument ``-c 1,2,3,4,5`` says which columns to display and sort by. This means
Environment, Sample, Description, Location, Freshwater:Marine. This sorting
groups the samples logically (although it would be preferable to have the mock
community listed first or last, and "10W" appears before "2W").

Argument ``-x 6`` or ``-x 7`` indicates the V4 or V8V9 filename stems can be
found in that column respectively.

We have not given a ``-g`` argument to assign colour bands in the Excel
reports, so it will default to the first column in ``-c``, meaning we get
coloured bands for "Freshwater", "Wastewater", "Mock community" and "Marine".

Other files
-----------

File ``mock_community.fasta`` contains the sequences from accession numbers
KU900218 to KU900229 (published with the paper), with the description line
holding just the assigned species.

Files ``mock_community.known.tsv``, ``mock_freshwater.known.tsv`` and
``mock_marine.known.tsv`` describe the expected 12 species in the mock
community (and six species mixes at the purely freshwater or purely marine
extremes).
