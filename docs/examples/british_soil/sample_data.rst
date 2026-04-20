.. _british_soil_sample_data:

Marker data
===========

Either clone the THAPBI PICT source code repository, or decompress the
latest source code release (``.tar.gz`` file). You should find it contains
a directory ``examples/british_soil/`` which is for this example.

Shell scripts ``setup.sh`` and ``run.sh`` should reproduce the analysis
discussed.

FASTQ data
----------

File ``PRJNA690943.tsv`` was download from the ENA and includes the FASTQ
checksums, URLs, and sample metadata for the ITS1 sequencing in Landa *et al.*
(2021). This was cross referenced with Table 1 from the paper to create a
:ref:`sample metadata file <metadata>` file ``metadata.tsv`` for use with
THAPBI PICT (see below).

Script ``setup.sh`` will download the raw FASTQ files from
https://www.ebi.ac.uk/ena/data/view/PRJNA690943

It will download 202 FASTQ files (101 pairs) of ITS1 data taking about 200MB
on disk.

If you have the ``md5sum`` tool installed (standard on Linux), verify the
FASTQ files downloaded correctly:

.. code:: console

    $ cd raw_data/
    $ md5sum -c MD5SUM.txt
    $ cd ..

There is no need to decompress the files.

Amplicon primers & reference sequences
--------------------------------------

The ITS primers used this this paper match the THAPBI PICT defaults, so our
default ITS1 database can be used. The paper also sequenced a region of COI,
but a meaningful analysis of that data with THAPBI PICT would require
considerable effort to make a reference database of full length sequences.

Metadata
--------

The provided file ``metadata.tsv`` has seven columns from the author's Table 1
and an additional column added with the ITS1 Illumina accessions:

1. Site Code (Site Type), e.g. "DS1 (Arboreta)"
2. Location, e.g. "DS1-1"
3. Surveyed Location and Nearest Hosts, e.g. "Mixed Pinus species and Abies"
4. Signs of Ill-Health Tree Status, e.g. "Healthy"
5. ITS, what the authors found, e.g. "P. austrocedri"
6. COI, what the authors found, e.g. "Clade 8 Phytophthora sp. uncultured 5b"
7. Isolation, what the authors found, e.g. "P. cryptogea"
8. ITS Illumina, the public archive accession, e.g. "SRR13393769"

Additionally rows have been added for the positive and negative controls.

When calling THAPBI PICT, the meta data commands are given as follows:

.. code:: console

    $ thapbi_pict ... -t metadata.tsv -x 8 -c 1,2,3,4,5,6,7

Argument ``-t metadata.tsv`` says to use this file for the metadata.

The ``-x 8`` argument indicates the filename stem can be found in column 8.

Argument ``-c 1,2,3,4,5,6`` says which columns to display and sort by (do
not include the indexed column again). If for example the accession was
listed first, it would be sorted on that, which is not helpful here.

We have not given a ``-g`` argument to assign colour bands in the Excel
reports, so it will default to the first column in ``-c``, meaning we get
coloured bands for the site.

Other files
-----------

File ``mock_community.known.tsv`` (plain text tab separated table) describes
the expected 10 species in the mock community positive controls.
