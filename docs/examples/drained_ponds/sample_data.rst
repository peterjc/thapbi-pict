.. _drained_ponds_sample_data:

Marker data
===========

Either clone the THAPBI PICT source code repository, or decompress the
latest source code release (``.tar.gz`` file). You should find it contains
a directory ``examples/drained_ponds/`` which is for this example.

Shell scripts ``setup.sh`` and ``run.sh`` should reproduce the analysis
discussed.

FASTQ data
----------

File ``PRJNA638011.tsv`` was download from the ENA and includes the FASTQ
checksums, URLs, and sample metadata. Related file ``metadata.tsv`` combines
this with metadata about the samples from the paper (see below).

Script ``setup.sh`` will download the raw FASTQ files for Muri *et al.*
(2020) from https://www.ebi.ac.uk/ena/data/view/PRJNA638011 - you could also
use https://www.ncbi.nlm.nih.gov/bioproject/PRJNA638011/

It will download 198 raw FASTQ files (99 pairs), taking about 550MB on disk

If you have the ``md5sum`` tool installed (standard on Linux), verify the
FASTQ files downloaded correctly:

.. code:: console

    $ cd raw_data/
    $ md5sum -c MD5SUM.txt
    $ cd ..

There is no need to decompress the files.

Amplicon primers & reference sequences
--------------------------------------

A region of 12S was amplified using a previously published primer pair
(``ACTGGGATTAGATACCCC`` and ``TAGAACAGGCTCCTCTAG``) described here:

    Kelly *et al.* (2014) Understanding PCR processes to draw meaningful
    conclusions from environmental DNA studies.
    https://doi.org/10.1038/s41598-019-48546-x

This primer amplifies not just the fish of interest, but also birds and
mammals (including human).

Rather than trying to use the same curated database from the University of
Hull Evolutionary and Environmental Genomics Group, who also wrote metaBEAT
(metaBarcoding and Environmental Analysis Tool) which was used in the paper,
we provide a crudely curated database culled from stringent BLASTN matches in
the NCBI NT database (search run with 100% query coverage and 99% identity,
see provided ``scripts/blast_to_fasta.py``), as file ``NCBI_12S.fasta``. The
``run.sh`` script starts by loading this into a new THAPBI PICT database.

Metadata
--------

The provided file ``metadata.tsv`` has ten columns, the first three are from
``PRJNA638011.tsv`` (ENA metadata) and the rest from the paper's Supplementary
Table S2 - cross referenced on the sample name/alias:

1. run_accession, from ENA metadata, e.g. "SRR11949879"
2. sample_alias, from ENA metadata, e.g. "Lib3-M3-1F1"
3. sample_title, from ENA metadata, e.g. "MCE-Sample 1 filter 1"
4. samples, from Supplementary Table S2, e.g. "M3-1F1"
5. lake, from Table S2, e.g. "Middle_lake"
6. filter, from Table S2, e.g. "MCE"
7. treatment, from Table S2, e.g. "F1"
8. extracted, from Table S2, e.g. "Filter"
9. control, from Table S2, either "", "blank", "negative", or "positive"
10. date, from Table S2, e.g. "17.02.2017"

When calling THAPBI PICT, the meta data commands are given as follows:

.. code:: console

    $ thapbi_pict ... -t metadata.tsv -x 1 -c 5,6,7,8,9,10,4,3

Argument ``-t metadata.tsv`` says to use this file for the metadata.

The ``-x 1`` argument indicates the filename stem can be found in column 1,
the ENA assigned run accession.

Argument ``-c 5,6,7,8,9,10,4,3`` says which columns to display and sort by (do
not include the indexed column again). If for example the accession was
listed first, it would be sorted on that, which is not helpful here. If you
prefer to sort on filter first, that change should be straightforward.

We have not given a ``-g`` argument to assign colour bands in the Excel
reports, so it will default to the first column in ``-c``, meaning we get
three coloured bands for "Middle_lake", "NA" (controls), and "New_lake".

Other files
-----------

Files ``cichlid_control.known.tsv`` and ``negative_control.known.tsv`` and are
used in ``setup.sh`` to create ``expected/*.known.tsv`` entries for the
positive and negative controls, including the blank controls.

12 fish species were translocated to New Lake, of which nine were also in the
middle lake. Referring to the results text and Figure 1B, and pooling the two
observed hybrids with a parent species, the expected species in the two lakes
are as follows.

Middle lake and new lake both had:

 - *Abramis brama*
 - *Barbus barbus*
 - *Carassius carassius*
 - *Cyprinus carpio*
 - *Perca fluviatilis*
 - *Rutilus rutilus*
 - *Scardinius erythrophthalmus*
 - *Squalius cephalus*
 - *Tinca tinca*

New lake only also had:

 - *Acipenser* spp.
 - *Ctenopharyngodon idella*
 - *Silurus glanis*

File ``middle_lake.known.tsv`` lists the 9 species found in the middle lake,
and ``new_lake.known.tsv`` lists the 12 species in the new lake (although not
all fish are expected at all sites within each lake), and these are assigned
to the remaining samples as ``expected/*.known.tsv`` by running ``setup.sh``.
