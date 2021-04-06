.. _custom_database_sample_data:

Marker data
===========

Either clone the THAPBI PICT source code repository, or decompress the
latest source code release (``.tar.gz`` file). You should find it contains
a directory ``examples/recycled_water/`` which is for this example.

Shell scripts ``setup.sh`` and ``run.sh`` should reproduce the analysis
discussed.

The documentation goes through running each step of the analysis gradually,
including building a custom database, before finally calling pipeline command
to do it all together. We provide script ``run.sh`` to do the final run-though
automatically, but encourage you to follow along the individual steps first.

FASTQ data
----------

File ``PRJNA417859.tsv`` was download from the ENA and includes the FASTQ
checksums, URLs, and sample metadata. With a little scripting to extract the
relevant :ref:`sample metadata <metadata>` for use with THAPBI PICT this was
reformatted as ``metadata.tsv`` (see below).

Script ``setup.sh`` will download the raw FASTQ files for Redekar *et al.*
(2019) from https://www.ebi.ac.uk/ena/data/view/PRJNA417859 - you could also
use https://www.ncbi.nlm.nih.gov/bioproject/PRJNA417859/

It will download 768 raw FASTQ files (384 pairs), taking about 4.8GB on disk

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

A region of ITS1 was amplified using the ITS6/ITS7 primer pair
(``GAAGGTGAAGTCGTAACAAGG`` and ``AGCGTTCTTCATCGATGTGC``) which bind the
5.8S rDNA, described here:

    Cooke *et al.* (2000) A molecular phylogeny of Phytophthora and related
    oomycetes. https://doi.org/10.1006/fgbi.2000.1202

The left primer (ITS6) matches the THAPBI PICT default, but their right primer
(ITS7) matches about 60bp further downstream in *Phytophthora*. This means we
can use THAPBI PICT default settings and get meaningful but blinkered results
(for the subset of the data which our narrower primer set would have amplified,
using a *Phytophthora* centric database).

In order to classify beyound *Phytophthora*, we need to build a THABPI PICT
database including *Pythium* and *Phytopythium*. Redekar *et al.* (2019)
Supplementary Table 3 provides a list of 1454 unique accessions and the
species they assigned to it (not always the same as that listed on the NCBI
record, as those annotations can change). Looking at those sequences, bar
a handful they extend though the right primer. However, only about 50 have
the left primer sequence included (depending how stringent you are), and
the rest are also missing the next 32bp.

The ITS6 primer is situated within a highly conserved region, and the next
32bp is highly conserved, usually ``TTTCCGTAGGTGAACCTGCGGAAGGATCATTA``.
Unfortunately, the majority of published *Oomycetes* ITS1 sequences omit
this. For the curated *Phytophthora* in the THAPBI PICT default database,
we have inserted the expected sequence - and have yet to find a counter
example. However, Redekar *et al.* (2019) took the other obvious choice, and
remove it from their reads:

    *trimming extra bases from read1: an additional 32 bases from the 5′ end
    of read1, which mapped to 18S segment, were trimmed as the oomycete ITS
    reference database does not include the 18S segment;*

We can do something similar in THAPBI PICT by treating this typically
conserved 32bp region as part of the left primer - requiring it be present
(while allowing some ambiguity) and removing it - leaving a shorter fragment
which can be matched to a database built of those 1454 accessions.

Metadata
--------

The provided file ``metadata.tsv`` has seven columns:

1. Source, "Reservoir", "River" or "Runoff"
2. Site,  "A", "B", "C", ..., "M"
3. Process, "Filtration" or "Leaf baiting"
4. Period, "01" to "28"
5. Year-Month, "2015-04" to "2016-05" (given as "YYYY-MM" for sorting)
6. Sample, author's sample name, e.g. "OSU484"
7. Accession, assigned by the public archive, e.g. "SRR6303585"

When calling THAPBI PICT, the meta data commands are given as follows:

.. code:: console

    $ thapbi_pict ... -t metadata.tsv -x 7 -c 1,2,3,4,5,6

Argument ``-t metadata.tsv`` says to use this file for the metadata.

The ``-x 7`` argument indicates the filename stem can be found in column 7,
Accession.

Argument ``-c 1,2,3,4,5,6`` says which columns to display and sort by (do
not include the indexed column again). If for example the accession was
listed first, it would be sorted on that, which is not helpful here. If you
prefer to sort on site first, or by date before process, this should be
straightforward.

We have not given a ``-g`` argument to assign colour bands in the Excel
reports, so it will default to the first column in ``-c``, meaning we get
three coloured bands for "Reservoir", "River" and "Runoff".

Other files
-----------

Files ``Redekar_et_al_2019_sup_table_3.tsv`` (plain text tab separated table)
and ``Redekar_et_al_2019_sup_table_3.fasta`` (FASTA format) are based on the
Excel format Supplementary Table 3 from the paper.
