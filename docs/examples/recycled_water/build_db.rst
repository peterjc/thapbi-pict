.. _custom_database_building:

Different primers
=================

This example is based on the following paper from another research group:

* Redekar *et al.* (2019) Diversity of *Phytophthora*, *Pythium*, and
  *Phytopythium* species in recycled irrigation water in a container nursery.
  https://doi.org/10.1094/PBIOMES-10-18-0043-R

The worked example starts by running the pipeline command with :ref:`default
settings <custom_database_defaults>`, which uses the default *Phytophthora*
centric database and primers. We can do that because the tool's default
database defines primers which target a subset of the longer amplicon
amplified in this example dataset.

Now we will change the primer settings. Using the actual right primer will
extend the *Phytophthora* FASTA sequences about 60bp (and accept many more
non-*Phytophthora*). The left primer is actually the same, but to match the
analysis and references from Redekar *et al.* (2019), we want to trim off the
typically conserved 32bp fragment ``TTTCCGTAGGTGAACCTGCGGAAGGATCATTA`` from
the start of each amplicon, which we can do by pretending this is part of the
left primer.

The up-shot is by cropping about 32bp off the start, and adding about 60bp
at the end, we will no longer get any matches against the default database
with the default classifier (it is too strict, the matches are too distant).

This means before we can run the entire pipeline, we will need to build a
custom database. We'll discuss the sequences which go into this database
next, but this use ``--marker ITS1-long`` to name this alternative marker,
and set ``--left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA``
and ``--right AGCGTTCTTCATCGATGTGC`` to declare the primers.

Building a custom database
==========================

We will build a database using the same public sequences as Redekar *et al.*
(2019), using the accessions given in Supplementary Table 3 in Excel format.

We have taken their list of accessions and species names (ignoring voucher or
isolate numbers), edited some punctuation to match the NCBI taxonomy, added
some missing accession version suffixes, deduplicated, and made a simple
tab-separated plain text table, with 1454 entries. In the setup instructions
for this example you should have got a copy of this file, named
``Redekar_et_al_2019_sup_table_3.tsv``, and a matching FASTA file
``Redekar_et_al_2019_sup_table_3.fasta`` which we will import into the new
database.

This table is sorted alphabetically by species then accession, and starts:

.. code:: console

    $ head Redekar_et_al_2019_sup_table_3.tsv
    <SEE TABLE EXCERPT BELOW>

You could also look at the TSV file in Excel:

========== ===================
HQ643082.1 Achlya ambisexualis
HQ643083.1 Achlya ambisexualis
HQ643084.1 Achlya americana
HQ643085.1 Achlya aquatica
HQ643086.1 Achlya bisexualis
HQ643087.1 Achlya bisexualis
HQ643088.1 Achlya bisexualis
HQ643089.1 Achlya caroliniana
HQ643090.1 Achlya colorata
HQ643091.1 Achlya colorata
========== ===================

Determining the species
-----------------------

Consider ``FJ666127.1`` which Redekar *et al.* (2019) listed as *Phytophthora
aquimorbida* - yet at the time of writing, the file downloaded from
https://www.ebi.ac.uk/ena/browser/api/fasta/FJ666127.1 is as follows, with
a species name of *Phytophthora* sp. CCH-2009b::

    >ENA|FJ666127|FJ666127.1 Phytophthora sp. CCH-2009b isolate 40A6 internal transcribed spacer 1, partial sequence; 5.8S ribosomal RNA gene, complete sequence; and internal transcribed spacer 2, partial sequence.
    CCACACCTAAAAACTTTCCACGTGAACTGTCTGTGATGTTGGGGGGCTGCTGCTGCTGCT
    TCGGTGGCGGCGTGCTCCCATCAAACGAGGCCCTGGGCTGCAAAGTCGGGGGTAGTAGTT
    ACTTTTTGTAAACCCTTTTCCTGTATTTTCTGAATATACTGGGGGGACGAAAGTCTCTGC
    TTTTAACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTCGCACATCGATGAAGAACG
    CTGCGAACTGCGATACGTAATGCGAATTGCAGGATTCAGTGAGTCATCGAAATTTTGAAC
    GCATATTGCACTTCCGGGTTATGCCTGGGAGTATGCCTGTATCAGTGTCCGTACATCAAT
    CTTGGCTTCCTTCCTTCCGTGTAGTCGGTGGCGGGAACGCGCAGACGTGAAGTGTCTTGC
    CTGTGGCTCCAGCTGTTGTTGGGGTGGTGTGGGCGAGTCCTTTGAAATGTAAGATACTGT
    TCTTCTCTTTGCTGGAAAAGCGTGCGCTGTGTGGTTGTGGAGGCTGCCGTGGTGGCCAGT
    CGGCGACTGACTTCGTGCTGATGCGTGTGGAGAGGCTCTGGATTCGCGGTATGGTTGGCT
    TCGGCTGAACTTCTGCTTATTTGGGTGTCTTTTCGCTGCGTTGGCGTGTCGGGGTTGGTG
    AACCGTAGTCATTTCGGCTTGGCTTTTGAACCGCGTGGCTGTAGCGCGAAGTATGGCGGC
    TGCCTTTGTGGCGGCCGAGAGGACGACCTATTTGGGACGATTGTGCGGCCTCGTGCTGCA
    TCTCAA

Notice that the species name runs into the general description, which
is problematic. Unless THAPBI PICT has a pre-loaded taxonomy to use
for validation, it has to use heuristics to split up this long string -
which is not fully reliable.

If we look at https://www.ncbi.nlm.nih.gov/nucleotide/FJ666127.1 on the
NCBI website, we see it in GenBank format which is a little different::

    LOCUS       FJ666127                 786 bp    DNA     linear   PLN 09-MAR-2009
    DEFINITION  Phytophthora sp. CCH-2009b isolate 40A6 internal transcribed spacer
                1, partial sequence; 5.8S ribosomal RNA gene, complete sequence;
                and internal transcribed spacer 2, partial sequence.
    ACCESSION   FJ666127
    VERSION     FJ666127.1
    KEYWORDS    .
    SOURCE      Phytophthora aquimorbida
      ORGANISM  Phytophthora aquimorbida
                Eukaryota; Stramenopiles; Oomycetes; Peronosporales;
                Peronosporaceae; Phytophthora.
    ...

The NCBI metadata has the species *Phytophthora aquimorbida* separate
from the author submitted description which starts with an older name,
"Phytophthora sp. CCH-2009b" - which is in fact listed as an alias on
the NCBI taxonomy database under `taxonomy ID 611798
<https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=611798>`_.

THAPBI PICT offers two solutions. By default the *entire* FASTA description
(after the identifier) is the species name, giving full control to the user.

However, ``-c ncbi`` switches on NCBI heuristics. This is best used with a
pre-loaded NCBI taxonomy in the database for validation purposes. This tries
as many words as possible from the NCBI style FASTA description in looking for
a match in the NCBI taxonomy, including synonyms. If that fails and lax mode
is used (``-x`` or ``--lax``), it falls back on heuristics to identify which
part of the description is the species.

Species validation
------------------

THAPBI PICT by default validates imports against the NCBI taxonomy, and
that includes support for known synonyms. This requires downloading the
taxonomy files and running the ``thapbi-pict load-tax`` command.

The NCBI currently provide their taxonomy dump in two formats, old and new.
THAPBI PICT supports both, we'll use the old format as the download is half
the size - we only need the ``names.dmp``, ``nodes.dmp`` and ``merged.dmp``
files:

.. code:: console

    $ curl -L -O https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2019-12-01.zip
    ...
    $ unzip -n -d taxdmp_2019-12-01 taxdmp_2019-12-01.zip
    ...
    $ ls -1 taxdmp_2019-12-01/*.dmp
    taxdmp_2019-12-01/citations.dmp
    taxdmp_2019-12-01/delnodes.dmp
    taxdmp_2019-12-01/division.dmp
    taxdmp_2019-12-01/gencode.dmp
    taxdmp_2019-12-01/merged.dmp
    taxdmp_2019-12-01/names.dmp
    taxdmp_2019-12-01/nodes.dmp

Building the database becomes a two-step process, first importing the
taxonomy, and second importing the sequences.

If you are working with different organisms you will also need to set the
``-a`` or ``--ancestors`` option which defaults to `NCBI taxonomy ID 4762
<https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4762>`_ for
*Oomycetes*.

Primer trimming
---------------

We have provided file ``Redekar_et_al_2019_sup_table_3.fasta`` which contains
primer trimmed versions of the full sequences of each accession, plus the
species name from ``Redekar_et_al_2019_sup_table_3.tsv`` which was based on
those given in Redekar *et al.* (2019) Supplementary Table 3 but with some
light curation to better match the NCBI usage. Note that matching sequences
have been combined into single FASTA records with a semi-colon separated
description.

The sequencing trimming ought to be very close to that used in the Redekar
*et al.* (2019) paper's database. This file was constructed with a short Python
script parsing the information in ``Redekar_et_al_2019_sup_table_3.tsv`` and
the downloaded full sequences.
Then ``cutadapt -g GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA ...``
found and removed 64 left prefixes. This was followed by running
``cutadapt -a GCACATCGATGAAGAACGCT ...`` which trimmed 1439 sequences (99.9%)
and warned that the "adapter" might be incomplete because the sequence
preceding it was highly conserved. That left 1451 sequences, but with many
duplicates. This was made non-redundant giving 841 unique sequences with
de-duplicated entries recorded with semi-colon separated FASTA title lines.

Now, let's load the FASTA file into a new THAPBI PICT database with the NCBI
taxonomy pre-loaded (which will enable synonym support), but not enforced
(``-x`` or ``--lax`` mode). We'll name the new marker "ITS1-long" and record
the left and right primers which will be used later when processing the reads:

.. code:: console

    $ rm -rf Redekar_et_al_2019_sup_table_3.sqlite  # remove it if already there
    $ thapbi_pict load-tax -d Redekar_et_al_2019_sup_table_3.sqlite -t taxdmp_2019-12-01/
    ...
    $ thapbi_pict import -d Redekar_et_al_2019_sup_table_3.sqlite \
      --lax --sep ";" -i Redekar_et_al_2019_sup_table_3.fasta \
      --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA \
      --right AGCGTTCTTCATCGATGTGC --marker ITS1-long
    File Redekar_et_al_2019_sup_table_3.fasta had 841 sequences, of which 838 accepted.
    Of 1451 potential entries, loaded 1451 entries, 0 failed parsing.

Just a few short sequences were rejected - giving in total 1451 entries.
The vast majority are recorded with an NCBI taxid, just four exceptions
(visible if you run the last command with ``-v`` or ``--verbose``):

- *Phytophthora taxon aquatilis* from
  `FJ666126.1 <https://www.ncbi.nlm.nih.gov/nucleotide/FJ666126.1>`_,
  which the NCBI say should be *Phytophthora* sp. CCH-2009a
- *Phytophthora fragaefolia* from
  `AB305065.1 <https://www.ncbi.nlm.nih.gov/nucleotide/AB305065.1>`_,
  which the NCBI say should be *Phytophthora fragariaefolia*.
- *Phytophthora citricola sensu stricto* from
  `FJ560913.1 <https://www.ncbi.nlm.nih.gov/nucleotide/FJ560913.1>`_,
  which the NCBI say should be just *Phytophthora citricola*.
- *Phytopythium sp. amazonianum* from
  `HQ261725.1 <https://www.ncbi.nlm.nih.gov/nucleotide/HQ261725.1>`_,
  which the NCBI say should be *Pythium* sp. 'amazonianum'.

None of these are clear cut (there were a lot more conflicts, mostly down to
differences in punctuation, already addressed in preparing the TSV and FASTA
file).

If you left off the ``-x`` (or ``--lax``) option, those four would not have
been imported into the database.

Taxonomic conflicts
-------------------

The ITS1 region is not ideal as a barcode sequence.  In the *Phytophthora*
there are many cases where the same marker is shared by multiple species.
The ``thapbi_pict conflicts`` command is provided to check for this, or
worse -- conflicts at genus level:

.. code:: console

    $ thapbi_pict conflicts -h
    ...

Let's run this on the custom database, with output to a file:

.. code:: console

    $ thapbi_pict conflicts -d Redekar_et_al_2019_sup_table_3.sqlite -o conflicts.tsv; echo "(Return code $?)"
    (Return code 3)

Command line tools use a non-zero return code by convention to indicate an
error. Here we return the number of genus level conflicts, three, as can be
seen by looking at the start of the plain text tab separated table output:

.. code:: console

    $ head -n 5 conflicts.tsv
    #MD5                              Level    Conflicts
    87e588784b04ba5f4538ff91acb50c0f  genus    Lagenidium;Pythium
    9bb2ab5b9f88256516f2ae618c16a62e  genus    Brevilegnia;Globisporangium
    ff35f216832110904cc6fd1c9def33fd  genus    Achlya;Saprolegnia
    077ae505c0ad210aa4c071417a4f2f9a  species  Saprolegnia monilifera;Saprolegnia unispora

There are lots species level conflicts, some of which might be subspecies etc.
However, more concerning is three genus level conflicts.

One way to see which accessions are a problem is filtering the dump command
output (introduced properly in :ref:`custom_database_examine`), e.g.

.. code:: console

    $ thapbi_pict dump -d Redekar_et_al_2019_sup_table_3.sqlite \
      | cut -f 2-6 | grep 87e588784b04ba5f4538ff91acb50c0f
    HQ643136.1  Lagenidium  caudatum   135481  87e588784b04ba5f4538ff91acb50c0f
    HQ643539.1  Pythium     flevoense  289620  87e588784b04ba5f4538ff91acb50c0f
    Wrote 1451 txt format entries

Some could be mislabelled, for ``9bb2ab5b9f88256516f2ae618c16a62e`` we see
the vast majority are *Globisporangium ultimum* with just one sequence
`HQ643127.1 <https://www.ncbi.nlm.nih.gov/nucleotide/HQ643127.1>`_ labelled
as *Brevilegnia gracilis*:

.. code:: console

    $ thapbi_pict dump -d Redekar_et_al_2019_sup_table_3.sqlite \
      | cut -f 3-6 | grep 9bb2ab5b9f88256516f2ae618c16a62e \
      | sort | uniq -c | sed 's/^ *//g'
    1 Brevilegnia       gracilis  944588   9bb2ab5b9f88256516f2ae618c16a62e
    42 Globisporangium  ultimum   2052682  9bb2ab5b9f88256516f2ae618c16a62e

Checking the current NCBI annotation of these accessions does not suggest
problems with recent taxonomy changes like *Phytopythium* vs *Pythium*.

Those assignments might have changed since this was written. Taxonomy is
fluid, so if using any single authority, make sure to document which version
(e.g. month and year for the NCBI taxonomy).
