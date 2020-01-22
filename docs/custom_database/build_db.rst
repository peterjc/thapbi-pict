.. _custom_database_building:

Building a custom database
==========================

This example is based on the following paper from another research group:

* Redekar *et al* (2019) Diversity of *Phytophthora*, *Pythium*, and
  *Phytopythium* species in recycled irrigation water in a container nursery.
  https://doi.org/10.1094/PBIOMES-10-18-0043-R

We will build a database using the same public sequences that they used,
accessions given in Supplementary Table 3 in Excel format. We will need to
collect these sequences in FASTA format.

We have taken their list of accessions and species names (ignoring voucher or
isolate numbers), added missing accession versions as needed, deduplicated,
and made a simple tab-separated plain text table, with 1454 entries. In the
setup instructions for this example you should have got a copy of this file,
named ``Redekar_et_al_2019_sup_table_3.tsv``.

This file is sorted alphabetically by species then accession, and starts:

========== ===================
Accession  Species
---------- -------------------
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

Using species given in annotation
---------------------------------

The following snippet of code will download the 1454 accessions in FASTA
format from the EBI - a similar approach could be used to download them
from the NCBI:

.. code:: console

    $ mkdir cache
    $ cd cache
    $ for ACC in `cut -f 1 ../Redekar_et_al_2019_sup_table_3.tsv`; do curl -L -o "$ACC.fasta" "https://www.ebi.ac.uk/ena/browser/api/fasta/$ACC"; done
    ...
    $ ls -1 *.fasta | wc -l
    1454
    $ cd ..

We will now turn that folder of FASTA files into a single combined file:

.. code:: console

    $ rm -rf Redekar_et_al_2019_sup_table_3_download.fasta
    $ for ACC in `cut -f 1 Redekar_et_al_2019_sup_table_3.tsv`; do \
      cat "cache/$ACC.fasta" >> Redekar_et_al_2019_sup_table_3_download.fasta; done
    $ grep -c "^>" Redekar_et_al_2019_sup_table_3_download.fasta
    1454

We can now load this 1454 sequence FASTA file into a new THAPBI PICT database:

.. code:: console

    $ rm -rf first_attempt.sqlite  # delete if already there
    $ thapbi_pict ncbi-import -d first_attempt.sqlite \
      -i Redekar_et_al_2019_sup_table_3_download.fasta -x \
      --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA \
      --right AGCGTTCTTCATCGATGTGC
    File Redekar_et_al_2019_sup_table_3_download.fasta had 1454 sequences, of which 1451 accepted.
    Of 1451 potential entries, loaded 1451 entries, 0 failed parsing.

So far so good - but only 1451 of the 1454 sequences were loaded.
If you run that again with ``-v`` or ``--verbose`` you will see the missing
three sequences are all *Pythium heterothallicum* entries which were shorter
than the default minimum length.

The ``-x`` (or ``--lax``) turns off species name validation, we'll come back
to this later.

Notice we specified the left and right primers (as discussed in :ref:`primers`).
In this context THAPBI PICT will look for and discard the left and right
primers in isolation - neither has to be present. This is to cope with
pre-trimmed entries, or cases like the *Oomycetes* markers where the
published sequence typically start immediately after our left "primer"
(recall ``GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA`` is
actually the PCR primer ``GAAGGTGAAGTCGTAACAAGG`` plus the conserved
32bp fragment ``TTTCCGTAGGTGAACCTGCGGAAGGATCATTA``).

To disable this you can alternatively set the primers to empty strings,
or use the similar ``thapbi_pict curated-import`` command which assumes
pre-trimmed inputs and handles determining the species differently.

Determining the species
-----------------------

Consider ``FJ666127.1`` which Redekar *et al.* listed as *Phytophthora
aquimorbida* - yet at the time of writing, the file downloaded from
https://www.ebi.ac.uk/ena/browser/api/fasta/FJ666127.1 is as follows::

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

THAPBI PICT offers two solutions. First, the import commands by default
expect a pre-loaded NCBI taxonomy in the database for validation purposes,
and this includes synoym support. This allows ``thapbi_pict ncbi-import``
to try as many words as possible from the FASTA description in looking
for a match in the NCBI taxonomy. In the example above we disabled this
with ``-x`` (or ``--lax``). Second, ``thapbi_pict curated-import`` takes
the *entire* FASTA description (after the identifier) as the species name.

Species validation
------------------

THAPBI PICT by default validates imports against the NCBI taxonomy, and
that includes support for known synonyms. This requires downloading the
taxonomy files and running the ``thapbi-pict load-tax`` command.

The NCBI currently provide their taxonomy dump in two formats, old and new.
THAPBI PICT supports both, we'll use the old format as the download is half
the size - we only need the ``names.dmp`` and ``nodes.dmp`` files:

.. code:: console

    $ curl -L -O https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/taxdmp_2020-01-01.zip
    $ unzip -d taxdmp_2020-01-01 taxdmp_2020-01-01.zip
    ...
    $ ls taxdmp_2020-01-01/n*.dmp
    taxdmp_2020-01-01/names.dmp  taxdmp_2020-01-01/nodes.dmp

Now building the database is a two-step process, first importing the
taxonomy, and second importing the sequences.

If you are working with different organisms you will also need to set the
``-a`` or ``--ancestors`` option which defaults to `NCBI taxonomy ID 4762
<https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=4762>`_ for
*Oomycetes*.

.. code:: console

    $ rm -rf with_validation.sqlite  # remove it if already there
    $ thapbi_pict load-tax -d with_validation.sqlite -t taxdmp_2020-01-01/
    WARNING: Treating species group 'Hyaloperonospora parasitica species group' (txid453155) as a species.
    WARNING: Genus Elongisporangium (1448050) has no children
    Loaded 81 new genera, and 1028 new species entries with 2636 synonyms into DB (0, 0 and 7 already there)
    $ thapbi_pict ncbi-import -d with_validation.sqlite \
      -i Redekar_et_al_2019_sup_table_3_download.fasta \
      --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA \
      --right AGCGTTCTTCATCGATGTGC
    File Redekar_et_al_2019_sup_table_3_download.fasta had 1454 sequences, of which 1444 accepted.
    Of 1451 potential entries, 0 unparsable, 7 failed sp. validation, 1444 OK.
    Could not validate 3 different species names

Notice this time we run ``thapbi_pict ncbi-import`` without the ``-x`` (``--lax``)
option, and it complains about three species names and seven entries - but which?
If you repeat this but add ``-v`` or ``--verbose`` to the import command you can
see:

- *Phytophthora sansomea* from five accessions including
  `HQ261667.1 <https://www.ncbi.nlm.nih.gov/nucleotide/HQ261667.1>`_,
  where it looks like it a typo for *Phytophthora sansomeana* (although the
  full story here is likely more complicated - both existed in the NCBI
  taxonomy of December 2019 as taxid 358102 and 555429 respectively, but were
  merged in the Janurary 2020 release *without* defining a synonym).
- *Phytophthora lagoariana* from
  `EF590256.1 <https://www.ncbi.nlm.nih.gov/nucleotide/EF590256.1>`_,
  which the NCBI says should be "*Phytophthora* sp. 'lagoariana'"
- *Phytophthora novaeguinee* from
  `EU035774.1 <https://www.ncbi.nlm.nih.gov/nucleotide/EU035774.1>`_,
  which the NCBI says should be "*Phytophthora* sp. *novaeguinee*"

Strict validation has its downsides when combined with uncurated metadata
and unrecorded synonyms.

In this case the problem wasn't actually splitting the species name from the
free text, but if it were, one fix would be to download the data in GenBank or
TinySeq XML format rather than FASTA, which would give the species separately.
Alternatively, THAPBI PICT will accept curated species data as described next.


Curated import
--------------

The ``thapbi_pict curated-import`` used above differs from the
``thapbi_pict ncbi-import`` command in two key points. First, by default it
expects the sequences to be pre-trimmed (it does not do primer trimming).
Second, it does not use heuristics but simply assumes the FASTA description
line is an identifier followed by the species name *only*.

We have provided file ``Redekar_et_al_2019_sup_table_3.fasta`` which contains
primer trimmed versions of the full sequences of each accession, plus the
species name as given in Redekar *et al.* 2019 Supplementary Table 3. This
ought to be very close to the trimmed sequences used in their database.

This was constructed with a short Python script parsing the information in
``Redekar_et_al_2019_sup_table_3.tsv`` and the downloaded full sequences.
Then ``cutadapt -g GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA ...``
found and removed only 64 left prefixes. This was followed by running
``cutadapt -a GCACATCGATGAAGAACGCT ...`` which trimmed 1439 sequences (99.9%)
and warned that the "adapter" might be incomplete because the sequence
preceeding it was highly conserved. That left 1451 sequences, but with many
duplicates. This was made non-redundant giving 841 unique sequences with
de-duplicated entries recorded with semi-colon separated FASTA title lines.

Now, let's load the FASTA file into a new THAPBI PICT database, without any
NCBI taxonomy validation:

.. code:: console

    $ rm -rf Redekar_et_al_2019_sup_table_3.sqlite  # remove it if already there
    $ thapbi_pict curated-import -x \
      -d Redekar_et_al_2019_sup_table_3.sqlite \
      -i Redekar_et_al_2019_sup_table_3.fasta
    File Redekar_et_al_2019_sup_table_3.fasta had 841 sequences, of which 838 accepted.
    Of 1451 potential entries, loaded 1451 entries, 0 failed parsing.

As above, three short sequences were rejected - giving in total 1451 entries.


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

    $ thapbi_pict conflicts -d Redekar_et_al_2019_sup_table_3.sqlite -o conflicts.tsv
    Loaded taxonomy for 838 sequences from DB

This produces a plain text tab separated table ``conflicts.tsv`` which starts
as follows:

================================ ======= =====================================
MD5                              Level   Conflicts
-------------------------------- ------- -------------------------------------
2a5b0b1f6404009e2987d004ff7e3a9c genus   Phytopythium;Pythium
3550a51c172b547e7626e8f99a942341 genus   Phytopythium;Pythium
87e588784b04ba5f4538ff91acb50c0f genus   Lagenidium;Pythium
9bb2ab5b9f88256516f2ae618c16a62e genus   Brevilegnia;Pythium
077ae505c0ad210aa4c071417a4f2f9a species Saprolegnia monilifera;Saprolegnia unispora
0966d89e2bcd49b6986db8231d7790bb species Phytophthora asparagi;Phytophthora taxon asparagi
...                              ...     ...
================================ ======= =====================================

There are 77 species level conflicts, some of which might be subspecies etc.
However, more concerning is four genus level conflicts - all involving *Pythium*.

One way to see which accessions are a problem is filtering the dump command
output, e.g.

.. code:: console

    $ thapbi_pict dump -d Redekar_et_al_2019_sup_table_3.sqlite \
      | cut -f 1-5 | grep 2a5b0b1f6404009e2987d004ff7e3a9c
    HQ643387.1  Phytopythium  litorale  2a5b0b1f6404009e2987d004ff7e3a9c
    HQ643677.1  Pythium       litorale  2a5b0b1f6404009e2987d004ff7e3a9c

In this case the NCBI has `HQ643387.1 <https://www.ncbi.nlm.nih.gov/nucleotide/HQ643387.1>`_
and `HQ643677.1 <https://www.ncbi.nlm.nih.gov/nucleotide/HQ643677.1>`_ both as
*Phytopythium litorale*, with *Pythium litoralis* as a homotypic synonym.

Likewise, using the same approach, ``3550a51c172b547e7626e8f99a942341`` is
from `HQ643393.1 <https://www.ncbi.nlm.nih.gov/nucleotide/HQ643393.1>`_ and
`HQ643394.1 <https://www.ncbi.nlm.nih.gov/nucleotide/HQ643394.1>`_ labelled
as *Phytopythium oedochilum* and
`HQ643712.1 <https://www.ncbi.nlm.nih.gov/nucleotide/HQ643712.1>`_ labelled
as *Pythium oedichilum* - which again the NCBI list as a homotypic synonym.

On the other hand, ``87e588784b04ba5f4538ff91acb50c0f`` is from
`HQ643136.1 <https://www.ncbi.nlm.nih.gov/nucleotide/HQ643136.1>`_ labelled
as *Lagenidium caudatum*
`HQ643539.1 <https://www.ncbi.nlm.nih.gov/nucleotide/HQ643539.1>`_ labelled
as *Pythium flevoense*
- and the NCBI still lists these as separate species in separate family.

While ``9bb2ab5b9f88256516f2ae618c16a62e`` is from multiple accessions for
*Pythium ultimum* or *Pythium ultimum var. ultimum* plus
`HQ643127.1 <https://www.ncbi.nlm.nih.gov/nucleotide/HQ643127.1>`_ labelled
as *Brevilegnia gracilis*. Again, currently listed as a separate species in
a separate family.

This kind of fluidity in the taxonomy makes it key to document which version
of the taxonomy you are using (month and year for the NCBI taxonomy).
