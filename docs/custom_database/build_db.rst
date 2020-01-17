.. _custom_database_building:

Introduction
============

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

    $ rm -rf Redekar_et_al_2019_sup_table_3_download.fasta
    $ for ACC in `cut -f 1 Redekar_et_al_2019_sup_table_3.tsv`; do cat "cache/$ACC.fasta" >> Redekar_et_al_2019_sup_table_3_download.fasta; done
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

So far so good - all 1454 sequences were loaded.

The ``-x`` (or ``--lax``) turns of species name validation, we'll come back
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

Curated import
--------------

The ``thapbi_pict curated-import`` mentioned above differs from the
``thapbi_pict ncbi-import`` command in two key points. First, it expects
the sequences to be pre-trimmed (it does not look for an remove primers).
Second, it does not use heuristics but simply assumes the FASTA description
line is an identifier followed by the species name *only*.
