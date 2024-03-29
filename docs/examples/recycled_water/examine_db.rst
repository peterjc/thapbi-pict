.. _custom_database_examine:

Examining the database
======================

This example follows on from :ref:`custom_database_building`, and assumes
you have used ``thapbi_pict import`` with the provided FASTA file
(based on Supplementary Table 3 in Redekar *et al.* 2019), and created a
THAPBI PICT database named ``Redekar_et_al_2019_sup_table_3.sqlite``.

As the extension might suggest, this is an Sqlite v3 database, and can be
examined directly at the command line if you are very curious. However,
we will briefly review the provided commands within THAPBI PICT for checking
a database.

Database export
---------------

The ``thapbi_pict dump`` command is intended for database export and/or
answering simple queries without needing to use SQL to query the database.
It defaults to giving plain text tab separated tables, but FASTA is also
supported:

.. code:: console

    $ thapbi_pict dump -h
    ...

By default it outputs all the sequences, but you can do simple taxonomic
filtering at genus or species level, for example:

.. code:: console

    $ thapbi_pict dump -d Redekar_et_al_2019_sup_table_3.sqlite \
      -g Phytophthora -s fallax -o P_fallax.tsv
    Wrote 5 txt format entries to 'P_fallax.tsv'
    $ cut -c 1-84 P_fallax.tsv
    <SEE TABLE BELOW>

This gives a short table, with the sequence truncated for display:

========= ========== ============ ======= ====== ================================ ========
#Marker   Identifier Genus        Species TaxID  MD5                              Sequence
========= ========== ============ ======= ====== ================================ ========
ITS1-long DQ297398.1 Phytophthora fallax  360399 693cf88b7f57bcc7a3532a6b7ff0268a CCA
ITS1-long HQ261557.1 Phytophthora fallax  360399 693cf88b7f57bcc7a3532a6b7ff0268a CCA
ITS1-long HQ261558.1 Phytophthora fallax  360399 693cf88b7f57bcc7a3532a6b7ff0268a CCA
ITS1-long HQ261559.1 Phytophthora fallax  360399 693cf88b7f57bcc7a3532a6b7ff0268a CCA
ITS1-long DQ297392.1 Phytophthora fallax  360399 da7ff4ae11bdb6cc2b8c2aea3937481f CCA
========= ========== ============ ======= ====== ================================ ========

The final columns give the amplicon marker sequence and its MD5 checksum.

Adding ``-m`` or ``--minimal`` to the command gives instead:

.. code:: console

    $ thapbi_pict dump -d Redekar_et_al_2019_sup_table_3.sqlite \
      -g Phytophthora -s fallax -o P_fallax.tsv -m
    Wrote 2 txt format entries to 'P_fallax.tsv'
    $ cut -c 1-56 P_fallax.tsv
    <SEE TABLE BELOW>

Now the table only has one data row per unique marker sequence, again showing
this with the sequence truncated:

================================ =================== ========
#MD5                             Species             Sequence
================================ =================== ========
693cf88b7f57bcc7a3532a6b7ff0268a Phytophthora fallax CCA
da7ff4ae11bdb6cc2b8c2aea3937481f Phytophthora fallax CCA
================================ =================== ========

Alternatively, we can ask for FASTA output:

.. code:: console

    $ thapbi_pict dump -d Redekar_et_al_2019_sup_table_3.sqlite \
      -g Phytophthora -s fallax -f fasta -o P_fallax.fasta
    Wrote 2 fasta format entries to 'P_fallax.fasta'

This produces a short FASTA file as follows (with line wrapping added
for display):

.. code:: console

    $ cat P_fallax.fasta
    >DQ297398.1 Phytophthora fallax taxid=360399;HQ261557.1 Phytophthora fallax
    taxid=360399;HQ261558.1 Phytophthora fallax taxid=360399;HQ261559.1 Phytophthora
    fallax taxid=360399
    CCACACCTAAAAAAATTCCACGTGAACTGTATTGTCAACCAAATTCGGGGATTCCTTGCTAGCGTGCCTTCGGGCGTGCC
    GGTAGGTTGAGACCCATCAAACGAAAACATCGGCTGAAAGGTCGGAGCCAGTAGTTACCTTTGTAAACCCTTTACTAAAT
    ACTGAAAAACTGTGGGGACGAAAGTCCTTGCTTTTACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC
    >DQ297392.1 Phytophthora fallax taxid=360399
    CCACACCTTAAAAAATTCCACGTGAACTGTATTGTCAACCAAATTCGGGGATTCCTTGCTAGCGTGCCTTCGGGCGTGCC
    GGTAGGTTGAGACCCATCAAACGAAAACATCGGCTGAAAGGTCGGAGCCAGTAGTTACCTTTGTAAACCCTTTACTAAAT
    ACTGAAAAACTGTGGGGACGAAAGTCCTTGCTTTTACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC

To be clear, each FASTA record is written as two potentially very long lines.
The first title line consists of the FASTA new record ``>`` marker and then
four semi-colon separated accessions with species. The sequence shared by those
four entries is given on the second line (without line breaks as markers tend
not to be overly long, and it facilitates command line analysis/debugging).

Using the optional ``-m`` or ``--minimal`` switch changes the FASTA output to:

.. code:: console

    $ thapbi_pict dump -d Redekar_et_al_2019_sup_table_3.sqlite \
      -g Phytophthora -s fallax -f fasta -o P_fallax_minimal.fasta -m
    Wrote 2 fasta format entries to 'P_fallax_minimal.fasta'
    $ cat P_fallax_minimal.fasta
    >693cf88b7f57bcc7a3532a6b7ff0268a Phytophthora fallax
    CCACACCTAAAAAAATTCCACGTGAACTGTATTGTCAACCAAATTCGGGGATTCCTTGCTAGCGTGCCTTCGGGCGTGCC
    GGTAGGTTGAGACCCATCAAACGAAAACATCGGCTGAAAGGTCGGAGCCAGTAGTTACCTTTGTAAACCCTTTACTAAAT
    ACTGAAAAACTGTGGGGACGAAAGTCCTTGCTTTTACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC
    >da7ff4ae11bdb6cc2b8c2aea3937481f Phytophthora fallax
    CCACACCTTAAAAAATTCCACGTGAACTGTATTGTCAACCAAATTCGGGGATTCCTTGCTAGCGTGCCTTCGGGCGTGCC
    GGTAGGTTGAGACCCATCAAACGAAAACATCGGCTGAAAGGTCGGAGCCAGTAGTTACCTTTGTAAACCCTTTACTAAAT
    ACTGAAAAACTGTGGGGACGAAAGTCCTTGCTTTTACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC

This discards the original accessions and instead uses ``>``, MD5 checksum,
space, semi-colon separated list of taxonomic assignments, new line, sequence,
new line. Again, there is deliberately no sequence line wrapping in the file
itself.

Edit graph
----------

In the worked example with the default database, we introduced the
``edit-graph`` command for use with CytoScape to examine the sequence space of
the samples. It can also be run on a database alone provided you include the
``-k`` or ``--marker`` switch:

.. code:: console

    $ thapbi_pict edit-graph -k ITS1-long \
      -d Redekar_et_al_2019_sup_table_3.sqlite \
      -o Redekar_et_al_2019_sup_table_3.xgmml
    Loaded 838 unique ITS1-long sequences from DB.
    Computed Levenshtein edit distances.
    Will draw 533 nodes with at least one edge (305 are isolated sequences).

Of the 838 unique sequences in the database, just over three hundred are
isolated sequences (over 3bp edits away from anything else). The remaining
five hundred plus give us an interesting edit distance graph.

Opening this in CytoScape the first thing that struck me was the largest two
components are both for *Pythium regulare* - suggesting if these are truly
all from one species that it has at least two distinct ITS1 markers in the
genome?

Another use of this view would be to consider the genus conflicts reported
by the ``thapbi_pict conflicts`` command - most of the handful of *Lagenidium*
and *Brevilegnia* nodes are isolated.
