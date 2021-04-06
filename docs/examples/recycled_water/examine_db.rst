.. _custom_database_examine:

Examining the database
======================

This example follows on from :ref:`custom_database_building`, and assumes
you have used ``thapbi_pict curated-import`` with the provided FASTA file
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

This gives a short table, summarised as follows:

========== ============ ======= ====== ================================ ======== ========
Identifier Genus        Species TaxID  ITS1-MD5                         ITS1-seq Sequence
---------- ------------ ------- ------ -------------------------------- -------- --------
DQ297398.1 Phytophthora fallax  360399 693cf88b7f57bcc7a3532a6b7ff0268a CCA...   CCA...
HQ261557.1 Phytophthora fallax  360399 693cf88b7f57bcc7a3532a6b7ff0268a CCA...   CCA...
HQ261558.1 Phytophthora fallax  360399 693cf88b7f57bcc7a3532a6b7ff0268a CCA...   CCA...
HQ261559.1 Phytophthora fallax  360399 693cf88b7f57bcc7a3532a6b7ff0268a CCA...   CCA...
DQ297392.1 Phytophthora fallax  360399 da7ff4ae11bdb6cc2b8c2aea3937481f CCA...   CCA...
========== ============ ======= ====== ================================ ======== ========

The final two columns give the marker sequence, and the full sequence it was
extracted from using any primer trimming -- in this case, they match as the
FASTA file used was pre-trimmed. The (shorter) marker sequence gives the MD5
checksum column.

Adding ``-m`` or ``--minimal`` to the command gives instead:

.. code:: console

    $ thapbi_pict dump -d Redekar_et_al_2019_sup_table_3.sqlite \
      -g Phytophthora -s fallax -o P_fallax.tsv -m
    Wrote 2 txt format entries to 'P_fallax.tsv'

Now the table only has one data row per unique marker sequences:

================================ =================== ========
MD5                              Species             Sequence
-------------------------------- ------------------- --------
693cf88b7f57bcc7a3532a6b7ff0268a Phytophthora fallax CCA...
da7ff4ae11bdb6cc2b8c2aea3937481f Phytophthora fallax CCA...
================================ =================== ========

Instead we can ask for FASTA output:

.. code:: console

    $ thapbi_pict dump -d Redekar_et_al_2019_sup_table_3.sqlite \
      -g Phytophthora -s fallax -f fasta -o P_fallax.fasta
    Wrote 2 fasta format entries to 'P_fallax.fasta'

This produces a short FASTA file as follows (with line wrapping added
for display)::

.. code:: console

    $ cat P_fallax.fasta
    >DQ297398.1 Phytophthora fallax;HQ261557.1 Phytophthora fallax;HQ261558.1
    Phytophthora fallax;HQ261559.1 Phytophthora fallax
    CCACACCTAAAAAAATTCCACGTGAACTGTATTGTCAACCAAATTCGGGGATTCCTTGCTAGCGTGCCTTCGGGCGTGCC
    GGTAGGTTGAGACCCATCAAACGAAAACATCGGCTGAAAGGTCGGAGCCAGTAGTTACCTTTGTAAACCCTTTACTAAAT
    ACTGAAAAACTGTGGGGACGAAAGTCCTTGCTTTTACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC
    >DQ297392.1 Phytophthora fallax
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
space, semi-colon separated list of taxonomic assignments, new line, sequences,
new line. Again, there is deliberately no sequence line wrapping in the file
itself.

Edit graph
----------

In the worked example with the default database, we introduced the
``edit-graph`` command for use with CytoScape to examine the sequence space of
the samples. It can also be run on a database alone provided you include the
``-s`` or ``--showdb`` switch:

.. code:: console

    $ thapbi_pict edit-graph -s \
      -d Redekar_et_al_2019_sup_table_3.sqlite \
      -o Redekar_et_al_2019_sup_table_3.xgmml
    Loaded 838 unique sequences from database
    Computed 350703 Levenshtein edit distances between 838 sequences.
    Will draw 533 nodes with at least one edge (305 are isolated sequences).

Of the 838 unique sequences in the database, just over three hundred are
isolated sequences (over 3bp edits away from anything else). The remaining
five hundred plus give us an interesting edit distance graph.

Opening this in CytoScape the first thing that struck me was the largest two
components are both for *Pythium regulare* - suggesting if these are truely
all from one species that it has at least two distinct ITS1 markers in the
genome?

Another use of this view would be to consider the genus conflicts reported
by the ``thapbi_pict conflicts`` command - most of the handful of *Lagenidium*
and *Brevilegnia* nodes are isolated.
