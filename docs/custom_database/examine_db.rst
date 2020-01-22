.. _custom_database_examine:

Examining the database
======================

This example follows on from :ref:`custom_database_building`, and assumes
you have used ``thapbi_pict curated-import`` with the provided FASTA file
(based on Supplementary Table 3 in Redekar *et al* 2019), and created a
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
    Wrote 8 txt format entries to 'P_fallax.tsv'

This would give a a short table, summarised as follows:

========== ============	======= ===== ================================ ======== ========
Identifier Genus        Species TaxID ITS1-MD5                         ITS1-seq Sequence
---------- ------------ ------- ----- -------------------------------- -------- --------
DQ297398.1 Phytophthora fallax        693cf88b7f57bcc7a3532a6b7ff0268a CCACA... CCACA...
HQ261557.1 Phytophthora fallax        693cf88b7f57bcc7a3532a6b7ff0268a CCACA... CCACA...
HQ261558.1 Phytophthora fallax        693cf88b7f57bcc7a3532a6b7ff0268a CCACA... CCACA...
HQ261559.1 Phytophthora fallax        693cf88b7f57bcc7a3532a6b7ff0268a CCACA... CCACA...
DQ297392.1 Phytophthora fallax        da7ff4ae11bdb6cc2b8c2aea3937481f CCACA... CCACA...
========== ============ ======= ===== ================================ ======== ========

The final two columns give the marker sequence, and the full sequence it was
extracted from using any primer trimming -- in this case, they match as the
FASTA file used was pre-trimmed. The (shorter) marker sequence gives the MD5
checksum column.

Adding ``-m`` or ``--minimal`` to the command gives instead:

.. code:: console

    $ thapbi_pict dump -d Redekar_et_al_2019_sup_table_3.sqlite \
      -g Phytophthora -s fallax -o P_fallax.tsv -m
    Wrote 1 txt format entries to 'P_fallax.tsv'

Now the table only has one data row per unique marker sequences:

================================ =================== ========
MD5                              Species             Sequence
-------------------------------- ------------------- --------
693cf88b7f57bcc7a3532a6b7ff0268a Phytophthora fallax CCACA...
da7ff4ae11bdb6cc2b8c2aea3937481f Phytophthora fallax CCACA...
================================ =================== ========

Instead we can ask for FASTA output:

.. code:: console

    $ thapbi_pict dump -d Redekar_et_al_2019_sup_table_3.sqlite \
      -g Phytophthora -s fallax -f fasta -o P_fallax.fasta
    Wrote 2 fasta format entries to 'P_fallax.fasta'

This produces a short FASTA file as follows (with line wrapping added
for display)::

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
five semi-colon separated accessions with species. The sequence shared by those
five entries is given on the second line (without line breaks as markers tend
not to be overly long, and it facilitates command line analysis/debugging).

Using the optional ``-m`` or ``--minimal`` switch changes the FASTA output to::

    >693cf88b7f57bcc7a3532a6b7ff0268a Phytophthora fallax
    CCACACCTAAAAAAATTCCACGTGAACTGTATTGTCAACCAAATTCGGGGATTCCTTGCTAGCGTGCCTTCGGGCGTGCC
    GGTAGGTTGAGACCCATCAAACGAAAACATCGGCTGAAAGGTCGGAGCCAGTAGTTACCTTTGTAAACCCTTTACTAAAT
    ACTGAAAAACTGTGGGGACGAAAGTCCTTGCTTTTACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC
    >da7ff4ae11bdb6cc2b8c2aea3937481f Phytophthora fallax
    CCACACCTTAAAAAATTCCACGTGAACTGTATTGTCAACCAAATTCGGGGATTCCTTGCTAGCGTGCCTTCGGGCGTGCC
    GGTAGGTTGAGACCCATCAAACGAAAACATCGGCTGAAAGGTCGGAGCCAGTAGTTACCTTTGTAAACCCTTTACTAAAT
    ACTGAAAAACTGTGGGGACGAAAGTCCTTGCTTTTACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC

This discards the original accessions and instead uses ``>``, MD5 checksum, space,
semi-colon separated list of taxonomic assignments, new line, sequences, new line.
Again, there is deliberatly no sequence line wrapping.
