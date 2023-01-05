.. _custom_database_primers:

Specifying custom primers
=========================

Running prepare-reads step
--------------------------

We first ran the pipeline command with :ref:`default settings
<custom_database_defaults>`, if you skipped that we can do just the reads now:

.. code:: console

    $ mkdir -p intermediate_defaults/
    $ thapbi_pict prepare-reads -i raw_data/ -o intermediate_defaults/
    ...
    $ ls -1 intermediate_defaults/ITS1/SRR*.fasta | wc -l
    384

We then created a database from the Redekar *et al.* (2019) reference
accessions with their primers. Now we can run the pipeline again with this,
which will start by applying the prepare-reads step to the FASTQ files in
``raw_data/``:

.. code:: console

    $ mkdir -p intermediate/
    $ thapbi_pict prepare-reads -i raw_data/ -o intermediate_long/ \
      --db Redekar_et_al_2019_sup_table_3.sqlite
    ...
    $ ls -1 intermediate_long/ITS1-long/SRR*.fasta | wc -l
    384

Here the database says the left primer is ``GAAGGTGAAGTCGTAACAAGG`` (same as
the THAPBI PICT default) plus ``TTTCCGTAGGTGAACCTGCGGAAGGATCATTA`` (conserved
32bp region), and that the right primer is ``AGCGTTCTTCATCGATGTGC``. This has
reverse complement ``GCACATCGATGAAGAACGCT`` and is found about 60bp downstream
of the default right primer in *Phytophthora*, and should also match *Pythium*
and *Phytopythium* species.

i.e. We should now find the *Phytophthora* FASTA sequences extracted are about
60 - 32 = 28bp longer, and many more non-*Phytophthora* are accepted.

Will now pick a couple of samples to compare and contrast with the first run.
For clarity these examples are deliberately from the less diverse samples.
The FASTA sequences have been line wrapped at 80bp for display.

Longer sequences
----------------

We will start with ``SRR6303586`` aka ``OSU483``, a leaf-baiting sample from
a reservoir. With the default primer trimming looking at the reads report, or
the simpler sally table, focusing on just the one sample and filtering out
non-zero counts:

.. code:: console

    $ tail -n +10 summary/recycled-water-defaults.ITS1.tally.tsv \
      | cut -f 3,386 | grep -v "^0"
    <SEE TABLE BELOW>

You could instead select and filter on this column in Excel:

========== ===================================================================================================================================================================================================================================
SRR6303586 Sequence
========== ===================================================================================================================================================================================================================================
35109      TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAACTAGTAGCTATCAATTTTAAACCCTTTCTTTAAATACTGAACATACT
10271      TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAACTAGTAGCTATCAATTTTAAACTCTTTCTTTAAATACTGAACATACT
580        TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAACTAGTAGCTATCAATTTTAAACCCTTTCTTTAAATACTGAACATACT
157        TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAACTAGTAGCTATCAATTTTAAACTCTTTCTTTAAATACTGAACATACT
========== ===================================================================================================================================================================================================================================

Four very similar sequences (differing in the length of the poly-A run, seven
is more common than six, and a ``C/T`` SNP towards the end), all matched to
*Phytophthora chlamydospora* with THAPBI PICT's default settings.

With the new primer setting, which you can see listed at the start of the
header, we again get four sequences passing the abundance threshold:

.. code:: console

    $ tail -n +10 summary/recycled-water-custom.ITS1-long.tally.tsv \
      | cut -f 3,386 | grep -v "^0"
    <SEE TABLE BELOW>

As before, you may prefer to open this as a spreadsheet:

========== =================================================================================================================================================================================================================================================================
SRR6303586 Sequence
========== =================================================================================================================================================================================================================================================================
33451      CCACACCTAAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAACTAGTAGCTATCAATTTTAAACCCTTTCTTTAAATACTGAACATACTGTGGGGACGAAAGTCTCTGCTTTTAACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC
9729       CCACACCTAAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAACTAGTAGCTATCAATTTTAAACTCTTTCTTTAAATACTGAACATACTGTGGGGACGAAAGTCTCTGCTTTTAACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC
545        CCACACCTAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAACTAGTAGCTATCAATTTTAAACCCTTTCTTTAAATACTGAACATACTGTGGGGACGAAAGTCTCTGCTTTTAACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC
143        CCACACCTAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAACTAGTAGCTATCAATTTTAAACTCTTTCTTTAAATACTGAACATACTGTGGGGACGAAAGTCTCTGCTTTTAACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC
========== =================================================================================================================================================================================================================================================================


Again four very similar sequences, each as before but with the starting
``TTTCCGTAGGTGAACCTGCGGAAGGATCATTA`` removed, and instead extended by
``GTGGGGACGAAAGTCTCTGCTTTTAACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC``.

The abundances are similar but slightly lower - there would have been some
minor variation in trimmed regions which would have been pooled, so with
less trimming we tend to get lower counts.

You can verify by NCBI BLAST online that the first and third (the
``C`` form) give perfect full length matches to published *Phytophthora
chlamydospora*, while an exact match to the ``T`` forms has not been
published at the time of writing (yet this occurs at good abundance in many of
these samples).

Losing sequences
----------------

If you examine ``SRR6303588`` you will see a similar example,
starting with five unique sequences (with one only just above the
default abundance threshold), dropping to four unique sequences.

Finding *Pythium*
-----------------

Now for a more interesting example, ``SRR6303596`` aka ``OSU121``, another
leaf baiting sample but from runoff water. With the defaults (using ``grep``
to omit the header):

.. code:: console

    $ tail -n +10 summary/recycled-water-defaults.ITS1.tally.tsv \
      | cut -f 13,386 | grep -v "^0"
    <SEE TABLE BELOW>

As a table,

========== ===========================================================================================================================================================================================================================================
SRR6303596 Sequence
========== ===========================================================================================================================================================================================================================================
953        TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAATCTTTCCACGTGAATTGTTTTGCTGTACCTTTGGGCTTCGCCGTTGTCTTGTTCTTTTGTAAGAGAAAGGGGGAGGCGCGGTTGGAGGCCATCAGGGGTGTGTTCGTCGCGGTTTGTTTCTTTTGTTGGAACTTGCGCGCGGATGCGTCCTTTTGTCAACCCATTTTTTGAATGAAAAACTGATCATACT
========== ===========================================================================================================================================================================================================================================

There was a single sequence, with no matches (NCBI BLAST suggests this is
*Phytopythium litorale*). Now with the revised primer settings this sequence
is still present but only the second most abundant sequence:

.. code:: console

    $ tail -n +10 summary/recycled-water-custom.ITS1-long.tally.tsv \
      | cut -f 13,386 | grep -v "^0"
    <SEE TABLE BELOW>

As a table, note this is sorted by global abundance:

========== =========================================================================================================================================================================================================================================================================
SRR6303596 Sequence
========== =========================================================================================================================================================================================================================================================================
40503      CCACACCAAAAAAACTTTCCACGTGAACCGTTGTAACTATGTTCTGTGCTCTCTTCTCGGAGAGAGCTGAACGAAGGTGGGCTGCTTAATTGTAGTCTGCCGATGTACTTTTAAACCCATTAAACTAATACTGAACTATACTCCGAAAACGAAAGTCTTTGGTTTTAATCAATAACAACTTTCAGCAGTGGATGTCTAGGCTC
878        CCACACCTAAAAATCTTTCCACGTGAATTGTTTTGCTGTACCTTTGGGCTTCGCCGTTGTCTTGTTCTTTTGTAAGAGAAAGGGGGAGGCGCGGTTGGAGGCCATCAGGGGTGTGTTCGTCGCGGTTTGTTTCTTTTGTTGGAACTTGCGCGCGGATGCGTCCTTTTGTCAACCCATTTTTTGAATGAAAAACTGATCATACTGTGGGGACGAAAGTCTCTGCTTTTAACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC
388        CCACACCAAAAAACTTTCCACGTGAACCGTTGTAACTATGTTCTGTGCTCTCTTCTCGGAGAGAGCTGAACGAAGGTGGGCTGCTTAATTGTAGTCTGCCGATGTACTTTTAAACCCATTAAACTAATACTGAACTATACTCCGAAAACGAAAGTCTTTGGTTTTAATCAATAACAACTTTCAGCAGTGGATGTCTAGGCTC
128        CCACACCAAAAAAACTTTCCACGTGAACCGTTGTAACTATGTTCTGTGCTCTCTTCTCGGAGAGAGCTGAACGAAGGTGGGCTGCTTAATTGTAGTCTGCCGATGTACTTTTAAACCCATTAAACTAATACTGAACTATACTCCGAAAACGAAAGTCTTTGGTTTTAATCAATAACAACTTTCAGCAGTGGATGTCTAGGCGC
102        CCACACCAAAAAAACTTTCCACGTGAACCGTTGTAACTATGTTCTGTGCTCTCTTCTCGGAGAGAGCTGAACGAAGGTGGGCTGCTTAATTGTAGTCTGCCGATGTACTTTTAAACCCATTAAACTAATACTGAACTATACTCCGAAAACGAAAGTCTTTGGTTTTAATCAATAACAACTTTCAGCAGTGGATGTCTAGGCCC
190        CCACACCAAAAAAACTTTCCACGTGAACCGTTGTAACTATGTTCTGTGCTCTCTTCTCGGAGAGAGCTGAACGAAGGTGGGCTGCTTAATTGTAGTCTGCCGATGTACTTTTAAACCCATTAAACTAATACTGAACTATACTCCGGAAACGAAAGTCTTTGGTTTTAATCAATAACAACTTTCAGCAGTGGATGTCTAGGCTC
========== =========================================================================================================================================================================================================================================================================

The probable *Phytopythium litorale* has been joined by five shorter
and very similar sequences (differing by a handful of SNPs and a
poly-A length change), which NCBI BLAST matches suggest are all
*Pythium coloratum/dissotocum*.

Finding more
------------

Another interesting example, ``SRR6303948`` aka ``OSU536.s203``,
from a runoff filtration sample. First with the default settings,
a single unique sequence matching *Phytophthora ramorum*:

.. code:: console

    $ tail -n +10 summary/recycled-water-defaults.ITS1.tally.tsv \
      | cut -f 365,386 | grep -v "^0"
    <SEE TABLE BELOW>

As a table,

========== ==========================================================================================================================================================================================================
SRR6303948 Sequence
========== ==========================================================================================================================================================================================================
1439       TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAACCCTTAGTTGGGGGCTTCTGTTCGGCTGGCTTCGGCTGGCTGGGCGGCGGCTCTATCATGGCGAGCGCTTGAGCCTTCGGGTCTGAGCTAGTAGCCCACTTTTTAAACCCATTCCTAAATACTGAATATACT
========== ==========================================================================================================================================================================================================

Now with the revised primer settings, we get a further nine sequences - and
the extended *Phytophthora ramorum* sequence drops to third most abundant:

.. code:: console

    $ tail -n +10 summary/recycled-water-custom.ITS1-long.tally.tsv \
      | cut -f 365,386 | grep -v "^0"
    <SEE TABLE BELOW>

As a table, note this is sorted by global abundance:

========== ================================================================================================================================================================================================================================================
SRR6303948 Sequence
========== ================================================================================================================================================================================================================================================
3287       CCACACCCGGGATCCTCGATCTTTCTCCTAGGTTAATTGTTGGGCCCTTTGAGGGTGGGCCTTAGGTGCGCTCAAGGATTTTTTCCTGTCCCATGTAGCTTTACTTATTTTTTTGCCTGGGTAAATGATGGATTATTTTTACAACTTTCAGCAATGGATGTCTAGGCTC
438        CCACACCAAAAAAACTTACCACGTGAATCTGTACTGTTTAGTTTTGTGCTGCGTTCGAAAGGATGCGGCTAAACGAAGGTTGGCTTGATTACTTCGGTAATTAGGCTGGCTGATGTACTCTTTTAAACCCCTTCATACCAAAATACTGATTTATACTGTGAGAATGAAAATTCTTGCTTTTAACTAGATAACAACTTTCAACAGTGGATGTCTAGGCTC
5329       CCACACCAAAAAAACACCCCACGTGAATTGTACTGTATGAGCTATGTGCTGCGGATTTCTGCGGCTTAGCGAAGGTTTCGAAAGAGACCGATGTACTTTTAAACCCCTTTACATTACTGTCTGATAAATTACATTGCAAACATTTAAAGTGGTTGCTCTTAATTTAACATACAACTTTCAACAGTGGATGTCTAGGCTC
144        CCACACCCGGGATCCTCGATCTTTCTCCTAGGTTAATTATTGGGCCCTTTGAGGGTGGGCCTTAGGTGCGCTCAAGGATTTTTTCCTGTCCCATGTAGCTTTACTTATTTTTTTGCCTGGGTAAATGATGGATTATTTTTACAACTTTCAGCAATGGATGTCTAGGCTC
230        AATCTATCACAATCCACACCTGTGAACTTGCTTGTTGGCCTCTGCATGTGCTTCGGTATGTGCAGGTTGAGCCGATCGGATTAACTTCTGGTCGGCTTGGGGCCTCAACCCAATCCTCGGATTGGTTTGGGGTCGGTCTCTATTAACAACCAACACCAAACCAAACTATAAAAAAACTGAGAATGGCTTAGAGCCAAACTCACTAACCAAGACAACTCTGAACAACGGATATCTTGGCTA
1319       CCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAACCCTTAGTTGGGGGCTTCTGTTCGGCTGGCTTCGGCTGGCTGGGCGGCGGCTCTATCATGGCGAGCGCTTGAGCCTTCGGGTCTGAGCTAGTAGCCCACTTTTTAAACCCATTCCTAAATACTGAATATACTGTGGGGACGAAAGTCTCTGCTTTTAACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC
224        CCACACCCGGGATCCTCGATCTTTCTCCTAGGTTAATTGTTTGGCCCTTTGAGGGTGGGCCTTAGGTGCGCTCAAGGATTTTTTCCTGTCCCATGTAGCTTTACTTATTTTTTTGCCTGGGTAAATGATGGATTATTTTTACAACTTTCAGCAATGGATGTCTAGGCTC
231        CCACACCCGGGATCCTCGATCTTTCTCCTAGGTTAATTGTTGGGCCCTTTGAGGGTGGGCCTTAGGTGCGCTCAAGGATTTTTTCCTGTCCCATGTAGCTTTACTTATTTTTTTGCCTGGGTAAATGATGGATTATTTTTACAACTTTCAGCAACGGATGTCTAGGCTC
102        CCACACCAAAAAACACCCCACGTGAATTGTACTGTATGAGCTATGTGCTGCGGATTTCTGCGGCTTAGCGAAGGTTTCGAAAGAGACCGATGTACTTTTAAACCCCTTTACATTACTGTCTGATAAATTACATTGCAAACATTTAAAGTGGTTGCTCTTAATTTAACATACAACTTTCAACAGTGGATGTCTAGGCTC
189        CCACACCTAAAAACTTTCCACGTGAATCGTTCTATATAGCTTTGTGCTTTGCGGAAACGCGAGGCTAAGCGAAGGATTAGCAAAGTAGTACTTCGGTGCGAAACACTTTTCCGATGTATTTTTCAAACCCTTTTACTTATACTGAACTATACTCTAAGACGAAAGTCTTGGTTTTAATCCACAACAACTTTCAGCAGTGGATGTCTAGGCTC
========== ================================================================================================================================================================================================================================================

NCBI BLAST suggests some of the new sequences could be *Oomycetes*, but there
are no very close matches - and some of the tenuous best matches include
uncultured fungus, diatoms, green algae, and even green plants.
