.. _custom_database_primers:

Specifying custom primers
=========================

Running prepare-reads step
--------------------------

We first ran the pipeline command with :ref:`default settings
<custom_database_defaults>`, but now we will change the primer settings.
Using the actual right primer will extend the *Phytophthora* FASTA sequences
about 60bp (and accept many more non-*Phytophthora*), while treating the
conserved 32bp fragment ``TTTCCGTAGGTGAACCTGCGGAAGGATCATTA`` as if it were
part of the left primer will trim the start of the sequenes.

The up-shot is by cropping about 32bp off the start, and adding about 60bp
at the end, we will no longer get any matches against the default database
with the default classifier (it is too strict, the matches are too distant).
This means before we can run the entire pipeline, we will need to build a
custom database. But first, we will examine the new FASTA intermediate files.

Again we assume you have setup the FASTQ files in ``raw_data/``, and just
run the prepare-reads step:

.. code:: console

    $ mkdir intermediate/
    $ thapbi_pict prepare-reads \
      -i raw_data/ -o intermediate/ \
      --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA \
      --right AGCGTTCTTCATCGATGTGC
    ...
    $ ls -1 intermediate/SRR*.fasta | wc -l
    384

Here we said the left primer is ``GAAGGTGAAGTCGTAACAAGG`` (same as the THAPBI
PICT default) plus ``TTTCCGTAGGTGAACCTGCGGAAGGATCATTA`` (conserved 32bp
region), and that the right primer is ``AGCGTTCTTCATCGATGTGC``. This has
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
a reservoir. Here is with the default primer trimming:

.. code:: console

    $ cat intermediate_defaults/SRR6303586.fasta
    #left_primer:GAAGGTGAAGTCGTAACAAGG
    #right_primer:GCARRGACTTTCGTCCCYRC
    #raw_fastq:70396
    #trimmomatic:70357
    #flash:67831
    #cutadapt:67111
    #abundance:46137
    #threshold:100
    >11b74237bf44899f24ace62be657172a_35127
    TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATT
    TGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAG
    CGTTTGGGCTTCGGCTCGAACTAGTAGCTATCAATTTTAAACCCTTTCTTTAAATACTGAACATACT
    >3cd14145747f26f7461cbba643388ff6_10275
    TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATT
    TGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAG
    CGTTTGGGCTTCGGCTCGAACTAGTAGCTATCAATTTTAAACTCTTTCTTTAAATACTGAACATACT
    >9ba3657d5ecb7d8f35ce7f5feacbc3f7_580
    TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATTT
    GGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAGC
    GTTTGGGCTTCGGCTCGAACTAGTAGCTATCAATTTTAAACCCTTTCTTTAAATACTGAACATACT
    >5bddf030b013e783a218734230c6f542_155
    TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATTT
    GGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAGC
    GTTTGGGCTTCGGCTCGAACTAGTAGCTATCAATTTTAAACTCTTTCTTTAAATACTGAACATACT

Four very similar sequences (differing in the length of the poly-A run, seven
is more common than six, and a ``C/T`` SNP towards the end), all matched to
*Phytophthora chlamydospora* with THAPBI PICT's default settings.

With the new primer setting, which you can see listed at the start of the
header, we again get four sequences passing the abundance threshold:

.. code:: console

    $ cat intermediate/SRR6303586.fasta
    #left_primer:GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA
    #right_primer:AGCGTTCTTCATCGATGTGC
    #raw_fastq:70396
    #trimmomatic:70357
    #flash:67831
    #cutadapt:67341
    #abundance:43923
    #threshold:100
    >e804f4fa9e197115c1f72b943e443dc7_33489
    CCACACCTAAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCC
    TGTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAACTAGTAGCTATC
    AATTTTAAACCCTTTCTTTAAATACTGAACATACTGTGGGGACGAAAGTCTCTGCTTTTAACTAGATAGCAACTTTCAGC
    AGTGGATGTCTAGGCTC
    >3804bc12d180cc145776cc3e77d50561_9746
    CCACACCTAAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCC
    TGTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAACTAGTAGCTATC
    AATTTTAAACTCTTTCTTTAAATACTGAACATACTGTGGGGACGAAAGTCTCTGCTTTTAACTAGATAGCAACTTTCAGC
    AGTGGATGTCTAGGCTC
    >0c2bc6a79b05e53d63636509e9ea8aba_545
    CCACACCTAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCT
    GTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAACTAGTAGCTATCA
    ATTTTAAACCCTTTCTTTAAATACTGAACATACTGTGGGGACGAAAGTCTCTGCTTTTAACTAGATAGCAACTTTCAGCA
    GTGGATGTCTAGGCTC
    >a09c91f2a4813209b3d22847e0b18482_143
    CCACACCTAAAAAACTTTCCACGTGAACCGTATCAACCCCTTAAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCT
    GTAATGGGTCGGCGTGCTGCTGCTGGGCAGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAACTAGTAGCTATCA
    ATTTTAAACTCTTTCTTTAAATACTGAACATACTGTGGGGACGAAAGTCTCTGCTTTTAACTAGATAGCAACTTTCAGCA
    GTGGATGTCTAGGCTC


Again four very similar sequences, each as before but with the starting
``TTTCCGTAGGTGAACCTGCGGAAGGATCATTA`` removed, and instead extended by
``GTGGGGACGAAAGTCTCTGCTTTTAACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC``.

The abundances are similar but slightly lower - there would have been
some minor variants timmed regions which would have been pooled, so with
less trimming we tend to get lower counts.

You can verify by NCBI BLAST online that the first and third (the
``C`` form) give perfect full length matches to published *Phytophthora
chlamydospora*, while an exact match to the ``T`` forms has not been
published at the time of writing (yet they occurs at good abundance in
many of these samples).

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

    $ grep -v "^#" intermediate_defaults/SRR6303596.fasta
    >3dd3b5989ee07ed2d2b3fac826dbb94f_954
    TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAATCTTTCCACGTGAATTGTTTTGCTGTACCTTTGGG
    CTTCGCCGTTGTCTTGTTCTTTTGTAAGAGAAAGGGGGAGGCGCGGTTGGAGGCCATCAGGGGTGTGTTCGTCGCGGTTT
    GTTTCTTTTGTTGGAACTTGCGCGCGGATGCGTCCTTTTGTCAACCCATTTTTTGAATGAAAAACTGATCATACT

There was a single sequence, with no matches (NCBI BLAST suggests this is
*Phytopythium litorale*). Now with the revised primer settings this sequence
is still present but only the second most abundant sequence:

.. code:: console

    $ grep -v "^#" intermediate/SRR6303596.fasta
    >23710597e30e5d95f1d94d6fe8848fb7_40569
    CCACACCAAAAAAACTTTCCACGTGAACCGTTGTAACTATGTTCTGTGCTCTCTTCTCGGAGAGAGCTGAACGAAGGTGG
    GCTGCTTAATTGTAGTCTGCCGATGTACTTTTAAACCCATTAAACTAATACTGAACTATACTCCGAAAACGAAAGTCTTT
    GGTTTTAATCAATAACAACTTTCAGCAGTGGATGTCTAGGCTC
    >b87b957d70d3681d0682103b0052c16b_882
    CCACACCTAAAAATCTTTCCACGTGAATTGTTTTGCTGTACCTTTGGGCTTCGCCGTTGTCTTGTTCTTTTGTAAGAGAA
    AGGGGGAGGCGCGGTTGGAGGCCATCAGGGGTGTGTTCGTCGCGGTTTGTTTCTTTTGTTGGAACTTGCGCGCGGATGCG
    TCCTTTTGTCAACCCATTTTTTGAATGAAAAACTGATCATACTGTGGGGACGAAAGTCTCTGCTTTTAACTAGATAGCAA
    CTTTCAGCAGTGGATGTCTAGGCTC
    >4007e1e8dedb33b5a3c5bc2cfe67c038_391
    CCACACCAAAAAACTTTCCACGTGAACCGTTGTAACTATGTTCTGTGCTCTCTTCTCGGAGAGAGCTGAACGAAGGTGGG
    CTGCTTAATTGTAGTCTGCCGATGTACTTTTAAACCCATTAAACTAATACTGAACTATACTCCGAAAACGAAAGTCTTTG
    GTTTTAATCAATAACAACTTTCAGCAGTGGATGTCTAGGCTC
    >f2a354f8c74387a287be1d08f31df143_190
    CCACACCAAAAAAACTTTCCACGTGAACCGTTGTAACTATGTTCTGTGCTCTCTTCTCGGAGAGAGCTGAACGAAGGTGG
    GCTGCTTAATTGTAGTCTGCCGATGTACTTTTAAACCCATTAAACTAATACTGAACTATACTCCGGAAACGAAAGTCTTT
    GGTTTTAATCAATAACAACTTTCAGCAGTGGATGTCTAGGCTC
    >818d3263599c4929cf2ef4a33f952949_128
    CCACACCAAAAAAACTTTCCACGTGAACCGTTGTAACTATGTTCTGTGCTCTCTTCTCGGAGAGAGCTGAACGAAGGTGG
    GCTGCTTAATTGTAGTCTGCCGATGTACTTTTAAACCCATTAAACTAATACTGAACTATACTCCGAAAACGAAAGTCTTT
    GGTTTTAATCAATAACAACTTTCAGCAGTGGATGTCTAGGCGC
    >48bcfd8e8daaa8351cb24b7deb63a4bc_102
    CCACACCAAAAAAACTTTCCACGTGAACCGTTGTAACTATGTTCTGTGCTCTCTTCTCGGAGAGAGCTGAACGAAGGTGG
    GCTGCTTAATTGTAGTCTGCCGATGTACTTTTAAACCCATTAAACTAATACTGAACTATACTCCGAAAACGAAAGTCTTT
    GGTTTTAATCAATAACAACTTTCAGCAGTGGATGTCTAGGCCC

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

    $ grep -v "^#" intermediate_defaults/SRR6303948.fasta
    >dcd6316eb77be50ee344fbeca6e005c7_1437
    TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAACCCTTAGTTG
    GGGGCTTCTGTTCGGCTGGCTTCGGCTGGCTGGGCGGCGGCTCTATCATGGCGAGCGCTTGAGCCTTCGGGTCTGAGCTA
    GTAGCCCACTTTTTAAACCCATTCCTAAATACTGAATATACT

Now with the revised primer settings, we get a further five sequences - and
the extended *Phytophthora ramorum* sequence drops to second most abundant:

.. code:: console

    $ grep -v "^#" intermediate/SRR6303948.fasta
    >f2d4b17eb421d8c52320c2bd883e77eb_5322
    CCACACCAAAAAAACACCCCACGTGAATTGTACTGTATGAGCTATGTGCTGCGGATTTCTGCGGCTTAGCGAAGGTTTCG
    AAAGAGACCGATGTACTTTTAAACCCCTTTACATTACTGTCTGATAAATTACATTGCAAACATTTAAAGTGGTTGCTCTT
    AATTTAACATACAACTTTCAACAGTGGATGTCTAGGCTC
    >283ae6bd5fb4ba9ec5fba94a3f67b43d_1315
    CCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAACCCTTAGTTGGGGGCTTCTGTTCGGCTGGCTTCGGCTGGCTG
    GGCGGCGGCTCTATCATGGCGAGCGCTTGAGCCTTCGGGTCTGAGCTAGTAGCCCACTTTTTAAACCCATTCCTAAATAC
    TGAATATACTGTGGGGACGAAAGTCTCTGCTTTTAACTAGATAGCAACTTTCAGCAGTGGATGTCTAGGCTC
    >9242fdd7b824838e583731161706caf1_437
    CCACACCAAAAAAACTTACCACGTGAATCTGTACTGTTTAGTTTTGTGCTGCGTTCGAAAGGATGCGGCTAAACGAAGGT
    TGGCTTGATTACTTCGGTAATTAGGCTGGCTGATGTACTCTTTTAAACCCCTTCATACCAAAATACTGATTTATACTGTG
    AGAATGAAAATTCTTGCTTTTAACTAGATAACAACTTTCAACAGTGGATGTCTAGGCTC
    >5d245b9970ea98e368afdd370a3dfae6_229
    ATCTATCACAATCCACACCTGTGAACTTGCTTGTTGGCCTCTGCATGTGCTTCGGTATGTGCAGGTTGAGCCGATCGGAT
    TAACTTCTGGTCGGCTTGGGGCCTCAACCCAATCCTCGGATTGGTTTGGGGTCGGTCTCTATTAACAACCAACACCAAAC
    CAAACTATAAAAAAACTGAGAATGGCTTAGAGCCAAACTCACTAACCAAGACAACTCTGAACAACGGATATCTTGGCTA
    >31bac939435fe6972e3e2d004937c876_189
    CCACACCTAAAAACTTTCCACGTGAATCGTTCTATATAGCTTTGTGCTTTGCGGAAACGCGAGGCTAAGCGAAGGATTAG
    CAAAGTAGTACTTCGGTGCGAAACACTTTTCCGATGTATTTTTCAAACCCTTTTACTTATACTGAACTATACTCTAAGAC
    GAAAGTCTTGGTTTTAATCCACAACAACTTTCAGCAGTGGATGTCTAGGCTC
    >44b31ed4182973c57683a561485745c4_103
    CCACACCAAAAAACACCCCACGTGAATTGTACTGTATGAGCTATGTGCTGCGGATTTCTGCGGCTTAGCGAAGGTTTCGA
    AAGAGACCGATGTACTTTTAAACCCCTTTACATTACTGTCTGATAAATTACATTGCAAACATTTAAAGTGGTTGCTCTTA
    ATTTAACATACAACTTTCAACAGTGGATGTCTAGGCTC

NCBI BLAST suggests the new sequences could all be *Oomycetes*, but there
are no very close matches - and some of the tenous best matches include
uncultured fungus, green algae, and even green plants.
