Presence and absence
====================

This example is a controlled setup where we know what the classifier
ought ideally to report for every single sample.

We have replicated negative controls (which should have no marker
sequences present), and plenty of positive controls (which should
have the expected 19 species only).

Of course, just as in the original author's analysis, not everything
we think was present is detected, and *vice versa*, lots of unwanted
things are detected. These kinds of controls are perfect for discussing
what to use as a minimum abundance threshold - how many reads do we
need to declare a species as present in a sample?

Negative controls
-----------------

If you have called the provided ``setup.py`` to download the FASTQ files
and ``run.py`` to call THAPBI PICT, it would have used an optimistic
minimum abundance threshold of 10 copies of each unique sequence (the
default is a far more pesimitic 100).

This is not a good idea in general, but for your negative controls it
can be interesting to deliberately set no threshold, and accept any
sequence even if only supported by one read.

(Be sure to remove the intermediate FASTA files if you try this, as
otherwise THAPBI PICT would not replace the older higher threshold files).

If you do this, just how bad are the contamination levels? These little
tables were extracted manually from the sample level reports run with
``-a 1`` (accepting even sequences seen in only one read). The counts
are the total number of reads in each sample, while max is the highest
single sequence's abundance.

Amplicon library one, ITS1 using the BITS/B58S3 primer pair, samples
replicated in duplicate:

==================================== ============ ========== ===== ====
Description                          MiSeq-name   Accession  Count Max
------------------------------------ ------------ ---------- ----- ----
negative control from DNA extraction NegDNAA_S163 SRR5314317    91   48
negative control from DNA extraction NegDNAB_S175 SRR5314316   117   90
negative control from PCR step       NegPCRA_S187 SRR5314315  1041  994
negative control from PCR step       NegPCRB_S104 SRR5314314  3860 3574
==================================== ============ ========== ===== ====

Amplicon library two, ITS1 using the ITS1f/ITS2 primer pair:

============================= ============== ========== ===== ===
Description                   MiSeq-name     Accession  Count Max
----------------------------- -------------- ---------- ----- ---
negative control with GoTaq   NegCtlGoTq_S20 SRR5838526     2   1
negative control with Phusion NegCtlPhGn_S30 SRR5838523     8   4
negative control with reAmp   NegCtlPrmp_S10 SRR5838524     9   1
============================= ============== ========== ===== ===

Amplicon library two, ITS2 using the ITS3‐KYO2 and ITS4‐KYO3 primer pair:

============================= ============== ========== ===== ===
Description                   MiSeq-name     Accession  Count Max
----------------------------- -------------- ---------- ----- ---
negative control with GoTaq   NegCtlGoTq_S20 SRR5838526    14   2
negative control with Phusion NegCtlPhGn_S30 SRR5838523    17   4
negative control with PreAmp  NegCtlPrmp_S10 SRR5838524     5   1
============================= ============== ========== ===== ===

Looking at these numbers the levels in amplicon library two are commendably
low, at most four copies of any unique sequence - suggesting using a minimum
threshold of 10 here is quite reasonable.

Hereafter we will assume the minimum abundance threshold of 10 was used, and
you are encouraged to look at the sample or read level reports (e.g. in Excel)
while following along with this discussion.

However, the levels in amplicon library one are cause for concern.
Starting with the negative control from the DNA extraction (given a green
background in the Excel reports), we see both replicates had two unwanted
sequences. At the command line:

.. code:: console

    $ cat intermediate/AL1_BITS_B58S3/SRR5314317.fasta
    #left_primer:ACCTGCGGARGGATC
    #right_primer:GAGATCCRTTGYTRAAAGTT
    #raw_fastq:12564
    #flash:11641
    #cutadapt:112
    #abundance:98
    #threshold:10
    >d51507f661ebee38a85bec35b70b7ee1_64
    ATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCTTATGTTGCCTCGGCGGATCAGCCCGCGCCCCGTAAAAAG
    GGACGGCCCGCCGCAGGAACCCTAAACTCTGTTTTTAGTGGAACTTCTGAGTATAAAAAACAAATAAATCAA
    >daadc4126b5747c43511bd3be0ea2438_34
    ATTAGTGATTAATTATTGAGTGTAAAAACTCATAATCTTCTACAAACCACTGTTTTAAAATGTTTACTAGAATGTCCAAG
    CAGTTAAGCTGTTGAAATATTAAAAGTTTTATAA

Using a minimum of 10 has excluded lots of singletons etc here.


.. code:: console

    $ cat intermediate/AL1_BITS_B58S3/SRR5314316.fasta
    #left_primer:ACCTGCGGARGGATC
    #right_primer:GAGATCCRTTGYTRAAAGTT
    #raw_fastq:16297
    #flash:15829
    #cutadapt:132
    #abundance:111
    #threshold:10
    >d51507f661ebee38a85bec35b70b7ee1_101
    ATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCTTATGTTGCCTCGGCGGATCAGCCCGCGCCCCGTAAAAAG
    GGACGGCCCGCCGCAGGAACCCTAAACTCTGTTTTTAGTGGAACTTCTGAGTATAAAAAACAAATAAATCAA
    >e5b7a8b5dc0da33108cc8a881eb409f5_10
    ATTAGTGATTAATTATTGAGTGATAAAAAACTCATAATCTTCTACAAACCACTGTTTTAAAATGTTTACTAGAATGTCCA
    AGCAGTTAAGCTGTTGAAATATTAAAAGTTTTATAA

Using a minimum of 10 has again excluded lots of singletons.

As elsewhere in the documentation, the FASTA sequences have been line
wrapped at 80bp for display.

Both have ``d51507f661ebee38a85bec35b70b7ee1`` as their more common unwanted
sequence, a perfect match to *Fusarium graminearum* in the mock community.

The lower abundance sequence ``daadc4126b5747c43511bd3be0ea2438`` gives
perfect NCBI BLAST matches to several accessions of fungus *Wallemia muriae*),
likewise ``e5b7a8b5dc0da33108cc8a881eb409f5`` gives perfect NCBI BLAST matches
to *Wallemia muriae* and *Wallemia sebi*.

Moving on to the worst case, the negative control from the PCR reaction (given
a pale blue background in the Excel reports). Working at the terminal:

.. code:: console

    $ cat intermediate/AL1_BITS_B58S3/SRR5314315.fasta
    #left_primer:ACCTGCGGARGGATC
    #right_primer:GAGATCCRTTGYTRAAAGTT
    #raw_fastq:19406
    #flash:12140
    #cutadapt:1153
    #abundance:1085
    #threshold:10
    >d51507f661ebee38a85bec35b70b7ee1_1085
    ATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCTTATGTTGCCTCGGCGGATCAGCCCGCGCCCCGTAAAAAG
    GGACGGCCCGCCGCAGGAACCCTAAACTCTGTTTTTAGTGGAACTTCTGAGTATAAAAAACAAATAAATCAA

The minimum abundance excuded lots of singletons etc.

.. code:: console

    $ cat intermediate/AL1_BITS_B58S3/SRR5314314.fasta
    #left_primer:ACCTGCGGARGGATC
    #right_primer:GAGATCCRTTGYTRAAAGTT
    #raw_fastq:7285
    #flash:6128
    #cutadapt:4343
    #abundance:4017
    #threshold:10
    >d51507f661ebee38a85bec35b70b7ee1_3961
    ATTACCGAGTTTACAACTCCCAAACCCCTGTGAACATACCTTATGTTGCCTCGGCGGATCAGCCCGCGCCCCGTAAAAAG
    GGACGGCCCGCCGCAGGAACCCTAAACTCTGTTTTTAGTGGAACTTCTGAGTATAAAAAACAAATAAATCAA
    >716f6111ac2ee192c23282e07d23078a_25
    ATTCACAATGAGTGTTTTATGGCACTCTTTAAAAATCCATATCCACCTTGTGTGCAATGTTTGTTGGGAGAGAGGCAACT
    CTTTCCCTTCACAAATATCAACCTATATCTTTAACAACATTCGTCTGATAACATATTATGAATATACTTAATTCAAAATA
    T
    >5194a4ae3a27d987892a8fee7b1669b9_17
    ATTAGTGATTGCCTTTTATAAGGCTTTAAACTATATCCACTTACACCTGTGAACTGTTCTACTACTTGACGCAAGTCGAG
    TATTTTTACAAACAATGTGTAATGAACGTCGTTTTATTATAACAAAATAA
    >702929cef71042156acb3a28270d8831_14
    ATTACTGATTTGCTTAATTGCACCACATGTGTTTTTTATTGAACAAATTTCTTTGGTGGCGGGAGCAATCCTACCGCCAG
    AGGTTATAACTAAACCAAACTTTTTATTTACAGTCAAACTTGATTTATTATTACAATAGTCAA

The minimum abundance excluded lots of singletons etc. The vast majority of
these are slight variants of the dominant sequence (you can call ``thapbi_pict
edit-graph`` on this single sample to visualise this), and can thus be
explained as PCR noise.

Again, both samples have ``d51507f661ebee38a85bec35b70b7ee1`` as their main
(or only) unwanted sequence above the threshold, a perfect match to *Fusarium
graminearum* in the mock community.
Additionally ``716f6111ac2ee192c23282e07d23078a`` matched *Mortierella
verticillata* from the mock community.

Then ``5194a4ae3a27d987892a8fee7b1669b9`` gives perfect NCBI BLAST matches to
fungus *Trichosporon asahii* and ``702929cef71042156acb3a28270d8831`` to fungus
*Candida tropicalis*, which are unexpected contamination.

I concur with the author that the high levels of *Fusarium graminearum* are
most likely cross-contamination from the mock-community samples:

    Negative control samples in this sequencing run displayed some
    contamination by *F. graminearum*. This taxon was represented at slightly,
    but not dramatically, higher than expected relative abundances in the mock
    community samples; some of the increase over expected relative abundance
    may have been related to cross‐sample contamination.

Looking at the DNA extraction control alone, the THAPBI PICT default threshold
of 100 seems reasonable. However, if we set that aside the likely *Fusarium
graminearum* contamination, then the next worst contamination in any of these
four controls is at 32 copies, so you might argue 100 is a little harsh?

Certainly I think for amplicon library one, a threshold of 10 is too low, but
it could be defended for amplicon library two (where the controls had up to
four copies of an unwanted sequence).

Missing positive controls
-------------------------

We will look at the ratios later, but were all 19 species in the mock community
found? Perhaps the quickest way to answer this is to look at the classification
assessment output. At the command line, looking at the BLAST based classifier
as the most fuzzy of the three:

.. code:: console

    $ cut -f 1-5 summary/AL1_BITS_B58S3.assess.blast.tsv
    <SEE TABLE BELOW>

Or open this in Excel. You should find:

======================== === == === ==
#Species                 TP  FP FN  TN
======================== === == === ==
OVERALL                  345 5  168 71
Alternaria alternata     26  0  1   4
Aspergillus flavus       25  0  2   4
Candida apicola          27  0  0   4
Chytriomyces hyalinus    0   0  27  4
Claviceps purpurea       27  0  0   4
Fusarium graminearum     27  4  0   0
Fusarium oxysporum       27  0  0   4
Fusarium verticillioides 0   0  27  4
Mortierella verticillata 27  1  0   3
Naganishia albida        27  0  0   4
Neosartorya fischeri     24  0  3   4
Penicillium expansum     22  0  5   4
Rhizoctonia solani       19  0  8   4
Rhizomucor miehei        0   0  27  4
Rhizophagus irregularis  13  0  14  4
Saccharomyces cerevisiae 0   0  27  4
Saitoella complicata     27  0  0   4
Trichoderma reesei       27  0  0   4
Ustilago maydis          0   0  27  4
======================== === == === ==

Or, open this plain text tab separated Excel.

Five expected species were never found (zero true positives) at the 10 reads
abundance threshold: *Chytriomyces hyalinus*, *Fusarium verticillioides*,
*Rhizomucor miehei*, *Saccharomyces cerevisiae* and *Ustilago maydis*.

The author wrote:

    Two of the expected 19 phylotypes, *Fusarium verticillioides* and
    *Saccharomyces cerevisiae*, were not detected in any of the samples.
    A large number of reads, presumably including many *F. verticillioides*
    reads, were binned into a phylotype as unclassified *Fusarium*. The
    primers used in ITS1 amplification for this sequencing library match
    the rRNA gene sequence of *S. cerevisiae*. However, the expected ITS1
    amplicon length is 402 bases for this taxon, compared to a range of
    141‐330 bases across the remaining taxa in the mock community. Examining
    the data at earlier stages of processing revealed that *S. cerevisiae*
    was originally represented in the data set, but was completely removed
    during quality screening (Table S3).

    *Chytriomyes hyalinus*, *Rhizomucor miehei* and *Ustilago maydis* were
    detected at dramatically lower abundances than expected. Each of these
    taxa possesses sequence mismatches compared to the PCR primers that were
    used. The number of mismatches to the forward and reverse primers was as
    follows: for *C. hyalinus*, 2 and 1; for *R. miehei*, 0 and 2; and for
    *U. maydis*, 2 and 1. Thus, selection against these taxa may have been
    due to primer annealing efficiency.

That's pretty consistent (we've talked about *Fusarium verticillioides*
earlier), and suggests using a minimum abudance threshold of 10 in THAPBI
PICT is a little stricter that the author's pipeline.

Moving on to the second amplicon library, the larger ITS1 marker using the
ITS1f/ITS2 primer is more successful:

.. code:: console

    $ cut -f 1-5 summary/AL2_ITS1f_ITS2.assess.blast.tsv
    <SEE TABLE BELOW>

Or open this in Excel. You should find:

======================== === == === ==
#Species                 TP  FP FN  TN
======================== === == === ==
OVERALL                  398 0  115 57
Alternaria alternata     23  0  4   3
Aspergillus flavus       27  0  0   3
Candida apicola          12  0  15  3
Chytriomyces hyalinus    25  0  2   3
Claviceps purpurea       27  0  0   3
Fusarium graminearum     27  0  0   3
Fusarium oxysporum       27  0  0   3
Fusarium verticillioides 12  0  15  3
Mortierella verticillata 27  0  0   3
Naganishia albida        27  0  0   3
Neosartorya fischeri     23  0  4   3
Penicillium expansum     24  0  3   3
Rhizoctonia solani       24  0  3   3
Rhizomucor miehei        4   0  23  3
Rhizophagus irregularis  11  0  16  3
Saccharomyces cerevisiae 9   0  18  3
Saitoella complicata     27  0  0   3
Trichoderma reesei       25  0  2   3
Ustilago maydis          17  0  10  3
======================== === == === ==

Everything was found, although *Rhizomucor miehei* in particular found rarely,
followed by *Saccharomyces cerevisiae*. The original author wrote:

    The ITS1 data set yielded 18 of the expected 19 taxa (Tables S3, S5); as
    in the first library, no reads were classified as *F. verticillioides*,
    although many reads were placed in unclassified Fusarium. *Rhizomucor
    miehei* and *S. cerevisiae* were substantially underrepresented. Compared
    to primers ITS1f and ITS2, *R. miehei* had three mismatches in the forward
    and two mismatches in the reverse. *Saccharomyces cerevisiae* had one
    mismatch in the forward primer and again likely suffered negative bias
    associated with amplicon length (Table 3) and low sequence quality
    (Table S3).

Again, broad agreement here, with the problem of *Fusarium verticillioides*
discussed earlier.

And finally, amplicon library two for ITS2 using the ITS3-KYO2 and ITS4-KYO3
primers:

.. code:: console

    $ cut -f 1-5 summary/AL2_ITS3-KYO2_ITS4-KYO3.assess.blast.tsv
    <SEE TABLE BELOW>

Or open this in Excel. You should find:

======================== === == === ==
#Species                 TP  FP FN  TN
======================== === == === ==
OVERALL                  313 0  200 57
Alternaria alternata     16  0  11  3
Aspergillus flavus       24  0  3   3
Candida apicola          0   0  27  3
Chytriomyces hyalinus    0   0  27  3
Claviceps purpurea       23  0  4   3
Fusarium graminearum     27  0  0   3
Fusarium oxysporum       27  0  0   3
Fusarium verticillioides 27  0  0   3
Mortierella verticillata 12  0  15  3
Naganishia albida        27  0  0   3
Neosartorya fischeri     16  0  11  3
Penicillium expansum     23  0  4   3
Rhizoctonia solani       11  0  16  3
Rhizomucor miehei        0   0  27  3
Rhizophagus irregularis  5   0  22  3
Saccharomyces cerevisiae 27  0  0   3
Saitoella complicata     26  0  1   3
Trichoderma reesei       22  0  5   3
Ustilago maydis          0   0  27  3
======================== === == === ==

This time we're missing *Candida apicola*, *Chytriomyces hyalinus*,
*Rhizomucor miehei* and *Ustilago maydis*.

This too is in board agreement with the original author, although
*Candida apicola* must have just dipped below our abundance threshold.

    Different amplification biases were evident between the ITS1 and ITS2
    loci. In the ITS2 data set, only 16 of the 19 taxa that were present
    could be detected; *C. hyalinus*, *R. miehei* and *U. maydis* were not
    observed (Tables S3, S6). ...
    *Rhizomucor miehei* has one mismatch to the forward primer and three
    mismatches to the reverse primer. While neither *C. hyalinus* nor
    *U. maydis* have sequence mismatches compared to the primers, these two
    taxa have longer ITS2 amplicons than any others in the mock community
    (Table 3). These two taxa were originally represented with a small number
    of reads in the raw data, but were completely removed during quality
    screening (Table S3). *Candida apicola*, which possesses two mismatches
    to the reverse primer for this amplicon, was detected at substantially
    lower than expected frequencies (Figure 7; Figures S5, S6).

So, using THAPBI PICT on these amplicon datasets with a minimum abundance
threshold of 10 gives broad agreement with the original analysis.
