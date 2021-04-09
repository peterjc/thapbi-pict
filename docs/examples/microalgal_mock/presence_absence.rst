Presence and absence
====================

This example includes mock communities which are a controlled setup where we
know what the classifier ought ideally to report for every sample - and all
their expected marker sequences are in the classification database.

Of course, just as in the original author's analysis, not everything we think
was present is detected. And *vice versa*, we see some things which are not
classified.

The experiment had a negative PCR control, but it was not sequenced. The two
different amplicons were sequenced separately, but since the FASTQ files are
provided pre-trimmed, we can't try using the other primer pair on each sample
as a negative control.

Nevertheless, this kind of data is important for discussing what to use as a
minimum abundance threshold - how many reads do we need to declare a species
as present in a sample?

18S rRNA V4 marker
------------------

If you have called the provided ``setup.py`` to download the FASTQ files and
``run.py`` to call THAPBI PICT, it would have used the default minimum
abundance threshold of 100 copies of each unique sequence.

Given we know the expected species makeup, the classification assessment output
is useful - this shows the expected 12 species:

.. code:: console

    $ cut -f 1-5 summary/V4.assess.onebp.tsv
    <SEE TABLE BELOW>

Or, open the file ``V4.assess.onebp.tsv`` in Excel. You should see:

================================================ === == == ==
#Species                                         TP  FP FN TN
================================================ === == == ==
OVERALL                                          134 0  82 36
Chlorella vulgaris strain CCAP 211/11B           15  0  3  3
Cryptomonas pyrenoidifera strain CCAP 979/61     15  0  3  3
Heterocapsa niei strain UTEX LB 2722             14  0  4  3
Isochrysis galbana strain CCAP 927/1             0   0  18 3
Nannochloropsis oculata strain CCMP225           7   0  11 3
Ochromonas sp. UTEX LB 2575                      12  0  6  3
Prymnesium parvum strain CCAP 946/6              0   0  18 3
Rhodomonas sp. CCAP 995/5                        15  0  3  3
Symbiodinium microadriaticum strain UTEX LB 2281 14  0  4  3
Tetradesmus obliquus strain CCAP 276/3A          15  0  3  3
Thalassiosira pseudonana strain CCAP 1085/12     12  0  6  3
Trebouxia sp. CCAP 213/3                         15  0  3  3
================================================ === == == ==

Notice 2 of the 12 species were not detected, *Isochrysis galbana* and
*Prymnesium parvum* exactly as per the authors' analysis, attributed to
mismatches in the reverse primer. I also note *Nannochloropsis oculata*
detection is comparatively low.

There are no false positives, but by the construction of the database
the only way that could happen would be a marine species in the pure
freshwater community (M1), or *vice versa* seeing a freshwater species
in the pure marine community (M7).

Please note this report does not include any unknown classifications,
those do appear in the main read and sample based reports.

We have here seven mock communities covering staggered ratios from six
freshwater species only though to six marine species only. Happily this
is reflected in the sample level reports. Open ``V4.samples.onebp.xlsx``
in a spreadsheet, or if you prefer the command line:

.. code:: console

    $ grep  -E "^(#|Mock)" summary/V4.samples.onebp.tsv | cut -f 2,5,12-23
    ...

With light editing for display you should get:

=== ====== ===== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====
ID  Ratio  Reads Chl. Cry. Het. Nan. Och. Rho. Sym. Tet. Tha. Tre. Unk.
--- ------ ----- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
MC1 100:0  17383 1856 1659    0    0 1301 5507    0 2629    0 4431    0
MC1 100:0  19382 2044 1807    0    0 1399 6302    0 2834    0 4996    0
MC1 100:0  18538 2448 1616    0    0 1269 5230    0 3348    0 4627    0
MC2 100:1  16069 2117 1347    0    0 1069 4341    0 3188    0 4007    0
MC2 100:1  11564 1511 1032    0    0  847 3381    0 2182    0 2611    0
MC2 100:1  11822 1450  986    0    0  815 3530    0 2270    0 2771    0
MC3 100:10  7822 1089  809  113    0  528 1973  111 1625    0 1574    0
MC3 100:10  9059 1086  771  180    0  683 3017  130 1397    0 1795    0
MC3 100:10 10591 1341 1037    0    0  865 3068    0 2051    0 2229    0
MC4 1:1    10383  785  711 1856    0  441 2271 1091 1070  498 1660    0
MC4 1:1    11107  711  592 1893    0  530 2023 1651 1348 1084 1275    0
MC4 1:1    10185  928  750 1036    0  685 2333  789 1416  612 1636    0
MC5 10:100 12187  147  117 4826  223    0  166 3200  175 3163  170    0
MC5 10:100 13477  195  134 5077  152    0  339 3598  225 3464  293    0
MC5 10:100 13247  131  159 5126  184    0  210 3710  211 3274  242    0
MC6 1:100  12482    0    0 5271  400    0    0 3690    0 3121    0    0
MC6 1:100  13254    0    0 5281  216    0    0 3961    0 3796    0    0
MC6 1:100  14128    0    0 6073  147    0    0 4335    0 3573    0    0
MC7 0:100   6415    0    0 2671    0    0    0 1943    0 1801    0    0
MC7 0:100   8848    0    0 3982    0    0    0 2542    0 2324    0    0
MC7 0:100  11205    0    0 4686  133    0    0 3323    0 3063    0    0
=== ====== ===== ==== ==== ==== ==== ==== ==== ==== ==== ==== ==== ====

In the Excel file, focus on the rows from the mock community only (yellow
background in the Excel file).

It should be clear that some species abundance increases with the ratio, while
others decrease. If we sort this into the six freshwater and four marine
species (the missing two species were both marine), the pattern is even clearer
as in this Excel screenshot (freshwater in pale blue, marine in dark blue):

.. image:: https://gist.githubusercontent.com/peterjc/3baeb3a648b8cdd7cfa970920eaf7f38/raw/4061bc7aafe7360e3c1b90ab82e92cffaaad02ff/V4.samples.onebp.svg?sanitize=true
   :alt: Excel screenshot from a cropped and sorted THAPBI PICT sample report

Notice there are no unknowns - this example seems to behave perfectly, and at
this threshold minority samples present at the 1:100 ratio are not found.

18S rRNA V8V9 marker
--------------------

Proceeding along the same basis:

.. code:: console

    $ cut -f 1-5 summary/V8V9.assess.onebp.tsv
    <SEE TABLE BELOW>

Or, open the file ``V8V9.assess.onebp.tsv`` in Excel. You should see:

================================================ === == == ==
#Species                                         TP  FP FN TN
================================================ === == == ==
OVERALL                                          131 0  85 36
Chlorella vulgaris strain CCAP 211/11B           12  0  6  3
Cryptomonas pyrenoidifera strain CCAP 979/61     14  0  4  3
Heterocapsa niei strain UTEX LB 2722             12  0  6  3
Isochrysis galbana strain CCAP 927/1             12  0  6  3
Nannochloropsis oculata strain CCMP225           5   0  13 3
Ochromonas sp. UTEX LB 2575                      12  0  6  3
Prymnesium parvum strain CCAP 946/6              12  0  6  3
Rhodomonas sp. CCAP 995/5                        14  0  4  3
Symbiodinium microadriaticum strain UTEX LB 2281 12  0  6  3
Tetradesmus obliquus strain CCAP 276/3A          12  0  6  3
Thalassiosira pseudonana strain CCAP 1085/12     0   0  18 3
Trebouxia sp. CCAP 213/3                         14  0  4  3
================================================ === == == ==

Very similar, but time *Thalassiosira pseudonana* is missing with the default
strict ``onebp`` classifier - and the authors did not report its absence. The
reason for this is discussed shortly.

Open ``V8V9.samples.onebp.xlsx`` and focus on the mock community rows (yellow
background). Again, in the following screen shot we have sorted the columns
into freshwater (six in pale blue) and marine (five in dark blue):

.. image:: https://gist.githubusercontent.com/peterjc/3baeb3a648b8cdd7cfa970920eaf7f38/raw/4061bc7aafe7360e3c1b90ab82e92cffaaad02ff/V8V9.samples.onebp.svg?sanitize=true
   :alt: Excel screenshot from a cropped and sorted THAPBI PICT sample report

Again, much the same picture *except* all the communities (although not all
the replicates at the freshwater end) report unknown sequences.

If you open ``V8V9.reads.onebp.xlsx`` you can see all the sequences not
assigned a species. If you delete the other samples and re-sort, you should
find this is the most abundance unassigned sequence (column 1 is the MD5
checksum, column 2 is blank for no classification, column 3 is the sequence)::

    >e586634fe2e532c7fd59b071493579c3
    TAGATGTCCTGGGCCGCACGCGCGCTACACTGATGCACTCAACGAGCATATAACCTTGGCCGAGAGGCCTGGGTAATCTT
    GTTAACATGCATCGTGATAGGGATAGATTATTGCAATTATTAATCTTGAACGAGGAATTCCTAGTAATCGCAGATCATCA
    ATCTGCAATGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCACCTACCGATTGGATGGTCCGGTGAGGAGTCGAG
    ATTGTGGCCTGGTTCCTTTATTGGGATTTGGCTACGAGAACTTCTCCAAACCTTATCATCTAGAGGAAGGTGAAGTCGTA
    ACAAGGTTTCC

Running an NCBI BLAST search online gives perfect full length hits to multiple
*Thalassiosira pseudonana* accessions (and *Thalassiosira oceanica* too). It
also gives a perfect but partial match to KU900218.1 which is the accession
from the authors for this mock community member, and is in our database:

.. code:: console

    $ thapbi_pict dump -d V8V9.sqlite -g Thalassiosira -f fasta
    >KU900218.1 Thalassiosira pseudonana strain CCAP 1085/12
    TAGATGTCCTGGGCCGCACGCGCGCTACACTGATGCACTCAACGAGCATATAACCTTGGCCGAGAGGCCTGGGTAATCTT
    GTTAACATGCATCGTGATAGGGATAGATTATTGCAATTATTAATCTTGAACGAGGAATTCCTAGTAATCGCAGATCATCA
    ATCTGCAATGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCACCTACCGATTGGATGGTCCGGTGAGGAGTCGAG
    ATTGTGGCCTGGTTCCTTTATTGGGATTTGGCTACGAGAACTTCTCCAAACCTTATCATCTAGAGGAAG
    Wrote 1 fasta format entries

Unfortunately it seems to be incomplete, missing the last 22bp,
``GTGAAGTCGTAACAAGGTTTCC`` - as can be seen on the full length record (e.g.
`KU900218.1 on NCBI <https://www.ncbi.nlm.nih.gov/nucleotide/KU900218.1>`_).

The second and third most common unexpected sequence are both marine
associated::

    >64f90363dd2c1f85645af55a92d4c376
    TAGATGTTCTGGGCTGCACGCGCGCTACACTGATGCGCTCAACGAGTTTATGACCTTGCCCGGAAGGGTTGGGTAATCTT
    CTTAAAACGCATCGTGATGGGGATAGATTATTGCAATTATTAATCTTCAACGAGGAATTCCTAGTAAGCGCGAGTCATCA
    GCTCGTGCTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTCCTACCGATTGAGTGATCCGGTGAATAATTCGG
    ACTGACGCAGTGCTCAGCTTCTGGACGTTGCGTTGGAAAGCTTCATGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAA
    CAAGGTTTCC
    >1dac8fc1b9b2736a190333d1b5a25056
    TAGATGTCCTGGGCTGCACGCGCGCTACACTGATGCGCTCAACGAGTTTTTGATCTTGCCTGAAATGGCTGGGTAATCTT
    TTTAAAATGCATCGTGATGGGGATAGATCATTGCAATTATTGATCTTCAACGAGGAATTCCTAGTAAGCGCGAGTCATCA
    GCTCGTGCTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTCCTACCGATTGAGTGATCCGGTGAATAATTCGG
    ACTGCAGCAGTGTTCGGTCACGAACGTTGCAGCGGAAAGTTTAGTGAACCTTATCACTTAGAGGAAGGAGAAGTCGTAAC
    AAGGTTTCC

Running an NCBI BLAST search online gives KU900226.1 *Symbiodinium
microadriaticum* and KU900227.1 *Heterocapsa niei* respectively amongst their
top hits - both accessions from the mock community - but at only 97% identify.
These could be a secondary variant copies in those genomes?

The fourth and fifth most common unexpected sequence are both freshwater
associated::

    >935f0cd55155d64af97bef8245b36f4d
    TAGATGTTCTGGGCCGCACGCGCGCTACACTGATGGATGCAACGAGCTCTACCCTTGACCGAAAGGCCCGGGTAAACTTG
    TCAAAATCCATCGTGATGGGGATAGATTATTGCAATTCTTGATCTTCAACGAGGAATTCCTAGTAAGCGCGAGTCATCAG
    CTCGCGTTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTCCTACCGATTGAATGGTCCGGTGAAATCTTCGGA
    TTGCTGACTTTGGCATTTATTTGTCTTAGTTGCGAGAAGTTGATTGAACCTTATCATTTAGAGGAAGGAGAAGTCGTAAC
    AAGGTTTCC
    >065d987074eac615c1a2292b2a6e2680
    TAGATGTTCTGGGCCGCACGCGCGCTACACTGATGAATGCAACGAGCTCCTCCCTTATTCGAAAGAATCGGGTAAACTTG
    TGAAAATTCATCGTGATGGGGATAGATTATTGCAATTATTAATCTTCAACGAGGAATTCCTAGTAAGCGCGAGTCATCAG
    CTCGCGTTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTCCTACCGATTGAATGGTCCGGCGAAATCTTCAGA
    TTGCTGGCGAGTTCTTCACGGTTCTCGCTGTGAGAAGTTGATTAAACCTTATCATTTAGAGGAAGGAGAAGTCGTAACAA
    GGTTTCC

Running an NCBI BLAST search online gives KU900223.1 *Rhodomonas* sp. CCAP
995/5, and KU900222.1 *Cryptomonas pyrenoidifera* respectively amongst their
top hits - again both accessions from the mock community - but at only 96%
identify. Likewise the sixth most common unexpected sequence has a 97% identity
match to KU900220.1 *Chlorella vulgaris*::

    >1610815dfa2c7d4b602e84114d15cd51
    TAGATGTTCTGGGCCGCACGCGCGCTACACTGATGCATTCAACAAGCCTATCCTTGACCGAAGGGTCTGGGTAATCTTTG
    AAACTGCATCGTGATGGGGATAGATTATTGCAATTATTAGTCTTCAACGAGGAATGCCTAGTAAGCGCAAGTCATCAGCT
    TGCGTTGATTACGTCCCTGCCCTTTGTACACACCGCCCGTCGCTCCTACCGATTGGGTGTGCTGGTGAAGTGTTCGGATT
    GGCGACCTGGGGCGGTCTCCGCTCTCGGCCGCCGAGAAGTTCATTAAACCCTCCCACCTAGAGGAAGGAGAAGTCGTAAC
    AAGGTTTCC

The remaining unexpected samples appeared in only one of the sequenced sample
replicates, but again could be BLAST matched to the mock community.

Conclusion
----------

Based on this initial examination, and looking at the edit-graph structure,
both markers should work with our default ``onebp`` classifier (looking for a
perfect match or at most 1bp away). For the V8V9 marker, it appears the
database would benefit from including secondary sequences from the reference
strains too.

For either marker, applying THAPBI PICT to any environmental samples will need
the database extended. For now, looking at ``V4.samples.onebp.xlsx``, the only
species assigned to the environmental samples was ``Tetradesmus obliquus`` in
the freshwater marsh (samples 4F and 5F), and wastewater from Urbana IL WWTP
(samples 8W, 9W and 10W). Likewise in ``V8V9.samples.onebp.xlsx``, but only in
samples 4F and 10W.

We refer you to the original paper for a much more detailed discussion of the
relative merits of these two primer sets for microalgae.
