Minimum Abundance Threshold
===========================

THAPBI PICT has a default minimum abundance threshold of 100 reads per sample
before accepting any unique sequence. Background contamination and PCR noise
levels will vary, so having multiple :ref:`negative_controls` will help set
this objectively.

In this dataset there is a single negative control for the MOL16S marker,
library ``BIM8M`` aka ``SRR5534986``. However, we can also treat all the
SPH16S libraries as negative controls for the MOL16S marker, and vice versa.
You could do this automatically within THAPBI PICT via the ``-n`` or
``--negctrls`` command line option, but as we shall see in this example it
will discard most of the data.

In order to examine an appropriate minimum abundance threshold, the ``run.sh``
script provided uses ``-a 10`` to accept any unique sequence seen in sample at
least ten times. This *does* allow unwanted noise though to the reports.

SPH16S
------

This was the more specific primer pair, expected to only amplify sphaeriid
mussel species, so in general we expect less unique sequences than with the
more general MOL16S primers.

Looking at some key columns in the sample report,

.. code:: console

    $ cut -f 1,2,4,6,9-10 SPH16S.samples.onebp.tsv
    ...

Or, open ``SPH16S.samples.onebp.xlsx`` in Excel. Focusing on the the left hand
columns, you should see:

======== =============== ============= ========== ========= =========
Marker   Group           library_name  Raw FASTQ  Cutadapt  Seq-count
-------- --------------- ------------- ---------- --------- ---------
MOL16S   Aquarium        BIR2M         306311     8         0
MOL16S   Aquarium        BIR6M         291954     24        0
MOL16S   Control         BIM8M         2433       0         0
MOL16S   Mock Community  SC3PRO1       689712     72        0
MOL16S   Mock Community  SC3PRO2       405048     70        0
MOL16S   Mock Community  SC3PRO3       402219     32        0
MOL16S   Mock Community  NFSC3PRO3     349590     41        10
MOL16S   Mock Community  SC3PRO4       671241     27        0
MOL16S   Mock Community  NFSC3PRO4     420015     18        0
MOL16S   Mock Community  SC3PRO5       480606     55        0
MOL16S   River           BIM6M         821849     0         0
MOL16S   River           BIM2M         1119271    40        0
MOL16S   River           BIM4M         709472     46        19
SPH16S   Aquarium        BIR2S         498926     251724    209358
SPH16S   Aquarium        BIR6S         240360     226083    191393
SPH16S   Mock Community  SPSC3PRO1     425271     318149    224510
SPH16S   Mock Community  SPSC3PRO2     341476     282623    204137
SPH16S   Mock Community  SPSC3PRO4     410780     304178    197340
======== =============== ============= ========== ========= =========

Things to note:

* In the "Raw FASTQ" column, the control has far fewer raw reads (good).
* The "Cutadapt" column shows reads after SPH16S primer trimming. There are
  hundreds of thousands for the final five samples amplified with these
  primers (good). The first 13 samples were amplified with the MOL16S primers,
  but still have low levels of sequences matching the SPH16S primers (bad).
* The "Seq-count" column is after applying the minimum abundance threshold
  (here 10). Two negative controls still have reads, lifting the threshold
  to 20 or more would fix this. These are *Sphaerium simile* in mock community
  ``NFSC3PRO3``, and an unknown in river sample ``BIM4M``.

So, using the MOL16S samples as negative controls suggests that for the SPH16S
the default minimum abundance threshold is perhaps overly harsh - but using
*at least* 20 would be wise.

MOL16S
------

We'll initially looking at the same key columns in the sample report,

.. code:: console

    $ cut -f 1,2,4,6,9-10 MOL16S.samples.onebp.tsv
    ...

Or, open ``MOL16S.samples.onebp.xlsx`` in Excel. Focusing on the the left hand
columns, you should see:

======== =============== ============= ========== ========= =========
Marker   Group           library_name  Raw FASTQ  Cutadapt  Seq-count
-------- --------------- ------------- ---------- --------- ---------
MOL16S   Aquarium        BIR2M         306311     279297    240659
MOL16S   Aquarium        BIR6M         291954     265476    238141
MOL16S   Control         BIM8M         2433       928       513
MOL16S   Mock Community  SC3PRO1       689712     247703    210045
MOL16S   Mock Community  SC3PRO2       405048     102735    85080
MOL16S   Mock Community  SC3PRO3       402219     116643    94799
MOL16S   Mock Community  NFSC3PRO3     349590     105425    85476
MOL16S   Mock Community  SC3PRO4       671241     168201    130340
MOL16S   Mock Community  NFSC3PRO4     420015     116666    75975
MOL16S   Mock Community  SC3PRO5       480606     152834    129045
MOL16S   River           BIM6M         821849     764830    673378
MOL16S   River           BIM2M         1119271    886766    767797
MOL16S   River           BIM4M         709472     342913    296602
SPH16S   Aquarium        BIR2S         498926     30        0
SPH16S   Aquarium        BIR6S         240360     30        0
SPH16S   Mock Community  SPSC3PRO1     425271     40        14
SPH16S   Mock Community  SPSC3PRO2     341476     2184      805
SPH16S   Mock Community  SPSC3PRO4     410780     102       16
======== =============== ============= ========== ========= =========

Looking at the same points, we see two problems:

* The control sample BIM8M (SRR5534986) had almost a thousand unwanted MOL16S
  matches, reduced to 513 with a minimum abundance threshold of 10.

* All the SPH16S mock community samples have unwanted MOS16S matches, the
  worst case being SPSC3PRO2 (SRR5534981) with over two thousand reduced to
  805 with the minimum abundance threshold of 10.

To see exactly what is in these two problematic samples, we can turn to the
read report - or look directly at the intermediate FASTA files:

.. code:: console

    $ head -n 14 intermediate/MOL16S/SRR5534986.fasta
    #left_primer:RRWRGACRAGAAGACCCT
    #right_primer:ARTCCAACATCGAGGT
    #raw_fastq:2433
    #trimmomatic:2306
    #flash:1837
    #cutadapt:928
    #abundance:513
    #threshold:10
    >20c0669e4c6f8436c9d42736df727c83_433
    ATCGAACTTAAATTATTTGTTTAAATTTTTAAATAGAAAAAGTTTAGTTGGGGAAACTTAAAGTAAAAGGTAACGCTTTA
    TTTTTTTGTCAGGAGCCTGTAGTATGGAAAAATGAAAAAGTTACCGTAGGGATAACAGCGCTTTCTTCTCTGAGAGGACT
    AATTAAAGAGTTGGTTGCG
    >a36d3f7291c173c4243f22c2afbd111e_47
    ATCGAACTTAAATTATTTGTTTAAATTTTTAAATAGAAAAAGTTTAGTTGGGGAAACTTAAAGTAAAAGGTAACGCTTTA
    TTTTTTTGTCAGGAGCCTGTAGTATGGAAAAATGAAAAAGTTACCGTAGGGATAACAGCGCTTTCTTCTCTGAGAGGATT
    AATTAAAGAGTTGGTTGCG
    >e1d838b4f39bffe88d8c0e79b52700f1_13
    ATCGAACTTAAATTATTTGTTTAAATTTTTAAATAGAAAAAGTTTAGTTGGGGAAACTTAAAGTAAAAGGTAACGCTTTA
    TTTTTTGTCAGGAGCCTGTAGTATGGAAAAATGAAAAAGTTACCGTAGGGATAACAGCGCTTTCTTCTCTGAGAGGACTA
    ATTAAAGAGTTGGTTGCG

The unwanted sequences in the control sample are dominated by a single
sequence (and variants of it; shown line wrapped at 80 characters), which was
matched to *Sphaerium simile*.

This is consistent with the original author's analysis - although our pipeline
has produced higher read counts:

    Finally, our water blank sample had 71 reads, eight of those being
    singletons with the remaining belonging to *Sphaerium striatinum*
    (Table 9), likely due to amplicon contamination in the lab.

What about the other problematic sample? Again, you can find this in the Excel
read report, or at the command line:

.. code:: console

    $ head -n 14 intermediate/MOL16S/SRR5534981.fasta
    #left_primer:RRWRGACRAGAAGACCCT
    #right_primer:ARTCCAACATCGAGGT
    #raw_fastq:341476
    #trimmomatic:341289
    #flash:314812
    #cutadapt:2184
    #abundance:805
    #threshold:10
    >abb4c9d82203b201ff91fc87b7c4e337_623
    ATCGAACTTGAATTGTGTGTTTTAGTTTTGGAATACAGAAAGTTTAGTTGGGGAAACTTAAAGTTAAGAAAAACGCTTTT
    TTGTTATAAAATGATCCTGTATTATAGAAAAATGAAAAAGTTACCGTAGGGATAACAGCGCTTTCTTCTCTGAGAGGACT
    AATCAAAGAGTTGGTTGCGACCTCGATGTTCGTACATCTAGT
    >dfb3668f028fad9ea3df1408f56c90b2_70
    ATCGAACTTGAATTGTGTGTTTTAGTTTTGGAATACAGAAAGTTTAGTTGGGGAAACTTAAAGTTAAGAAAAACGCTTTT
    TTGTTATAAAATGATCCTGTATTATAGAAAAATGAAAAAGTTACCGTAGGGATAACAGCGCTTTCTTCTCTGAGAGGACT
    AATCAAAGAGTTGGTTGCGACCTCGATGTTCGTATATCTAGT
    >4ffe4f9f031bea2734d75e8b6e55a5d5_29
    ATCGAACTTGAATTGTGTGTTTTAGTTTTGGAATACAGAAAGTTTAGTTGGGGAAACTTAAAGTTAAGAAAAACGCTTTT
    TTGTTATAAAATGATCCTGTATTATAGAAAAATGAAAAAGTTACCGTAGGGATAACAGCGCTTTCTTCTCTGAGAGGACT
    AATCAAAGAGTTGGTTGCGACATCGATGTTCGTACATCTAGT

The unwanted mock community sample is again dominated by a single sequence,
which was not matched in the database constructed for this example. NCBI BLAST
identifies it as *Pisidium compressum*, giving a perfect match if we discard
the final 12bp. This is one of the control species in the mock community, but
recall the amplified regions of the MOL16S and SPH16S primers overlap...

This sequence appears be part of a longer unwanted product of the SPH16S_F
primer (``TAGGGGAAGGTATGAATGGTTTG`` - should be present here) and MOL16S_R
primer (``ARTCCAACATCGAGGT`` - not not be present here), which can be trimmed
to look like either a SPH16S product *or* a MOL16S product.

.. code:: console

    $ head intermediate/large/SRR5534981.fasta
    #left_primer:TAGGGGAAGGTATGAATGGTTTG
    #right_primer:ARTCCAACATCGAGGT
    #raw_fastq:341476
    #trimmomatic:341289
    #flash:314812
    #cutadapt:2237
    #abundance:584
    #threshold:10
    >c40a4b99f05302d2fecdbc3b5f619c54_462
    ACGTGGGAAAAGCTGTCTCTTTTATATAGAAAGAAGTTTATTTTTGAGTGAAAAAGCTTAAATATTTGTAAAAGACGAGA
    AGACCCTATCGAACTTGAATTGTGTGTTTTAGTTTTGGAATACAGAAAGTTTAGTTGGGGAAACTTAAAGTTAAGAAAAA
    CGCTTTTTTGTTATAAAATGATCCTGTATTATAGAAAAATGAAAAAGTTACCGTAGGGATAACAGCGCTTTCTTCTCTGA
    GAGGACTAATCAAAGAGTTGGTTGCGACCTCGATGTTCGTACATCTAGT

This longer sequence (shown here with line wrapping at 80 characters) again
matches *Pisidium compressum* (ignoring the last 12 bases).

Running THAPBI PICT with this primer pair (as done in the ``run.sh`` script)
reveals that the only other sample with this kind of primer mixing is
SRR5534978 aka SPSC3PRO1, with an unwanted long sequence seen 10 times.

.. code:: console

    $ cat intermediate/large/SRR5534978.fasta
    #left_primer:TAGGGGAAGGTATGAATGGTTTG
    #right_primer:ARTCCAACATCGAGGT
    #raw_fastq:425271
    #trimmomatic:425236
    #flash:395494
    #cutadapt:149
    #abundance:10
    #threshold:10
    >f520da824d259a518c08d2f4ec46eaf3_10
    ACGTGGAAAAAACTGTCTCTTTTGTATAAAAAGAAGTTTATTTTTAAGTGAAAAAGCTTGAATGTTTATAAAAGACGAGA
    AGACCCTATCGAACTTAAATTATTTGTTTAAATTTTTAAATAGAAAAAGTTTAGTTGGGGAAACTTAAAGTAAAAGGTAA
    CGCTTTATTTTTTTGTCAGGAGCCTGTAGTATGGAAAAATGAAAAAGTTACCGTAGGGATAACAGCGCTTTCTTCTCTGA
    GAGGACTAATTAAAGAGTTGGTTGCG

Note this is shown with the sequence line wrapped at 80 characters.

Minimum threshold
-----------------

Clearly using a minimum abundance threshold of 10 is too low, and it should be
increased to at least 20 based on this. However, we have the two exceptional
sequences present at over 500 copies. Setting the minimum that high seems
excessive - but perhaps the THAPBI PICT default of 100 is more reasonable?
