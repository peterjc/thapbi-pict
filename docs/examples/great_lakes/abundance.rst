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

    $ cut -f 1,2,4,7,9-10 summary/SPH16S.samples.onebp.tsv
    <SEE TABLE BELOW>

Or, open ``SPH16S.samples.onebp.xlsx`` in Excel. Focusing on the the left hand
columns, you should see:

======= ============== ============ ========= ======== ==========
#Marker Group          Library-name Raw FASTQ Cutadapt Read count
======= ============== ============ ========= ======== ==========
MOL16S  Aquarium       BIR2M        306311    22       0
MOL16S  Aquarium       BIR6M        291954    29       0
MOL16S  Control        BIM8M        2433      0        0
MOL16S  Mock Community SC3PRO1      689712    2693     0
MOL16S  Mock Community SC3PRO2      405048    2063     0
MOL16S  Mock Community SC3PRO3      402219    1007     0
MOL16S  Mock Community NFSC3PRO3    349590    853      10
MOL16S  Mock Community SC3PRO4      671241    974      0
MOL16S  Mock Community NFSC3PRO4    420015    533      0
MOL16S  Mock Community SC3PRO5      480606    1396     0
MOL16S  River          BIM6M        821849    2        0
MOL16S  River          BIM2M        1119271   447      0
MOL16S  River          BIM4M        709472    50       19
SPH16S  Aquarium       BIR2S        498926    251731   209435
SPH16S  Aquarium       BIR6S        240360    226084   191469
SPH16S  Mock Community SPSC3PRO1    425271    318150   224718
SPH16S  Mock Community SPSC3PRO2    341476    282637   204289
SPH16S  Mock Community SPSC3PRO4    410780    304194   197522
======= ============== ============ ========= ======== ==========

Things to note:

* In the "Raw FASTQ" column, the control has far fewer raw reads (good).
* The "Cutadapt" column shows reads after SPH16S primer trimming. There are
  hundreds of thousands for the final five samples amplified with these
  primers (good). The first 13 samples were amplified with the MOL16S primers,
  but still have low levels of sequences matching the SPH16S primers (bad).
* The "Read count" column is after applying the minimum abundance threshold
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

    $ cut -f 1,2,4,7,9-10 summary/MOL16S.samples.onebp.tsv
    <SEE TABLE BELOW>

Or, open ``MOL16S.samples.onebp.xlsx`` in Excel. Focusing on the the left hand
columns, you should see:

======= ============== ============ ========= ======== ==========
#Marker Group          Library-name Raw FASTQ Cutadapt Read count
======= ============== ============ ========= ======== ==========
MOL16S  Aquarium       BIR2M        306311    297738   256466
MOL16S  Aquarium       BIR6M        291954    286488   256527
MOL16S  Control        BIM8M        2433      1010     551
MOL16S  Mock Community SC3PRO1      689712    656795   550336
MOL16S  Mock Community SC3PRO2      405048    377068   297912
MOL16S  Mock Community SC3PRO3      402219    380395   304641
MOL16S  Mock Community NFSC3PRO3    349590    328997   262983
MOL16S  Mock Community SC3PRO4      671241    628747   494294
MOL16S  Mock Community NFSC3PRO4    420015    364291   262739
MOL16S  Mock Community SC3PRO5      480606    459043   383892
MOL16S  River          BIM6M        821849    799497   703741
MOL16S  River          BIM2M        1119271   954741   823977
MOL16S  River          BIM4M        709472    367498   317391
SPH16S  Aquarium       BIR2S        498926    33       0
SPH16S  Aquarium       BIR6S        240360    35       0
SPH16S  Mock Community SPSC3PRO1    425271    64       14
SPH16S  Mock Community SPSC3PRO2    341476    2322     833
SPH16S  Mock Community SPSC3PRO4    410780    432      108
======= ============== ============ ========= ======== ==========

Looking at the same points, we see two problems:

* The control sample BIM8M (SRR5534986) had almost a thousand unwanted MOL16S
  matches, reduced to 551 with a minimum abundance threshold of 10.

* All the SPH16S mock community samples have unwanted MOS16S matches, the
  worst case being SPSC3PRO2 (SRR5534981) with over two thousand reduced to
  833 with the minimum abundance threshold of 10.

To see exactly what is in these two problematic samples, we can turn to the
read report - or look directly at the intermediate FASTA files:

.. code:: console

    $ head -n 13 intermediate/MOL16S/SRR5534986.fasta
    #left_primer:RRWRGACRAGAAGACCCT
    #right_primer:ARTCCAACATCGAGGT
    #raw_fastq:2433
    #flash:1963
    #cutadapt:1010
    #abundance:551
    #threshold:10
    >20c0669e4c6f8436c9d42736df727c83_478
    ATCGAACTTAAATTATTTGTTTAAATTTTTAAATAGAAAAAGTTTAGTTGGGGAAACTTAAAGTAAAAGGTAACGCTTTA
    TTTTTTTGTCAGGAGCCTGTAGTATGGAAAAATGAAAAAGTTACCGTAGGGATAACAGCGCTTTCTTCTCTGAGAGGACT
    AATTAAAGAGTTGGTTGCG
    >a36d3f7291c173c4243f22c2afbd111e_49
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

    $ head -n 13 intermediate/MOL16S/SRR5534981.fasta
    #left_primer:RRWRGACRAGAAGACCCT
    #right_primer:ARTCCAACATCGAGGT
    #raw_fastq:341476
    #flash:314983
    #cutadapt:2322
    #abundance:833
    #threshold:10
    >abb4c9d82203b201ff91fc87b7c4e337_624
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
the final 12bp, ``CGTACATCTAGT``.

This is one of the control species in the mock community, and I initially
considered this might be from the mixing of the primers to give a longer
product which could be trimmed to look a SPH16S product *or* a MOL16S product.

However, on closer examination it appears to be a chimera, and if we were to
check for and require the author's bespoke adapter sequences, we could have
excluded it automatically.

The most common merged paired-read sequence for this sample contains this
subsequence of interest:

.. code:: console

    $ cat tmp_merged/SRR5534981.fasta.gz | gunzip | head -n 4
    #raw_fastq:341476
    #flash:314983
    >aa2a73352f76a966b72155d74f0f0e5c_93368
    TCCTATGTAGGGGAAGGTATGAATGGTTTGACGTGGGAAAAGCTGTCTCTTTTATATAGAAAGAAGTTTATTTTTGAGTG
    AAAAAGCTTAAATATTTGTAAAAGACGAGAAGACCCTATCGAACTTGAATTGTGTGTTTTAGTTTTGGAATACAGAAAGT
    TTAGTTGGGGAAACTTAAAGTTAAGAAAAACGCTTTTTTGTTATAAAATGATCCTGTATTATAGAAAAATGAAAAAGTTA
    CCGTAGGGATAACAGCGCTTTCTTCTCTGAGAGGACTAATCAAAGAGTTGGTTGCGACCTCGATGTTCGTACATCTAGTA
    CG

Breaking this down, the most common merged read starts ``TCCTATG`` (author's
forward adapter E), then the SPH16S_F primer (``TAGGGGAAGGTATGAATGGTTTG``),
then 277bp of sequence matching *Pisidium compressum*
(``ACGTGGGAAAAGCTGTCTCTTTTATATAGA...AATCAAAGAGTTGGTTGCGACCTCGATGTT``) and
finally ``CGTACATCTAGTACG``. Neither this nor its reverse complement are
immediately recognisable, but it is certainly not the full length reverse
primer, nor does it include the author's expected reverse adapter sequence -
and thus these reads would have been excluded in the author's analysis.

The vast majority of these reads were rejected by THAPBI PICT too. However, in
over six hundred cases, variation in those final three bases was enough for it
to be taken as a fragment of the right primer, giving an apparent amplicon
product of mostly  *Pisidium compressum* but ending ``CGTACATCTAGT``.

Minimum threshold
-----------------

Clearly using a minimum abundance threshold of 10 is too low, and it should be
increased to at least 20 based on this. However, we have the two exceptional
sequences present at over 500 copies. Setting the minimum that high seems
excessive - but perhaps the THAPBI PICT default of 100 is more reasonable?
