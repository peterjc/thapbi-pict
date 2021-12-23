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
script provided uses ``-a 10 -f 0`` to accept any unique sequence seen in
sample at least ten times (regardless the fraction of the sample read total).
This *does* allow unwanted noise though to the reports.

SPH16S
------

This was the more specific primer pair, expected to only amplify sphaeriid
mussel species, so in general we expect less unique sequences than with the
more general MOL16S primers.

Looking at some key columns in the sample report,

.. code:: console

    $ cut -f 1,2,4,7,9-10 summary/SPH16S.samples.onebp.tsv
    <SEE TABLE BELOW>

Or, open ``SPH16S.samples.onebp.xlsx`` in Excel. Focusing on the left hand
columns, you should see:

======= ============== ============ ========= ======== ==========
#Marker Group          Library-name Raw FASTQ Cutadapt Read count
======= ============== ============ ========= ======== ==========
MOL16S  Aquarium       BIR2M        306311    2        0
MOL16S  Aquarium       BIR6M        291954    14       0
MOL16S  Control        BIM8M        2433      0        0
MOL16S  Mock Community SC3PRO1      689712    17       0
MOL16S  Mock Community SC3PRO2      405048    0        0
MOL16S  Mock Community SC3PRO3      402219    16       0
MOL16S  Mock Community NFSC3PRO3    349590    33       10
MOL16S  Mock Community SC3PRO4      671241    6        0
MOL16S  Mock Community NFSC3PRO4    420015    7        0
MOL16S  Mock Community SC3PRO5      480606    13       0
MOL16S  River          BIM6M        821849    0        0
MOL16S  River          BIM2M        1119271   0        0
MOL16S  River          BIM4M        709472    40       19
SPH16S  Aquarium       BIR2S        498926    251148   209401
SPH16S  Aquarium       BIR6S        240360    226012   191446
SPH16S  Mock Community SPSC3PRO1    425271    317961   224690
SPH16S  Mock Community SPSC3PRO2    341476    282516   204248
SPH16S  Mock Community SPSC3PRO4    410780    303957   197505
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

Or, open ``MOL16S.samples.onebp.xlsx`` in Excel. Focusing on the left hand
columns, you should see:

======= ============== ============ ========= ======== ==========
#Marker Group          Library-name Raw FASTQ Cutadapt Read count
======= ============== ============ ========= ======== ==========
MOL16S  Aquarium       BIR2M        306311    297656   256386
MOL16S  Aquarium       BIR6M        291954    286427   256471
MOL16S  Control        BIM8M        2433      1010     551
MOL16S  Mock Community SC3PRO1      689712    656660   550295
MOL16S  Mock Community SC3PRO2      405048    377024   297877
MOL16S  Mock Community SC3PRO3      402219    380345   304629
MOL16S  Mock Community NFSC3PRO3    349590    328952   262964
MOL16S  Mock Community SC3PRO4      671241    628640   494270
MOL16S  Mock Community NFSC3PRO4    420015    364229   262727
MOL16S  Mock Community SC3PRO5      480606    458894   383868
MOL16S  River          BIM6M        821849    799333   703594
MOL16S  River          BIM2M        1119271   954520   823795
MOL16S  River          BIM4M        709472    367468   317366
SPH16S  Aquarium       BIR2S        498926    25       0
SPH16S  Aquarium       BIR6S        240360    27       0
SPH16S  Mock Community SPSC3PRO1    425271    35       0
SPH16S  Mock Community SPSC3PRO2    341476    168      27
SPH16S  Mock Community SPSC3PRO4    410780    420      108
======= ============== ============ ========= ======== ==========

Looking at the same points, I see two problems:

* The control sample BIM8M (SRR5534986) had almost a thousand unwanted MOL16S
  matches, reduced to 551 with a minimum abundance threshold of 10.

* All the SPH16S mock community samples have unwanted MOS16S matches, the
  worst case being SPSC3PRO4 (SRR5534980) with over four hundred reads reduced
  to 108 with the minimum abundance threshold of 10.

To see exactly what is in these two problematic samples, we can turn to the
read report - or look directly at the intermediate FASTA files:

.. code:: console

    $ cat intermediate/MOL16S/SRR5534986.fasta
    #marker:MOL16S
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
    >778e3dace4b993135e11d450e6c559ff_11
    ATCGAACTTAAATTATTTGTTTAAATTTTAAATAGAAAAAGTTTAGTTGGGGAAACTTAAAGTAAAAGGTAACGCTTTAT
    TTTTTTGTCAGGAGCCTGTAGTATGGAAAAATGAAAAAGTTACCGTAGGGATAACAGCGCTTTCTTCTCTGAGAGGACTA
    ATTAAAGAGTTGGTTGCG

The unwanted sequences in the control sample are dominated by a single
sequence (and variants of it; shown line wrapped at 80 characters), which was
matched to *Sphaerium simile*:

.. code:: console

    $ grep 20c0669e4c6f8436c9d42736df727c83 summary/MOL16S.reads.onebp.tsv | cut -f 1,2
    20c0669e4c6f8436c9d42736df727c83  Sphaerium simile

Or look at this TSV or ``summary/MOL16S.reads.onebp.xlsx`` in Excel.

This is consistent with the original author's analysis - although our pipeline
has produced higher read counts:

    Finally, our water blank sample had 71 reads, eight of those being
    singletons with the remaining belonging to *Sphaerium striatinum*
    (Table 9), likely due to amplicon contamination in the lab.

What about the other problematic sample? Again, you can find this in the Excel
read report, or at the command line:

.. code:: console

    $ cat intermediate/MOL16S/SRR5534980.fasta
    #marker:MOL16S
    #left_primer:RRWRGACRAGAAGACCCT
    #right_primer:ARTCCAACATCGAGGT
    #raw_fastq:410780
    #flash:375539
    #cutadapt:420
    #abundance:108
    #threshold:10
    >20c0669e4c6f8436c9d42736df727c83_46
    ATCGAACTTAAATTATTTGTTTAAATTTTTAAATAGAAAAAGTTTAGTTGGGGAAACTTAAAGTAAAAGGTAACGCTTTA
    TTTTTTTGTCAGGAGCCTGTAGTATGGAAAAATGAAAAAGTTACCGTAGGGATAACAGCGCTTTCTTCTCTGAGAGGACT
    AATTAAAGAGTTGGTTGCG
    >ecdaa082b70f5e268f76128693531760_45
    ATCGAACTTAATCATTTTTAAAAGAACTCATTTAAAAAGAATTTTTACTGGGGCAGTAAGAAGAAAAAAATAACTCTTCC
    TTATAAAAAAAAGATCCCTCAATGAGGAGAAAAAGAAAAAGTTACCGTAGGGATAACAGCGTTATCGTTTTTAAGAGATC
    TAATCGAAGAAACGGTTTGCG
    >98dc259e48de3e258cb93a34c38a9484_17
    ATCGAACTTAATAGTTTTTAAGAGAAATAGCTTAGAAAGAAGTTTTACTGGGGCAGTAAGAAGAAAAAAATAATTCTTCC
    TTGAAAAAAAGATCCCTTATTAAGGACAAAAGAAAAAGTTACCGTAGGGATAACAGCGTTATCGTTTTTAAGAGAACTAA
    TCGAAGAAACGGTTTGCG
    $ grep -E "(MD5|20c0669e4c6f8436c9d42736df727c83|ecdaa082b70f5e268f76128693531760|98dc259e48de3e258cb93a34c38a9484)" \
      summary/MOL16S.reads.onebp.tsv | cut -f 1,2
    <SEE TABLE BELOW>

Giving:

================================ =========================================
#Marker-MD5                      onebp-predictions
================================ =========================================
ecdaa082b70f5e268f76128693531760 Dreissena bugensis;Dreissena rostriformis
98dc259e48de3e258cb93a34c38a9484 Dreissena polymorpha
20c0669e4c6f8436c9d42736df727c83 Sphaerium simile
================================ =========================================

The unwanted mock community sample content is split between *Sphaerium* and
*Dreissena*, and suggest using a minimum threshold of perhaps 50 reads?

Minimum threshold
-----------------

Clearly using a minimum abundance threshold of 10 is too low, and it should be
increased to at perhaps 50 based on this. However, we have one exceptional
sequence present at almost 500 copies. Setting the minimum that high seems
excessive - but perhaps the THAPBI PICT default of 100 is more reasonable?
