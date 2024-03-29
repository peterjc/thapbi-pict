Presence and absence
====================

Controls
--------

Quoting the Muri *et al.* (2020) paper:

    A low-frequency noise threshold of 0.001 (0.1%) was applied across the
    dataset to reduce the probability of false positives arising from
    cross-contamination or tag-jumping (De Barba et al. 2014; Hänfling et al.
    2016). Based on the level of contamination found in sampling/filtration
    blanks and PCR negatives, a second arbitrary threshold was applied and all
    records occurring with less than 50 reads assigned were removed.

To match the paper, this example uses ``-a 50`` for an absolute threshold of
50, and ``-f 0.001`` for a 0.1% sample specific factional threshold.

At this threshold, the 4 cichlid "positive" samples, 6 PCR "negative", and 8
"blank" controls are perfect - as far as the fish go. We do see unexpected
human and chicken reads in the PCR negatives, and also ducks, cattle and pigs
in the field "blanks":

.. code:: console

    $ grep -E "(^#|positive|negative|blank)" summary/drained_ponds.12S.samples.onebp.tsv | cut -f 5,10-11,16,19
    <SEE TABLE BELOW>

Or, filter/search ``summary/drained_ponds.12S.samples.onebp.tsv`` in Excel:

======== ================= =========================================================================== ========= ========
control  Sequencing sample Classification summary                                                      Threshold Accepted
======== ================= =========================================================================== ========= ========
blank    SRR11949861       -                                                                           50        0
blank    SRR11949885       -                                                                           50        0
blank    SRR11949884       (Off-target) Homo sapiens, (Off-target) Sus scrofa                          50        544
blank    SRR11949883       (Off-target) Bos taurus, (Off-target) Homo sapiens, (Off-target) Sus scrofa 50        1629
blank    SRR11949882       (Off-target) Anatidae (waterfowl)                                           50        61
blank    SRR11949881       (Off-target) Homo sapiens                                                   50        56
blank    SRR11949880       (Off-target) Anatidae (waterfowl), (Off-target) Homo sapiens                50        436
blank    SRR11949834       (Off-target) Homo sapiens                                                   50        175
negative SRR11949908       -                                                                           50        0
negative SRR11949907       (Off-target) Gallus gallus, (Off-target) Homo sapiens                       50        606
negative SRR11949851       -                                                                           50        0
negative SRR11949850       -                                                                           50        0
negative SRR11949838       (Off-target) Homo sapiens                                                   50        71
negative SRR11949837       (Off-target) Homo sapiens                                                   50        356
positive SRR11949836       Astatotilapia calliptera(*), Maylandia zebra(*)                             50        39748
positive SRR11949835       Astatotilapia calliptera(*), Maylandia zebra(*)                             50        39244
positive SRR11949906       Astatotilapia calliptera(*), Maylandia zebra(*)                             65        62249
positive SRR11949849       Astatotilapia calliptera(*), Maylandia zebra(*)                             50        24566
======== ================= =========================================================================== ========= ========

Only in one sample (``SRR11949906``, a positive control) was the percentage
based abundance threshold stricter than the absolute threshold (65 not 50),
and it still gives the highest number of reads.

Note that the positive samples only yield a single unique sequence (MD5
checksum ``17dbc1c331d17cd075aabd6f710a039b``) which matches both the cichlid
control species *Astatotilapia calliptera* and *Maylandia zebra*.

High level overview
-------------------

Looking over ``summary/drained_ponds.12S.samples.onebp.xlsx`` in Excel, or the
TSV equivalent of the sample report, there are some general trends visible.

First, as noted above the controls are extremely clean with just the expected
cichlid species for the positive controls, and a few off-target matches in the
PCR negatives and field blanks as noted above.

Next, both the Middle Lake Sterivex Ethanol Buffer, and Middle Lake Sterivex
RNAlater Buffer are very clean. There are traces of human and waterfowl, and
but only three of these buffer samples shows any fish (eg ``SRR11949911`` aka
``5RNB``). In contrast, most of the Middle Lake Sterivex Longmire Buffer
samples do give fish reads.

Controls and buffers aside, all the field samples gave Anatidae (waterfowl)
matches, most had human. There were traces of other birds and mammals such as
pig and dog - with most of the new lake samples showing sheep (*Ovis*).

As to the fish, we see strong signal in most samples for *Abramis brama*,
*Carassius carassius*, *Cyprinus carpio*, *Rutilus rutilus* and *Tinca tinca*.

Expected Fish
-------------

This paper was selected as an example because something is known about the
expected content of all the biological samples - the lakes were drained and
all the fish identified, counted and weighed. However, we cannot expect all
the species present to have left DNA at all the sampling points within their
lake, but that is a useful approximation for assessing the classifier:

.. code:: console

    $ cut -f 1-5,9,11 summary/drained_ponds.12S.assess.onebp.tsv
    <SEE TABLE BELOW>

You might prefer to open this in Excel:

=================================== === === === ==== ==== ===========
#Species                            TP  FP  FN  TN   F1   Ad-hoc-loss
=================================== === === === ==== ==== ===========
OVERALL                             433 388 331 5877 0.55 0.624
(Off-target) Anatidae (waterfowl)   0   70  0   29   0.00 1.000
(Off-target) Apodemus               0   4   0   95   0.00 1.000
(Off-target) Ardea cinerea          0   11  0   88   0.00 1.000
(Off-target) Bos taurus             0   3   0   96   0.00 1.000
(Off-target) Canis lupus familiaris 0   7   0   92   0.00 1.000
(Off-target) Capra hircus           0   1   0   98   0.00 1.000
(Off-target) Columba                0   47  0   52   0.00 1.000
(Off-target) Gallinula chloropus    0   50  0   49   0.00 1.000
(Off-target) Gallus gallus          0   13  0   86   0.00 1.000
(Off-target) Homo sapiens           0   83  0   16   0.00 1.000
(Off-target) Ovis aries             0   17  0   82   0.00 1.000
(Off-target) Ovis dalli             0   1   0   98   0.00 1.000
(Off-target) Phalacrocorax carbo    0   25  0   74   0.00 1.000
(Off-target) Sturnus                0   3   0   96   0.00 1.000
(Off-target) Sus scrofa             0   16  0   83   0.00 1.000
(Off-target) Turdus                 0   7   0   92   0.00 1.000
Abramis brama                       65  0   16  18   0.89 0.198
Acipenser spp.                      0   0   9   90   0.00 1.000
Alburnus mossulensis                0   1   0   98   0.00 1.000
Astatotilapia calliptera            4   0   0   95   1.00 0.000
Barbus barbus                       46  0   35  18   0.72 0.432
Carassius carassius                 64  0   17  18   0.88 0.210
Ctenopharyngodon idella             3   15  6   75   0.22 0.875
Cyprinus carpio                     61  0   20  18   0.86 0.247
Maylandia zebra                     4   0   0   95   1.00 0.000
Perca fluviatilis                   40  0   41  18   0.66 0.506
Pseudorasbora parva                 0   2   0   97   0.00 1.000
Rutilus rutilus                     63  0   18  18   0.88 0.222
Scardinius erythrophthalmus         6   0   75  18   0.14 0.926
Silurus glanis                      9   0   0   90   1.00 0.000
Spinibarbus denticulatus            0   11  0   88   0.00 1.000
Squalidus gracilis                  0   1   0   98   0.00 1.000
Squalius cephalus                   6   0   75  18   0.14 0.926
Tinca tinca                         62  0   19  18   0.87 0.235
OTHER 37 SPECIES IN DB              0   0   0   3663 0.00 0.000
=================================== === === === ==== ==== ===========

False positives
---------------

We touched on the assorted "false positives" from the off-target 12S PCR
amplification above. What is more interesting is the fish false positives.
Let's look at these starting with the most false positives.

*Ctenopharyngodon idella*
~~~~~~~~~~~~~~~~~~~~~~~~~

First, many middle lake samples unexpectedly have *Ctenopharyngodon idella*
(this is expected in the new lake samples). Why? They all stem from sequence
``285edce3d193c92b1959e60bc130b518`` which was matched to both *C. idella*
and *Tinca tinca* (expected in both lakes)::

    >285edce3d193c92b1959e60bc130b518
    ACTATGCTCAGCCATAAACCTAGACATCCACCTACAATTAAACGTCCGCCCGGGTACTACGAGCATTAGCTTGAAACCCA
    AAGGACCTGACGGTGCCTTAGACCCCC

This is both a one base pair edit away from AY897013.1 etc as *C. idella*, and
from AB218686.1 etc as *T. tinca*. Reviewing the NCBI BLAST matches both sets
of species are supported from multiple complete mitochondrion genomes and a
range of research groups. In the context of this experiment, we could infer
for the four middle lake samples this sequence was *T. tinca*.

*Spinibarbus denticulatus*
~~~~~~~~~~~~~~~~~~~~~~~~~~

Next, we see 16 samples with unexpected cyprinid fish *Spinibarbus
denticulatus*. Referring to the read report, all are from a single sequence
``4c53f6ed1ecdad3af2299999ec83d756`` which has been matched perfectly to both
this unexpected species and expected species *Carassius carassius*::

    >4c53f6ed1ecdad3af2299999ec83d756
    ACTATGCTCAGCCGTAAACTTAGACATCCTACTACAATAGATGTCCGCCAGGGTACTACGAGCATTAGCTTAAAACCCAA
    AGGACCTGACGGTGTCTCAGACCCCC

Given the actual fish in these lakes have been taxonomically identified, we
can safely dismiss this - and perhaps drop AP013335.1 *S. denticulatus* from
the ad-hoc DB?

A similar choice was made in compiling the *ad hoc* database, dropping all the
*Sander* sp. entries for the following sequence in favour of just *Perca
fluviatilis* as the sole expected Percidae::

    >7e88b1bdeff6b6a361cc2175f4f630fd
    ACTATGCCTAGCCATAAACATTGGTAGCACACTACACCCACTACCCGCCTGGGAACTACGAGCATCAGCTTGAAACCCAA
    AGGACTTGGCGGTGCTTTAGATCCAC

This was based on the authors' choice:

    All fish OTUs were identified to species level with the exceptions of
    records matching the family Percidae. Percidae records were manually
    assigned to *P. fluviatilis* as this was the only species of the family
    identified in the study area during fish relocation.

*Pseudorasbora parva*
~~~~~~~~~~~~~~~~~~~~~

We see two samples containing *Pseudorasbora parva*, the invasive species
which prompted these fish ponds to be drained as a control measure. You can
find this in the read report, at the command line:

.. code:: console

    $ grep -E "(Pseudorasbora parva|samples|predictions)" \
      summary/drained_ponds.12S.reads.onebp.tsv | cut -f 2,3,7,48,59
                                                           samples          2LMB         3LMF
    MD5                               onebp-predictions    Total-abundance  SRR11949854  SRR11949925
    e819f3c222d6493572534fb6a5b7cda7  Pseudorasbora parva  520              323          197

Specifically we saw 323 reads in ``SRR11949854`` aka ``2LMB`` and 197 reads in
``SRR11949925`` aka ``3LMF`` - both middle lake Sterivex (STX) samples.
Quoting the paper:

    *P. parva* reads found in two Middle Lake-STX samples (279 and 148 reads)
    were also excluded from further analyses as after eradication this species
    was not physically present at the site surveyed.

The exact counts differ, but referring to the paper's supplementary data the
sample names match.

Other Fish
~~~~~~~~~~

We also see one false positive for each of the two fish species *Alburnus
mossulensis*, and *Squalidus gracilis*:

.. code:: console

    $ grep -E "(Alburnus mossulensis|samples|predictions)" \
      summary/drained_ponds.12S.reads.onebp.tsv | cut -f 2,3,7,25
                                                                          samples          M3-MF2
    MD5                               onebp-predictions                   Total-abundance  SRR11949859
    916da937dccfd5d29502e83713e5d998  Abramis brama;Alburnus mossulensis  98               98

This sequence is ambiguous with equally good matches to expected species
*Abramis brama*. Again, we might remove *Alburnus mossulensis* from the DB?

.. code:: console

    $ grep -E "(Squalidus gracilis|samples|predictions)" \
      summary/drained_ponds.12S.reads.onebp.tsv | cut -f 2,3,7,20
                                                          samples          M3-4F2
    MD5                               onebp-predictions   Total-abundance  SRR11949871
    c0d532d1c6f8ffff9c72ac4a1873151c  Squalidus gracilis  82               82

This sequence match is with AP011393.1 in the provided reference set.

False negatives
---------------

The classifier assessment shown above expected all the fish in each lake to be
found at all the sites within that lake - an overly strong assertion which
could explain many of the reported false negatives.

However, there is one clear false negative - neither this nor the original
analysis found any *Acipenser* spp.

True positives
--------------

Rather than reviewing all of the true positives, I will note that in some
cases we found more reads and thus declared a result in more samples.
For example, we report *Barbus barbus* in 49 samples, versus:

    In addition, *Barbus barbus* was detected at two sites (202 reads), ...

We found *Scardinius erythrophthalmus* in six samples:

.. code:: console

    $ grep -E "(Scardinius erythrophthalmus|samples|predictions)" \
      summary/drained_ponds.12S.reads.onebp.tsv | cut -f 7,8,12,13,83,84,85
    samples          M3-1F1       M3-5F1       M3-6F1       7RNF         8RNF         MRNF
    Total-abundance  SRR11949879  SRR11949870  SRR11949868  SRR11949893  SRR11949886  SRR11949852
    761              156          120          147          136          76           126

Quoting the original paper:

    The presence of *Scardinius erythrophthalmus* was found at two sites with
    a low number of reads (38 and 25 reads) and, therefore, removed after
    applying the filter threshold

In these cases at least, we are seeing much higher read counts. Given the
supplementary data provided, it could be possible to plot the read counts from
the two methods against each other.

Conclusion
----------

While not in-depth, this hopefully demonstrates the THAPBI PICT could be
meaningfully applied to this 12S dataset which was originally analysed with
metaBEAT v0.97.11.
