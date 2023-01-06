Presence and absence
====================

As discussed in the paper, the recognised species recovered from the mock
community varied dramatically by marker. This example has been setup with
the same list of 23 species expected for all the markers.

Note that none of the four reference sets include a known sequence for
*Laimaphelenchus penardi*, and most are missing more than just that species.

The ``run.sh`` script runs a classifier assessment over all the samples which
is meaningful for the pooled results. There is then a loop to assess each
marker individually on the four relevant samples only.

We can compare these results to Ahmed *et al.* (2019) Table 9.

NF1-18Sr2b
----------

This marker as the best database coverage, lacking only any entry for
*Laimaphelenchus penardi*, which is therefore the only expected false negative.

.. code:: console

    $ cut -f 1-5,9,11 summary/NF1-18Sr2b.assess.onebp.tsv
    <SEE TABLE BELOW>

Or open this in Excel. You should find:

================================== == == == === ==== ===========
#Species                           TP FP FN TN  F1   Ad-hoc-loss
================================== == == == === ==== ===========
OVERALL                            38 51 31 200 0.48 0.683
Acrobeles sp.                      0  0  3  1   0.00 1.000
Acrobeloides sp.                   1  0  2  1   0.50 0.667
Alaimus sp.                        0  0  3  1   0.00 1.000
Anaplectus sp.                     0  0  3  1   0.00 1.000
Anatonchus tridentatus             3  0  0  1   1.00 0.000
Aphelenchoides sp.                 0  0  3  1   0.00 1.000
Aporcelaimellus sp.                3  0  0  1   1.00 0.000
Criconema sp.                      0  0  3  1   0.00 1.000
Ditylenchus dipsaci                3  0  0  1   1.00 0.000
Ditylenchus weischeri              0  3  0  1   0.00 1.000
Globodera achilleae                0  3  0  1   0.00 1.000
Globodera artemisiae               0  3  0  1   0.00 1.000
Globodera mexicana                 0  3  0  1   0.00 1.000
Globodera pallida                  0  3  0  1   0.00 1.000
Globodera rostochiensis            3  0  0  1   1.00 0.000
Globodera sp.                      0  3  0  1   0.00 1.000
Globodera tabacum                  0  3  0  1   0.00 1.000
Hemicycliophora sp.                0  0  3  1   0.00 1.000
Laimaphelenchus penardi            0  0  3  1   0.00 1.000
Longidorus caespiticola            3  0  0  1   1.00 0.000
Meloidogyne hapla                  3  0  0  1   1.00 0.000
Plectus sp.                        3  0  0  1   1.00 0.000
Prionchulus cf. punctatus TSH-2005 0  2  0  2   0.00 1.000
Prionchulus muscorum               0  2  0  2   0.00 1.000
Prionchulus punctatus              2  0  1  1   0.80 0.333
Pristionchus sp.                   3  0  0  1   1.00 0.000
Rhabditis sp.                      3  0  0  1   1.00 0.000
Steinernema carpocapsae            3  0  0  1   1.00 0.000
Steinernema monticolum             0  3  0  1   0.00 1.000
Steinernema sp.                    0  3  0  1   0.00 1.000
Steinernema websteri               0  3  0  1   0.00 1.000
Trichodorus primitivus             3  0  0  1   1.00 0.000
Tripyla daviesae                   0  3  0  1   0.00 1.000
Tripyla glomerans                  0  0  3  1   0.00 1.000
Tripyla sp.                        0  3  0  1   0.00 1.000
Tylenchus sp.                      0  0  3  1   0.00 1.000
Urtica sp.                         0  1  0  3   0.00 1.000
Xiphinema bakeri                   0  2  0  2   0.00 1.000
Xiphinema coxi europaeum           0  2  0  2   0.00 1.000
Xiphinema diversicaudatum          2  0  1  1   0.80 0.333
Xiphinema japonicum                0  2  0  2   0.00 1.000
Xiphinema pseudocoxi               0  2  0  2   0.00 1.000
Xiphinema vuittenezi               0  2  0  2   0.00 1.000
OTHER 37 SPECIES IN DB             0  0  0  148 0.00 0.000
================================== == == == === ==== ===========

We have explainable false positives as within genus conflicts in
*Ditylenchus*, *Globodera*, *Prionchulus*, and *Xiphinema*.
Additionally there are unexplained FP from *Tripyla daviesae*, *Tripyla* sp.
and *Urtica* sp.

We also have false negatives. This is not performing as well as the authors'
analysis:

> The NF1-18Sr2b had the highest coverage, producing 100% recovery of the
> sampled taxa (Table 9). All 23 taxa were detected in all three replicates,
> apart from *Acrobeles* and *Criconema*. They both failed to appear in one
> of the replicates.

Perhaps our abundance threshold is too high?

SSUF04-SSUR22
-------------

The assess command here warns the DB lacks 10 of the expected species in
the mock community, which are therefore false negatives.

.. code:: console

    $ cut -f 1-5,9,11 summary/SSUF04-SSUR22.assess.onebp.tsv
    <SEE TABLE BELOW>

Or open this in Excel. You should find:

========================= == == == == ==== ===========
#Species                  TP FP FN TN F1   Ad-hoc-loss
========================= == == == == ==== ===========
OVERALL                   30 4  39 39 0.58 0.589
Acrobeles sp.             0  0  3  1  0.00 1.000
Acrobeloides sp.          2  0  1  1  0.80 0.333
Alaimus sp.               3  0  0  1  1.00 0.000
Anaplectus sp.            3  0  0  1  1.00 0.000
Anatonchus tridentatus    3  0  0  1  1.00 0.000
Aphelenchoides sp.        0  0  3  1  0.00 1.000
Aporcelaimellus sp.       3  0  0  1  1.00 0.000
Blastocystis sp.          0  1  0  3  0.00 1.000
Criconema sp.             0  0  3  1  0.00 1.000
Ditylenchus dipsaci       0  0  3  1  0.00 1.000
Globodera rostochiensis   0  0  3  1  0.00 1.000
Hemicycliophora sp.       0  0  3  1  0.00 1.000
Laimaphelenchus penardi   0  0  3  1  0.00 1.000
Longidorus caespiticola   3  0  0  1  1.00 0.000
Meloidogyne hapla         0  0  3  1  0.00 1.000
Plectus sp.               3  0  0  1  1.00 0.000
Prionchulus muscorum      0  2  0  2  0.00 1.000
Prionchulus punctatus     2  0  1  1  0.80 0.333
Prionchulus sp.           0  1  0  3  0.00 1.000
Pristionchus sp.          0  0  3  1  0.00 1.000
Rhabditis sp.             0  0  3  1  0.00 1.000
Steinernema carpocapsae   3  0  0  1  1.00 0.000
Trichodorus primitivus    3  0  0  1  1.00 0.000
Tripyla glomerans         0  0  3  1  0.00 1.000
Tylenchus sp.             0  0  3  1  0.00 1.000
Xiphinema diversicaudatum 2  0  1  1  0.80 0.333
OTHER 2 SPECIES IN DB     0  0  0  8  0.00 0.000
========================= == == == == ==== ===========

There are false positives within the genus *Prionchulus*, and also from
*Blastocystis* sp.

We have TP for 11 species only. The original analysis reported recovering 15
out of 23 species with this marker (Table 9), and wrote:

> In the case of the SSUF04-SSUR22 marker, eight taxa were missing from all
> three assignment methods. The taxa that were recovered occurred in all three
> replicates. With all three methods of taxonomy assignment combined, the
> number of correctly assigned OTUs improved to 56.

Many of our false negatives are likely due to the database coverage, with
the Table 9 noting the majority of their reference sequences from NCBI RefSeq
were partial - our pipeline requires full length reference amplicons.

D3Af-D3Br
---------

The assess command here warns the DB lacks three of the expected species in
the mock community, *Criconema* sp., *Laimaphelenchus penardi*, and
*Steinernema carpocapsae* - which are therefore false negatives.

.. code:: console

    $ cut -f 1-5,9,11 summary/D3Af-D3Br.assess.onebp.tsv
    <SEE TABLE BELOW>

Or open this in Excel. You should find:

================================= == == == === ==== ===========
#Species                          TP FP FN TN  F1   Ad-hoc-loss
================================= == == == === ==== ===========
OVERALL                           29 11 40 104 0.53 0.637
Acrobeles sp.                     0  0  3  1   0.00 1.000
Acrobeloides sp.                  0  0  3  1   0.00 1.000
Alaimus sp.                       0  0  3  1   0.00 1.000
Anaplectus sp.                    0  0  3  1   0.00 1.000
Anatonchus tridentatus            0  0  3  1   0.00 1.000
Aphelenchoides sp.                0  0  3  1   0.00 1.000
Aporcelaimellus sp.               3  0  0  1   1.00 0.000
Cercomonas sp.                    0  1  0  3   0.00 1.000
Criconema sp.                     0  0  3  1   0.00 1.000
Ditylenchus dipsaci               3  0  0  1   1.00 0.000
Globodera pallida                 0  3  0  1   0.00 1.000
Globodera rostochiensis           3  0  0  1   1.00 0.000
Globodera sp.                     0  3  0  1   0.00 1.000
Hemicycliophora sp.               0  0  3  1   0.00 1.000
Laimaphelenchus penardi           0  0  3  1   0.00 1.000
Longidorus caespiticola           0  0  3  1   0.00 1.000
Meloidogyne hapla                 3  0  0  1   1.00 0.000
Plectus sp.                       3  0  0  1   1.00 0.000
Prionchulus punctatus             3  0  0  1   1.00 0.000
Pristionchus sp.                  3  0  0  1   1.00 0.000
Rhabditis sp.                     3  0  0  1   1.00 0.000
Sphaerularioidea gen. sp. EM-2016 0  1  0  3   0.00 1.000
Steinernema carpocapsae           0  0  3  1   0.00 1.000
Trichodorus primitivus            3  0  0  1   1.00 0.000
Tripyla glomerans                 0  0  3  1   0.00 1.000
Tylenchus sp.                     0  0  3  1   0.00 1.000
Xiphinema bakeri                  0  1  0  3   0.00 1.000
Xiphinema diversicaudatum         2  0  1  1   0.80 0.333
Xiphinema japonicum               0  1  0  3   0.00 1.000
Xiphinema sp.                     0  1  0  3   0.00 1.000
OTHER 16 SPECIES IN DB            0  0  0  64  0.00 0.000
================================= == == == === ==== ===========

Most of the false positives are within the genus *Globodera* or *Xiphinema*,
but additionally *Cercomonas* sp. and *Sphaerularioidea* gen. sp. EM-2016.

We have 10 species correctly identified (from all three samples, save
*Xiphinema diversicaudatum* in just two samples), which exceeds authors'
analysis with UTAX but falls short of their consensus:

> The 28S rDNA-based D3Af-D3Br marker assigned 70 OTUs to nematodes and
> recovered all taxa except *Criconema* in the consensus taxonomy. Amongst the
> recovered taxa, *Hemicycliophora* occurred in one of the replicates,
> *Acrobeles* in two, while the rest were found in all three replicates.

JB3-JB5GED
----------

The assess command here warns the DB lacks 20 of the expected species in the
mock community, which puts the results into perspective:

.. code:: console

    $ cut -f 1-5,9,11 summary/JB3-JB5GED.assess.onebp.tsv
    <SEE TABLE BELOW>

Or open this in Excel. You should find:

========================= == == == == ==== ===========
#Species                  TP FP FN TN F1   Ad-hoc-loss
========================= == == == == ==== ===========
OVERALL                   3  0  66 27 0.08 0.957
Acrobeles sp.             0  0  3  1  0.00 1.000
Acrobeloides sp.          0  0  3  1  0.00 1.000
Alaimus sp.               0  0  3  1  0.00 1.000
Anaplectus sp.            0  0  3  1  0.00 1.000
Anatonchus tridentatus    0  0  3  1  0.00 1.000
Aphelenchoides sp.        0  0  3  1  0.00 1.000
Aporcelaimellus sp.       0  0  3  1  0.00 1.000
Criconema sp.             0  0  3  1  0.00 1.000
Ditylenchus dipsaci       0  0  3  1  0.00 1.000
Globodera rostochiensis   0  0  3  1  0.00 1.000
Hemicycliophora sp.       0  0  3  1  0.00 1.000
Laimaphelenchus penardi   0  0  3  1  0.00 1.000
Longidorus caespiticola   0  0  3  1  0.00 1.000
Meloidogyne hapla         3  0  0  1  1.00 0.000
Plectus sp.               0  0  3  1  0.00 1.000
Prionchulus punctatus     0  0  3  1  0.00 1.000
Pristionchus sp.          0  0  3  1  0.00 1.000
Rhabditis sp.             0  0  3  1  0.00 1.000
Steinernema carpocapsae   0  0  3  1  0.00 1.000
Trichodorus primitivus    0  0  3  1  0.00 1.000
Tripyla glomerans         0  0  3  1  0.00 1.000
Tylenchus sp.             0  0  3  1  0.00 1.000
Xiphinema diversicaudatum 0  0  3  1  0.00 1.000
OTHER 1 SPECIES IN DB     0  0  0  4  0.00 0.000
========================= == == == == ==== ===========

This has performed perfectly on *Meloidogyne hapla*, but despite having
entries in the database gives false negatives for *Globodera rostochiensis*,
and *Steinernema carpocapsae*.

There are sequences close to the DB entry for *Steinernema carpocapsae*, but
more than one base pair away so not matched by the default classifier.

This is similar to the authors analysis:

> For the COI-based JB3-JB5GED marker, even the consensus taxonomy drawn from
> all three assignment methods could only recover two taxa, namely Meloidogyne
> and Steinernema.

Pooled
------

The pipeline is setup to assess the pooled results expecting all 23 species in
each mock community, regardless of which marker was being sequenced. i.e. This
is handicapped by adding up to 9 false negatives per species.

.. code:: console

    $ cut -f 1-5,9,11 summary/pooled.assess.onebp.tsv
    <SEE TABLE BELOW>

Or open this in Excel. You should find:

================================== === == === ==== ==== ===========
#Species                           TP  FP FN  TN   F1   Ad-hoc-loss
================================== === == === ==== ==== ===========
OVERALL                            100 66 176 1162 0.45 0.708
Acrobeles sp.                      0   0  12  4    0.00 1.000
Acrobeloides sp.                   3   0  9   4    0.40 0.750
Alaimus sp.                        3   0  9   4    0.40 0.750
Anaplectus sp.                     3   0  9   4    0.40 0.750
Anatonchus tridentatus             6   0  6   4    0.67 0.500
Aphelenchoides sp.                 0   0  12  4    0.00 1.000
Aporcelaimellus sp.                9   0  3   4    0.86 0.250
Blastocystis sp.                   0   1  0   15   0.00 1.000
Cercomonas sp.                     0   1  0   15   0.00 1.000
Criconema sp.                      0   0  12  4    0.00 1.000
Ditylenchus dipsaci                6   0  6   4    0.67 0.500
Ditylenchus weischeri              0   3  0   13   0.00 1.000
Globodera achilleae                0   3  0   13   0.00 1.000
Globodera artemisiae               0   3  0   13   0.00 1.000
Globodera mexicana                 0   3  0   13   0.00 1.000
Globodera pallida                  0   6  0   10   0.00 1.000
Globodera rostochiensis            6   0  6   4    0.67 0.500
Globodera sp.                      0   6  0   10   0.00 1.000
Globodera tabacum                  0   3  0   13   0.00 1.000
Hemicycliophora sp.                0   0  12  4    0.00 1.000
Laimaphelenchus penardi            0   0  12  4    0.00 1.000
Longidorus caespiticola            6   0  6   4    0.67 0.500
Meloidogyne hapla                  9   0  3   4    0.86 0.250
Plectus sp.                        9   0  3   4    0.86 0.250
Prionchulus cf. punctatus TSH-2005 0   2  0   14   0.00 1.000
Prionchulus muscorum               0   4  0   12   0.00 1.000
Prionchulus punctatus              7   0  5   4    0.74 0.417
Prionchulus sp.                    0   1  0   15   0.00 1.000
Pristionchus sp.                   6   0  6   4    0.67 0.500
Rhabditis sp.                      6   0  6   4    0.67 0.500
Sphaerularioidea gen. sp. EM-2016  0   1  0   15   0.00 1.000
Steinernema carpocapsae            6   0  6   4    0.67 0.500
Steinernema monticolum             0   3  0   13   0.00 1.000
Steinernema sp.                    0   3  0   13   0.00 1.000
Steinernema websteri               0   3  0   13   0.00 1.000
Trichodorus primitivus             9   0  3   4    0.86 0.250
Tripyla daviesae                   0   3  0   13   0.00 1.000
Tripyla glomerans                  0   0  12  4    0.00 1.000
Tripyla sp.                        0   3  0   13   0.00 1.000
Tylenchus sp.                      0   0  12  4    0.00 1.000
Urtica sp.                         0   1  0   15   0.00 1.000
Xiphinema bakeri                   0   3  0   13   0.00 1.000
Xiphinema coxi europaeum           0   2  0   14   0.00 1.000
Xiphinema diversicaudatum          6   0  6   4    0.67 0.500
Xiphinema japonicum                0   3  0   13   0.00 1.000
Xiphinema pseudocoxi               0   2  0   14   0.00 1.000
Xiphinema sp.                      0   1  0   15   0.00 1.000
Xiphinema vuittenezi               0   2  0   14   0.00 1.000
OTHER 46 SPECIES IN DB             0   0  0   736  0.00 0.000
================================== === == === ==== ==== ===========

As expected from the per-marker results, the false positives are largely due
to species level difficulties within the genera *Globodera*, *Steinernema*,
and *Xiphinema*.

While many of the number of false negatives may be down to database coverage,
it would also be worth exploring dropping the minimum abundance threshold.
