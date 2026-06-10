.. _custom_database_pipeline:

Pipeline with custom database
=============================

Running thapbi-pict pipeline
----------------------------

Having created the dual-marker database ``pooled.sqlite``, running the pipeline
is quite straightforward - but I have modified a lot of the default settings:

.. code:: console

    $ mkdir -p intermediate/ summary/
    $ thapbi_pict pipeline -d pooled.sqlite --synthetic '' \
        -m 1s3g --denoise unoise-l --unoise_alpha 6 \
        -i raw_data/ expected/ \
        --merged-cache tmp_merged/ \
        -s intermediate/ -o summary/ \
        -t metadata.tsv -x 1 -c 3,7,4,6
    ...
    Running 1s3g classifier on summary/rps10.tally.tsv, 186 sequences
    1s3g classifier assigned species/genus to 62 of 186 unique sequences from 1 files
    ...

In addition to setting the database ``-d pooled.sqlite`` we also use
``--synthetic ''`` tell the pipeline this dataset lacks any synthetic controls.
Minimal denoising has been enabled with ``--denoise unoise-l --unoise_alpha 6``
to remove the long tail of unclassified unique amplicon sequence variants.

The classifier method was set explicitly with ``-m 1s3g`` to allow fuzzier
matching than the more conservative default reflecting the relatively sparse
rps10 database coverage.

Results
-------

Quoting Foster *et al.* (2021):

    The taxonomic classifications of ASVs produced by the rps10 method included
    23 of the 24 mock community species whereas the ITS1 classifications included
    17 (Table 3). The species missing in the rps10 classifications was
    *Phytophthora ipomoeae* and the species missing in the ITS1 classifications
    were *P. citrophthora*, *P. himalsilva, *P. ipomoeae*, *P. quercina*,
    *Pythium dissotocum*, *P. oligandrum*, and *P. undulatum.*

Here are the results for the ITS1 marker on the 24 species mock community,
sequenced twice so the true positive (TP) and false negative (FN) columns
sum to 48:

.. code:: console

    $ cut -f 1-4 summary/ITS1.assess.1s3g.tsv
    <SEE TABLE BELOW>

We have lots of false negatives with ITS1, often for both samples:

=============================== == == ==
#Species                        TP FP FN
=============================== == == ==
OVERALL                         10 6  38
Aphanomyces euteiches           0  0  2
Elongisporangium undulatum      0  0  2
Globisporangium apiculatum      0  0  2
Globisporangium irregulare      0  0  2
Peronospora farinosa            0  0  2
Peronospora schachtii           0  0  2
Phytophthora andina             0  1  0
Phytophthora castanetorum       0  1  0
Phytophthora cinnamomi          1  0  1
Phytophthora citrophthora       1  0  1
Phytophthora himalsilva         0  0  2
Phytophthora hydrogena          1  0  1
Phytophthora infestans          1  0  1
Phytophthora ipomoeae           1  0  1
Phytophthora kernoviae          1  0  1
Phytophthora pluvialis          1  0  1
Phytophthora pseudocitrophthora 0  1  0
Phytophthora quercina           1  0  1
Phytophthora ramorum            1  0  1
Phytophthora rosacearum         1  0  1
Phytophthora terminalis         0  1  0
Phytophthora urerae             0  1  0
Phytophthora versiformis        0  1  0
Phytopythium citrinum           0  0  2
Plasmopara halstedii            0  0  2
Plasmopara obducens             0  0  2
Pseudoperonospora cubensis      0  0  2
Pythium dissotocum              0  0  2
Pythium oligandrum              0  0  2
Saprolegnia diclina             0  0  2
OTHER 303 SPECIES IN DB         0  0  0
=============================== == == ==

This isn't actually a shock if you spotted this warning while running the pipeline:

    WARNING: 13 expected species were not a possible prediction:
    Aphanomyces euteiches;Elongisporangium undulatum;Globisporangium irregulare;Peronospora farinosa;Peronospora schachtii;Phytopythium citrinum;Plasmopara halstedii;Plasmopara obducens;Pseudoperonospora cubensis;Pythium apiculatum;Pythium dissotocum;Pythium oligandrum;Saprolegnia diclina

We can ignore the 13 non-*Phytophthora* and the associated 26 false negatives
because the THAPBI PICT ITS1 database intentionally doesn't cover them.

In fact, if you pooled the two samples (which I assume was done for the quote
above), then all the *Phytophthora* are found, other than *Phytophthora himalsilva*.
So for the *Phytophthora* at least, this ITS1 analysis is more successful.

.. code:: console

    $ cut -f 1-4 summary/rps10.assess.1s3g.tsv
    <SEE TABLE BELOW>

We have less false negatives with rps10:

========================== == == ==
#Species                   TP FP FN
========================== == == ==
OVERALL                    17 2  31
Aphanomyces euteiches      1  0  1
Elongisporangium undulatum 0  0  2
Globisporangium apiculatum 1  0  1
Globisporangium irregulare 1  0  1
Peronospora farinosa       0  0  2
Peronospora schachtii      0  0  2
Phytophthora andina        0  1  0
Phytophthora betacei       0  1  0
Phytophthora cinnamomi     1  0  1
Phytophthora citrophthora  1  0  1
Phytophthora himalsilva    0  0  2
Phytophthora hydrogena     1  0  1
Phytophthora infestans     1  0  1
Phytophthora ipomoeae      1  0  1
Phytophthora kernoviae     1  0  1
Phytophthora pluvialis     1  0  1
Phytophthora quercina      1  0  1
Phytophthora ramorum       1  0  1
Phytophthora rosacearum    1  0  1
Phytopythium citrinum      1  0  1
Plasmopara halstedii       0  0  2
Plasmopara obducens        1  0  1
Pseudoperonospora cubensis 0  0  2
Pythium dissotocum         0  0  2
Pythium oligandrum         1  0  1
Saprolegnia diclina        1  0  1
OTHER 317 SPECIES IN DB    0  0  0
========================== == == ==

This time the database is only missing one species, the warning was:

    WARNING: 1 expected species were not a possible prediction:
    Peronospora farinosa

This is because the version of the OomyceteDB rps10 reference sequences used
lacks *Peronospora farinosa*.

As to the rest, the THAPBI PICT default abundance thresholds are cautious and
set to minimise false positives with MiSeq data volumes. It is likely that
lowering these would report the missing entries of *Elongisporangium undulatum*,
*Peronospora schachtii*, *Phytophthora himalsilva*, *Plasmopara halstedii*,
*Pseudoperonospora cubensis*, and *Pseudoperonospora cubensis* (at the cost of
more false positives).

Note however that unlike the original analysis, this already finds *Phytophthora
ipomoeae* in one of the mock community samples.
