Pooling
=======

This is a nice example to show the pooling script included with THAPBI PICT,
here pooling on the first two columns of the sample report:

.. code:: console

    $ ../../scripts/pooling.py -i summary/430_bats.COI.samples.onebp.tsv -c 1,2
    <SEE TABLE BELOW>

You can specify an output stem like ``-o pooled`` and get ``pooled.tsv`` and
matching ``pooled.xlsx`` files, but by default the plain text table is printed
to the terminal:

==== ===== ================= ======================= ================ ===================== =======
Rare Ratio Samples-sequenced Corynorhinus townsendii Eptesicus fuscus Tadarida brasiliensis Unknown
==== ===== ================= ======================= ================ ===================== =======
COTO 1:192 10                58948                   99888            82587                 19059
COTO 1:64  10                45632                   51977            0                     148446
EPFU 1:192 10                99840                   9668             103545                21191
EPFU 1:64  10                91018                   52574            21507                 65809
TABR 1:192 10                149636                  73958            1563                  52279
TABR 1:64  10                128019                  106581           773                   50833
==== ===== ================= ======================= ================ ===================== =======

As discussed earlier, where *Corynorhinus townsendii* (COTO) is the rare
species at a 1:64 ratio there is no *Tadarida brasiliensis* matched with the
initial database, but it is found with the extended database:

.. code:: console

    $ ../../scripts/pooling.py -i summary/ext_bats.COI.samples.onebp.tsv  -c 1,2
    <SEE TABLE BELOW>

Again, shown as a table:

==== ===== ================= ======================= ================ ===================== =======
Rare Ratio Samples-sequenced Corynorhinus townsendii Eptesicus fuscus Tadarida brasiliensis Unknown
==== ===== ================= ======================= ================ ===================== =======
COTO 1:192 10                61727                   100185           92815                 5755
COTO 1:64  10                70121                   68495            101333                6106
EPFU 1:192 10                100822                  9668             108264                15490
EPFU 1:64  10                91242                   68322            67690                 3654
TABR 1:192 10                154907                  98791            1563                  22175
TABR 1:64  10                133876                  140456           773                   11101
==== ===== ================= ======================= ================ ===================== =======

One of the options in this script is ``-b`` or ``--boolean`` for a yes/no
summary rather than showing the sum of the reads:

.. code:: console

    $ ../../scripts/pooling.py -i summary/ext_bats.COI.samples.onebp.tsv  -c 1,2 -b
    <SEE TABLE BELOW>

All three species (and unknowns) are found in at least one of the 10 samples
sequenced in each of the six groups:

==== ===== ================= ======================= ================ ===================== =======
Rare Ratio Samples-sequenced Corynorhinus townsendii Eptesicus fuscus Tadarida brasiliensis Unknown
==== ===== ================= ======================= ================ ===================== =======
COTO 1:192 10                Y                       Y                Y                     Y
COTO 1:64  10                Y                       Y                Y                     Y
EPFU 1:192 10                Y                       Y                Y                     Y
EPFU 1:64  10                Y                       Y                Y                     Y
TABR 1:192 10                Y                       Y                Y                     Y
TABR 1:64  10                Y                       Y                Y                     Y
==== ===== ================= ======================= ================ ===================== =======

In the Excel output the species labels are rotated 90 degrees allowing a very
compact display.
