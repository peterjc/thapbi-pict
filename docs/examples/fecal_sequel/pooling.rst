Pooling
=======

This is a nice example to show the pooling script included with THAPBI PICT,
here pooling on the first two columns of the sample report:

.. code:: console

    $ ../../scripts/pooling.py -i summary/mock-community.COI_430_bats.samples.onebp.tsv -c 1,2
    <SEE TABLE BELOW>

You can specify an output stem like ``-o pooled`` and get ``pooled.tsv`` and
matching ``pooled.xlsx`` files, but by default the plain text table is printed
to the terminal:

===== ===== ================= ======================= ================ ===================== =======
#Rare Ratio Samples-sequenced Corynorhinus townsendii Eptesicus fuscus Tadarida brasiliensis Unknown
===== ===== ================= ======================= ================ ===================== =======
COTO  1:192 10                58936                   100360           82575                 19392
COTO  1:64  10                45619                   51949            0                     148505
EPFU  1:192 10                99837                   9666             103539                21111
EPFU  1:64  10                90994                   52559            21503                 65797
TABR  1:192 10                149564                  73931            1563                  52292
TABR  1:64  10                127945                  106498           773                   50740
===== ===== ================= ======================= ================ ===================== =======

As discussed earlier, where *Corynorhinus townsendii* (COTO) is the rare
species at a 1:64 ratio there is no *Tadarida brasiliensis* matched with the
initial database, but it is found with the extended database:

.. code:: console

    $ ../../scripts/pooling.py -i summary/mock-community.COI_ext_bats.samples.onebp.tsv  -c 1,2
    <SEE TABLE BELOW>

Again, shown as a table:

===== ===== ================= ======================= ================ ===================== =======
#Rare Ratio Samples-sequenced Corynorhinus townsendii Eptesicus fuscus Tadarida brasiliensis Unknown
===== ===== ================= ======================= ================ ===================== =======
COTO  1:192 10                61715                   100877           92913                 5758
COTO  1:64  10                70102                   68460            101303                6208
EPFU  1:192 10                100819                  9666             108257                15411
EPFU  1:64  10                91219                   68302            67677                 3655
TABR  1:192 10                154832                  98754            1563                  22201
TABR  1:64  10                133698                  140355           773                   11130
===== ===== ================= ======================= ================ ===================== =======

One of the options in this script is ``-b`` or ``--boolean`` for a yes/no
summary rather than showing the sum of the reads:

.. code:: console

    $ ../../scripts/pooling.py -i summary/mock-community.COI_ext_bats.samples.onebp.tsv  -c 1,2 -b
    <SEE TABLE BELOW>

All three species (and unknowns) are found in at least one of the 10 samples
sequenced in each of the six groups:

===== ===== ================= ======================= ================ ===================== =======
#Rare Ratio Samples-sequenced Corynorhinus townsendii Eptesicus fuscus Tadarida brasiliensis Unknown
===== ===== ================= ======================= ================ ===================== =======
COTO  1:192 10                Y                       Y                Y                     Y
COTO  1:64  10                Y                       Y                Y                     Y
EPFU  1:192 10                Y                       Y                Y                     Y
EPFU  1:64  10                Y                       Y                Y                     Y
TABR  1:192 10                Y                       Y                Y                     Y
TABR  1:64  10                Y                       Y                Y                     Y
===== ===== ================= ======================= ================ ===================== =======

In the Excel output the species labels are rotated 90 degrees allowing a very
compact display.
