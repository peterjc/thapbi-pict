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
COTO  1:192 10                44234                   70926            72542                 13937
COTO  1:64  10                33297                   37207            0                     120021
EPFU  1:192 10                73437                   6424             89956                 15930
EPFU  1:64  10                66160                   37024            18702                 50973
TABR  1:192 10                99987                   48159            1289                  33098
TABR  1:64  10                88233                   71607            610                   29588
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
COTO  1:192 10                45582                   71101            81249                 3707
COTO  1:64  10                50986                   46798            88085                 4656
EPFU  1:192 10                73785                   6424             92849                 12689
EPFU  1:64  10                66160                   46125            57761                 2813
TABR  1:192 10                100567                  61193            1289                  19484
TABR  1:64  10                90603                   90455            610                   8370
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
