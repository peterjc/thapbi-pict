High level overview
===================

The high level summary is that all the samples have high coverage, much higher
than most of the examples we have used. Some of the samples yield over a
million reads for the COI and 12S amplicons, which with the default fractional
minimum abundance threshold of 0.1% (``-f 0.001``) would mean using over 1000
reads as the threshold. This was too stringent, so the worked example reduces
this to 0.01% (with ``-f 0.0001``) matching the author's analysis, and dropped
the default absolute abundance threshold of 100 to 50 (with ``-a 50``).

Note that the rarest members of the mock communities are expected from 1 in
500 individuals (0.02%) or 1 in 1000 individuals (0.01%), which is ten times
higher than the fractional abundance threshold.

Sequence yield
--------------

We'll start by looking at the number of read-pairs found for each marker.
After calling ``./run.sh`` you should be able to inspect these report files
at the command line or in Excel.

.. code:: console

    $ cut -f 3,6-8,10,12-14 summary/COI.samples.onebp.tsv
    <SEE TABLE BELOW>

Or open the Excel version ``summary/COI.samples.onebp.xlsx``, and focus
on those early columns:

============ ========= ======= ======== ========= ========== ======== ======
sample_alias Raw FASTQ Flash   Cutadapt Threshold Singletons Accepted Unique
============ ========= ======= ======== ========= ========== ======== ======
100-Pool-1   478705    474621  109233   50        11074      86402    3
250-Pool-1   1845819   1829913 157310   50        23119      118383   3
500-Pool-1   647776    643030  51092    50        6446       36718    3
1000-Pool-1  855997    848914  66002    50        7967       49058    3
100-Pool-2   737998    732014  432826   50        29168      368249   4
250-Pool-2   2037475   2022814 1250718  126       85718      1042562  4
500-Pool-2   1908370   1895715 1231908  124       59702      1042441  4
1000-Pool-2  1068715   1060596 584017   59        33060      498955   4
100-Pool-3   950692    940342  249422   50        24964      189156   7
250-Pool-3   1631700   1615113 274422   50        39974      192944   7
500-Pool-3   923807    916621  358429   50        32221      284819   6
1000-Pool-3  1773647   1758637 468361   50        42263      374487   7
100-Pool-4   634017    628523  117499   50        14596      74799    7
250-Pool-4   2501145   2480381 441558   50        61904      324512   6
500-Pool-4   572779    568565  144488   50        18279      96537    7
1000-Pool-4  1198812   1189853 294607   50        30130      220678   7
100-Pool-5   1817929   1800594 434739   50        45224      329015   5
250-Pool-5   1632786   1617219 440995   50        58159      328842   5
500-Pool-5   807060    801471  321428   50        30944      247519   5
1000-Pool-5  1423279   1411512 332286   50        32751      255309   5
Trap-1       1759819   1719671 110882   50        19740      73024    7
Trap-10      2445993   2420303 308371   50        58670      204842   8
Trap-2       1127739   1107970 110856   50        24385      55757    7
Trap-3       2422054   2366037 161686   50        30631      110043   8
Trap-4       742893    732907  63107    50        11933      35225    5
Trap-5       3437292   3346620 346696   50        71464      208989   10
Trap-6       697389    689125  91284    50        17153      57037    7
Trap-7       2853448   2820200 223330   50        31011      169121   6
Trap-8       2196646   2161966 146646   50        28814      92632    9
Trap-9       2065455   2049024 70591    50        14636      40131    7
============ ========= ======= ======== ========= ========== ======== ======

The marker specific tables show the threshold applied was usually 50, the
default absolute value set via ``-a 50`` at the command line. Occasionally
this has been increased to 0.1% of the sequences matching the primers for this
marker, set via ``-f 0.0001`` at the command line.

The numbers are similar for the 12S and 18S markers, or pooling them all:

.. code:: console

    $ cut -f 3,6,7,9,10 summary/pooled.samples.onebp.tsv
    <SEE TABLE BELOW>

Again, alternatively open Excel file ``summary/pooled.samples.onebp.xlsx``,
and focus on those early columns:

============ ========= ======= ======== ======
sample_alias Raw FASTQ Flash   Accepted Unique
============ ========= ======= ======== ======
100-Pool-1   478705    474621  371045   5
250-Pool-1   1845819   1829913 1508292  5
500-Pool-1   647776    643030  522396   5
1000-Pool-1  855997    848914  692639   5
100-Pool-2   737998    732014  587902   5
250-Pool-2   2037475   2022814 1243165  5
500-Pool-2   1908370   1895715 1551757  7
1000-Pool-2  1068715   1060596 863574   5
100-Pool-3   950692    940342  684297   10
250-Pool-3   1631700   1615113 1158575  9
500-Pool-3   923807    916621  697552   8
1000-Pool-3  1773647   1758637 1366298  9
100-Pool-4   634017    628523  451801   10
250-Pool-4   2501145   2480381 1867605  8
500-Pool-4   572779    568565  416456   9
1000-Pool-4  1198812   1189853 918004   10
100-Pool-5   1817929   1800594 1369274  8
250-Pool-5   1632786   1617219 1128901  7
500-Pool-5   807060    801471  603390   7
1000-Pool-5  1423279   1411512 1104412  8
Trap-1       1759819   1719671 392775   13
Trap-10      2445993   2420303 492325   19
Trap-2       1127739   1107970 129956   13
Trap-3       2422054   2366037 427533   14
Trap-4       742893    732907  232800   10
Trap-5       3437292   3346620 486282   18
Trap-6       697389    689125  80003    13
Trap-7       2853448   2820200 1158684  9
Trap-8       2196646   2161966 683669   19
Trap-9       2065455   2049024 1352408  10
============ ========= ======= ======== ======

The "Accepted" column is the number of reads matching the primer pairs and
passing our abundance thresholds. The fraction accepted varies from 61% to
82% for the mock community samples, but is considerably lower for the
environmental traps, varying from 11% to 65%. Much of that would be noise and
trace level environmental DNA.

The "Unique" column is the number of accepted unique sequences. For the mock
communities this should be up to 18 with at most six species each, and three
markers. The observed counts are lower, so perhaps the abundance threshold is
still too harsh? On the other hand, dropping it further does raise the false
positive rate inferred from the mock communities.
