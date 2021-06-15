High level overview
===================

The high level summary is that all the samples have high coverage, much higher
than most of the examples we have used. There is minimal off target signal
(from the other primer sets), and the no template blanks have lower yields.
The read counts in the blanks are high, but happily do not appear to contain
nematode sequence.

Per-marker yield
----------------

We'll start by looking at the number of read-pairs found for each marker.
After calling ``./run.sh`` you should be able to inspect these report files
at the command line or in Excel.

.. code:: console

    $ cut -f 1,2,5-8 summary/NF1-18Sr2b.samples.onebp.tsv
    <SEE TABLE BELOW>

Or open the Excel version ``summary/NF1-18Sr2b.samples.onebp.xlsx``, and focus
on those early columns:

============= ====== ========= ======= ======== ==========
#marker       sample Raw FASTQ Flash   Cutadapt Read count
============= ====== ========= ======= ======== ==========
D3Af-D3Br     Blank  1193593   1039205 0        0
D3Af-D3Br     MC1    3897994   3317661 0        0
D3Af-D3Br     MC2    4228233   3685150 0        0
D3Af-D3Br     MC3    4309817   3864130 0        0
JB3-JB5GED    Blank  69641     62060   0        0
JB3-JB5GED    MC1    1236201   1157824 0        0
JB3-JB5GED    MC2    2160885   2058441 1        0
JB3-JB5GED    MC3    1204900   1139777 0        0
NF1-18Sr2b    Blank  260778    218813  187776   116145
NF1-18Sr2b    MC1    2483453   2126062 2109488  1062526
NF1-18Sr2b    MC2    2349364   1985981 1972923  1060055
NF1-18Sr2b    MC3    2435278   2088185 2070379  1108276
SSUF04-SSUR22 Blank  57199     46879   0        0
SSUF04-SSUR22 MC1    3162379   2633321 77       0
SSUF04-SSUR22 MC2    2790363   2370732 280      0
SSUF04-SSUR22 MC3    1953138   1640045 52       0
============= ====== ========= ======= ======== ==========

You should find the raw FASTQ numbers match the author's Table 5, although
that omits the blanks - which happily are all much lower.

The "Flash" column reports how many of those raw FASTQ read pairs could be
overlap merged into a single sequence - and our numbers range from 82% to 95%
(it is easy to add this calculation in Excel). This is very different from the
author's results in Table 6, although we agree that the best yield was with
the JB3-JB5GED markers. Exploring the flash settings here, using ``-O`` or
``--allow-outies`` was important here to maximize yield, but that alone does
not explain this discrepancy.

The "Cutadapt" column reports how many of those merged reads could be primer
trimmed with the NF1-18Sr2b primers, and happily we get high numbers only for
the NF1-18Sr2b samples, but low levels from the other samples. That could be
barcode leakage in the demultiplexing, or actual unwanted DNA in the samples.

The final column highlighted here is the "Read count" after applying our
minimum abundance threshold - and now we only get reads from the NF1-18Sr2b
samples.

We can repeat this for the other three primer sets, and the same pattern is
observed - strong signal only for the matching samples (with the blanks giving
strong but lower counts), and all non-matching samples zero after the minimum
abundance threshold is applied.

Blank controls
--------------

The excellent news is at the default minimum abundance threshold there are no
recognisable nematode sequences in any of the blanks.

Looking at the same sample reports (or the more detailed read reports), we
see that while the blank samples with no PCR template control give lots of
reads, where they can be identified the organisms are not seen in the mock
communities. Quoting the paper:

  *Blank samples only yielded sequences of fungi and streptophyta.*

In our case, we found lots of fungi and also the genus *Urtica* (which is a
green plant under streptophyta), but also some *Blastocystis* (Stramenopiles),
*Cercomonas* (Rhizaria) and *Sphaerularioidea* (Opisthokonta).

.. code:: console

    $ for MARKER in NF1-18Sr2b SSUF04-SSUR22 D3Af-D3Br JB3-JB5GED; do \
      grep $MARKER.Blank summary/$MARKER.samples.onebp.tsv | cut -f 1,2,4; \
      done
    <SEE TABLE EXCERPT BELOW>

Or manually looking at the four separate files - where column 4 is a text
summary of the classifier output:

============= ===== ===================================================================================
NF1-18Sr2b    Blank Fungi (unknown species), Urtica sp., Unknown
SSUF04-SSUR22 Blank Blastocystis sp., Fungi (unknown species), Unknown
D3Af-D3Br     Blank Cercomonas sp., Fungi (unknown species), Sphaerularioidea gen. sp. EM-2016, Unknown
JB3-JB5GED    Blank Unknown
============= ===== ===================================================================================

It should stressed that all the blank samples have unknown sequences (indeed
the JB3-JB5GED blank sequences are *all* reported as unknown).
