Minimum Abundance Threshold
===========================

THAPBI PICT has a default minimum abundance threshold of 100 reads per sample
before accepting any unique sequence. Background contamination and PCR noise
levels will vary, so having multiple :ref:`negative_controls` and positive
controls will help set this objectively.

In this dataset there are three positive controls (mock communities) and two
negative controls for the ITS1 sequencing:

.. code:: console

    $ grep Control metadata.tsv | cut -f 1-2,8
    CP (Control Phytophthora)  CP1  SRR13393837
    CP (Control Phytophthora)  CP2  SRR13393813
    CP (Control Phytophthora)  CP3  SRR13393802
    CN (Control Negative)      CN1  SRR13393836
    CN (Control Negative)      CN2  SRR13393824

Alternatively, search ``PRJNA690943.tsv`` for "Phytophthora" and "Negative".

Negative Controls
-----------------

We can deliberately try running ``thapbi_pict prepare-reads`` on the two
negative control samples with lower and lower minimum abundance thresholds. It
turns out we do see ITS1 reads, but only two copies of even the most abundant:

.. code:: console

    $ rm -rf /tmp/ITS1/SRR13393836.fasta
    $ thapbi_pict prepare-reads -i raw_data/SRR13393836_* -a 1 -o /tmp/
    $ head -n 9 /tmp/ITS1/SRR13393836.fasta
    #left_primer:GAAGGTGAAGTCGTAACAAGG
    #right_primer:GCARRGACTTTCGTCCCYRC
    #raw_fastq:10310
    #flash:59
    #cutadapt:14
    #abundance:14
    #threshold:1
    >688ca1f52e3ffeb2417adceff1011655_2
    TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAACCCTTTATTG
    GGGGCTTCTGTCTGGTCTGGCTTCGGCTGGGCTGGGTGGCGGCTCTATCATGGCGACCGCTCTGGGCTTCGGCCTGGAGT
    TAGTAGCTCACTTTTAAACCCATTCTTAATTACTGAACATACT

Doing the same with the second negative control ``SRR13393824`` shows only
singletons. These negative controls are exceptionally clean, so on this basis
alone you might consider setting a very relaxed threshold, perhaps as low as
``-a 10``?

Positive Controls
-----------------

In this dataset we have three positive controls, mock communities of ten
species. We can try different abundance thresholds - too high and we'd expect
false negatives (real signal being discarded), but too low and we'd get false
positives (unwanted noise or contamination).

After running ``setup.py`` there should be 3 files ``expected/*.known.tsv``
setup with the 10 species mock community as expected classifier output. We
can use these with the ``thapbi_pict assess`` and/or ``pipeline`` command.

The provided ``run.sh`` script settled on ``-a 50`` to balance these needs.

.. code:: console

    $ ./run.sh
    ...
    $ cut -f 1-5 summary/british_soil.ITS1.assess.onebp.tsv
    <SEE TABLE BELOW>

As a table,

========================== == == == ===
#Species                   TP FP FN TN
========================== == == == ===
OVERALL                    26 7  4  958
Phytophthora agathidicida  0  3  0  2
Phytophthora boehmeriae    0  0  3  2
Phytophthora capsici       3  0  0  2
Phytophthora castaneae     3  0  0  2
Phytophthora fallax        3  0  0  2
Phytophthora foliorum      3  0  0  2
Phytophthora glovera       0  3  0  2
Phytophthora idaei         2  0  1  2
Phytophthora obscura       3  0  0  2
Phytophthora plurivora     3  0  0  2
Phytophthora rubi          3  0  0  2
Phytophthora siskiyouensis 3  0  0  2
Phytophthora syringae      0  1  0  4
OTHER 186 SPECIES IN DB    0  0  0  930
========================== == == == ===

False positives
~~~~~~~~~~~~~~~

Looking at the reports, we do see unwanted species in the controls. First we
can easily dismiss the unwanted predictions of *Phytophthora agathidicida* as
being indistinguishable from positive control *P. castaneae* (this is explicit
in the read report). Next, we have the same problem with the unwanted
*Phytophthora glovera* being indistinguishable from positive control
*P. capsici*.

There are other more interesting false positives (FP), the most prominent and
only case passing our chosen minimum abundance threshold of 50 is
*Phytophthora syringae* in ``SRR13393813`` (Control Phytophthora 2), readily
identifable in ``summary/british_soil.ITS1.samples.onebp.xlsx`` and which from
``summary/british_soil.ITS1.reads.onebp.xlsx`` is seen to be from 86 copies of
a single unique sequence ``32159de6cbb6df37d084e31c37c30e7b``:

.. code:: console

    $ grep -A 1 32159de6cbb6df37d084e31c37c30e7b intermediate/ITS1/SRR13393813.fasta
    >32159de6cbb6df37d084e31c37c30e7b_86
    TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAACCCTTTTATT
    GGGGGCTTCTGTCTGGTCTGGCTTCGGCTGGATTGGGTGGCGGCTCTATCATGGCGACCGCTCTGAGCTTCGGCCTGGAG
    CTAGTAGCCCACTTTTTAAACCCATTCTTAATTACTGAACAAACT
    $ grep 32159de6cbb6df37d084e31c37c30e7b summary/british_soil.ITS1.all_reads.onebp.tsv
    32159de6cbb6df37d084e31c37c30e7b_48894  67594  Phytophthora syringae

Interestingly before removing the primers this sequence came from a range of
unique sequences, none seen more than ten times:

.. code:: console

    $ cat tmp_merged/SRR13393813.fasta.gz | gunzip \
      | grep -B 1 "TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAACCCTTTTATTGGGGGCTTCTGTCTGGTCTGGCTTCGGCTGGATTGGGTGGCGGCTCTATCATGGCGACCGCTCTGAGCTTCGGCCTGGAGCTAGTAGCCCACTTTTTAAACCCATTCTTAATTACTGAACAAACT" \
      | grep "^>" | head
    >590e14c00cacf04bc580415ad7cca33f_10
    >85570853bb0a4f6ff59a2dc0cf1535e6_10
    >695cd92b7e12e48c0250b4eee7c6c3a1_10
    >34fc5250c9913e7fa250cbb1dd17ecb9_9
    >c6d536bd0a2b169c02564d944aa52ca3_9
    >a2de7264ea9e0fb4ffb3429253961852_6
    >ade1678898b72a487635a4d0ee729dab_5
    >519ddc9abe8ca622664debcc91b4ade0_3
    >64e9b619d178752bf8b7f7e2ec0b2998_3
    >7314c74a9e7d5e8b6bca1297473838d0_2

There is variation at the allowed ambiguities in the right primer, which we
can again show with grep. Here ``GAAGGTGAAGTCGTAACAAGG`` is the left primer,
and ``GYRGGGACGAAAGTCYYTGC`` is the reverse complement of the right primer. We
are using a regular expression wildcard in place of the ambiguous bases:

.. code:: console

    $ cat tmp_merged/SRR13393813.fasta.gz | gunzip \
      | grep -B 1 "GAAGGTGAAGTCGTAACAAGG${SEQ}G..GGGACGAAAGTC..TGC" \
      | grep "^>" | head
    >590e14c00cacf04bc580415ad7cca33f_10
    >85570853bb0a4f6ff59a2dc0cf1535e6_10
    >695cd92b7e12e48c0250b4eee7c6c3a1_10
    >34fc5250c9913e7fa250cbb1dd17ecb9_9
    >c6d536bd0a2b169c02564d944aa52ca3_9
    >a2de7264ea9e0fb4ffb3429253961852_6
    >ade1678898b72a487635a4d0ee729dab_5
    >519ddc9abe8ca622664debcc91b4ade0_3
    >64e9b619d178752bf8b7f7e2ec0b2998_3
    >7314c74a9e7d5e8b6bca1297473838d0_2

The THAPBI PICT pipeline drops low abundance sequences *after* removing the
primers when pooling unique sequences. It seems possible the authors did not
find this *Phytophthora syringae* matching sequence because their pipeline
removed these sequences as being low abundance *prior* to primer trimming?

False negatives
~~~~~~~~~~~~~~~

We find *Phytophthora boehmeriae* is absent at this minimum threshold of 50.
To double check the less abundant sequences you may wish to try running this
again specifically on just the three positive controls:

.. code:: console

    $ rm -rf controls && mkdir controls
    $ thapbi_pict pipeline -i raw_data/SRR13393802_* raw_data/SRR13393813_* \
          raw_data/SRR13393837_* expected/ -o controls/controls-only -a 1
    ...
    $ grep -E "(predictions|boehmeriae)" controls/controls-only.ITS1.reads.onebp.tsv | cut -f 1,2,5
    #Marker-MD5                       onebp-predictions        Max-sample-abundance
    7ac50609279c89c7fc3d88ffed426dac  Phytophthora boehmeriae  1
    869fb51182270e82dc07e19401f2f8c0  Phytophthora boehmeriae  1

Looking in the Excel file ``controls/controls-only.reads.onebp.xlsx`` or at
the tabular file, we find *Phytophthora boehmeriae* is practically absent,
appearing at trace level only (single reads). This false negative (FN) matches
the authors' results and was observed in Riddell *et al.* (2019) - see our
:ref:`wooded hosts worked example <woody_hosts>` - and put down to a poor
primer match for this species in competitive PCR.

.. code:: console

    $ grep -E "(predictions|idaei)" controls/controls-only.ITS1.reads.onebp.tsv | cut -f 1,2,7- | head
    #Marker-MD5                       onebp-predictions   SRR13393802  SRR13393813  SRR13393837
    fe1bd3a42e730f95c9fde798e32f8478  Phytophthora idaei  135          71           41
    23529b55e483660b4aa4b61d49002695  Phytophthora idaei  1            3            2
    5ecb44ee3586c65fdb758f0e25a92bef  Phytophthora idaei  1            2            0
    993e56c425b8651e8871fe63b33a640e  Phytophthora idaei  2            1            0
    c9a456bd908038ec7d94f75fe69d7b2a  Phytophthora idaei  2            1            0
    f43547ee36b8fbcbce121235caeec266  Phytophthora idaei  1            1            0
    6b6ac3a5d175beed741750ee044ce374  Phytophthora idaei  4            0            0
    fda75c109fad4f0878d4ad445244cde5  Phytophthora idaei  4            0            0
    122b5f2fedd9653ce0d2174f8adf7db9  Phytophthora idaei  3            0            0

The next lowest abundant species in the 3 controls and potential false
negative (FN) is *Phytophthora idaei*, with the most abundant sequence
appearing at 135, 71 and 41 copies. That might suggest setting the threshold
up to 40 reads to ensure all the positive control sequences (bar *Phytophthora
boehmeriae*) come through, while setting it higher (e.g. 50) gives at least
one false negative (FN).

Minimum threshold
-----------------

The provided ``run.sh`` script settled on ``-a 50`` to balance the desire to
minimize false positives (which suggest using a threshold over 86 to exclude
the *P. syringae*) and minimize false negatives (which suggest using a
threshold under 41 to include the *P. idaei*).
