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

    $ rm -rf /tmp/SRR13393836.fasta
    $ thapbi_pict prepare-reads -i raw_data/SRR13393836_* -a 1 -o /tmp/
    $ head -n 9 /tmp/SRR13393836.fasta
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
    $ cut -f 1-5 summary/british_soil_its1.assess.onebp.tsv
    <SEE TABLE BELOW>

As a table,

========================== == == == ===
#Species                   TP FP FN TN
========================== == == == ===
OVERALL                    26 7  4  833
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
OTHER 161 SPECIES IN DB    0  0  0  805
========================== == == == ===

Looking at the reports, we do see unwanted species in the controls. First we
can easily dismiss the unwanted predictions of *Phytophthora agathidicida* as
being indistinguishable from positive control *P. castaneae* (this is explicit
in the read report). Next, we have the same problem with the unwanted
*Phytophthora glovera* being indistinguishable from positive control
*P. capsici*.

There are other more interesting false positives (FP), the most prominent and
only case passing our chosen minimum abundance threshold of 50 is
*Phytophthora syringae* in ``SRR13393813`` (Control Phytophthora 2), with the
most abundance sequence appearing 86 copies.

.. code:: console

    $ grep "Phytophthora syringae" intermediate/SRR13393813.onebp.tsv
    32159de6cbb6df37d084e31c37c30e7b_86  67594  Phytophthora syringae

We find *Phytophthora boehmeriae* is absent at this minimum threshold of 50.
To double check the less abundant sequences you may wish to try running this
again specifically on just the three positive controls:

.. code:: console

    $ rm -rf controls && mkdir controls
    $ thapbi_pict pipeline -i raw_data/SRR13393802_* raw_data/SRR13393813_* \
          raw_data/SRR13393837_* expected/ -o controls/ -r controls-only -a 1
    ...
    onebp classifier assigned species/genus to 16467 of 21259 sequences from 3 files
    ...
    $ grep -E "(predictions|boehmeriae)" controls/controls-only.reads.onebp.tsv | cut -f 1,2,5
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

    $ grep -E "(predictions|idaei)" controls/controls-only.reads.onebp.tsv | cut -f 1,2,7- | head
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
