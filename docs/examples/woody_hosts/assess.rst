.. _assess:

Assessing the classifier
========================

This sample dataset includes two positive control mock communities.
We know the species which went into the two different DNA mixes used,
so for each sequenced positive control sample we can compare the
expected list of species with the predicted list of species, and thus
count true positives, false positives, false negatives, etc.

We will first do this by hand, and then explore the tool's own built
in assessment framework.

Counting species for one sample by hand
---------------------------------------

The woody hosts dataset had two positive control mixes. From the first
plate, a set of 15 *Phytophthora* species (listed here alphabetically):

- *Phytophthora austrocedri*
- *Phytophthora boehmeriae*
- *Phytophthora cactorum*
- *Phytophthora cambivora*
- *Phytophthora chlamydospora*
- *Phytophthora cinnamomi*
- *Phytophthora gonapodyides*
- *Phytophthora ilicis*
- *Phytophthora kernoviae*
- *Phytophthora lateralis*
- *Phytophthora obscura*
- *Phytophthora plurivora*
- *Phytophthora pseudosyringae*
- *Phytophthora ramorum*
- *Phytophthora syringae*

Quoting from the :ref:`sample summary report <summary_reports>`, using the
default settings for classification of ``DNA15MIX``, we got:

.. code:: console

    $ grep DNA15MIX summary/thapbi-pict.ITS1.samples.onebp.tsv | cut -f 2
    Phytophthora aleatoria(*), Phytophthora alpina(*), Phytophthora austrocedri, Phytophthora cactorum(*), Phytophthora gonapodyides, Phytophthora ilicis, Phytophthora kernoviae, Phytophthora obscura, Phytophthora pseudosyringae, Phytophthora ramorum

Or, as a list:

- *Phytophthora aleatoria* (uncertain/ambiguous)
- *Phytophthora alpina* (uncertain/ambiguous)
- *Phytophthora austrocedri*
- *Phytophthora cactorum* (uncertain/ambiguous)
- *Phytophthora gonapodyides*
- *Phytophthora ilicis*
- *Phytophthora kernoviae*
- *Phytophthora obscura*
- *Phytophthora pseudosyringae*
- *Phytophthora ramorum*

The good news is that eight are correct classifications (eight true
positives, 8 TP), but two false positives (2 FP). Those false positives
*Phytophthora alpina* and *P. aleatoria* are indistinguishable from
*P. cactorum*, a problem flagged via the ``conflicts`` command:

.. code:: console

    $ thapbi_pict conflicts | grep cactorum
    Loaded taxonomy for 1329 sequences from DB
    f27df8e8755049e831b1ea4521ad6eb3  species  Phytophthora aleatoria;Phytophthora alpina;Phytophthora cactorum
    $ grep f27df8e8755049e831b1ea4521ad6eb3 intermediate/ITS1/DNA15MIX.fasta
    >f27df8e8755049e831b1ea4521ad6eb3_981

The bad news is we are missing seven expected species (seven false negatives,
7 FN):

- *Phytophthora boehmeriae*
- *Phytophthora cambivora*
- *Phytophthora chlamydospora*
- *Phytophthora cinnamomi*
- *Phytophthora lateralis*
- *Phytophthora plurivora*
- *Phytophthora syringae*

We will return to interpretation after showing how to get the tool to compute
these FP, FP and FN values.

Running thapbi_pict assess for one sample
-----------------------------------------

Comparing a few samples like this by hand is one thing, but doing it at scale
requires automation. For assessing changes to the classifier method and
database, we mainly run ``thapbi_pict assess`` against a set of single isolate
positive controls. This requires a computer readable files listing the
expected species in a particular format.

.. code:: console

    $ thapbi_pict assess -h
    ...

The simplest way to run the assess command is to tell it two TSV input
filenames, named ``<sample_name>.known.tsv`` (the expected results) and
``<sample_name>.<method>.tsv`` (from running `thapbi_pict classify`` on
``<sample_name>.fasta``).

The "known" file uses the same column based layout as the intermediate TSV
files, but while you can provide the expected species for each unique sequence
in the sample, this can be simplified to a single wildcard ``*`` line
followed by all the NCBI taxids (optional), and species names using semi-colon
separators.

The assess command will default to printing its tabular output to screen -
shown here abridged after piping through the ``cut`` command to pull out just
the first five columns:

.. code:: console

    $ rm -rf DNA15MIX.onebp.tsv
    $ thapbi_pict classify -i intermediate/ITS1/DNA15MIX.fasta -o .
    $ thapbi_pict assess -i expected/DNA15MIX.known.tsv DNA15MIX.onebp.tsv | cut -f 1-5
    Assessed onebp vs known in 1 files (230 species)
    #Species                     TP  FP  FN  TN
    OVERALL                      8   2   7   213
    Phytophthora aleatoria       0   1   0   0
    Phytophthora alpina          0   1   0   0
    Phytophthora austrocedri     1   0   0   0
    Phytophthora boehmeriae      0   0   1   0
    Phytophthora cactorum        1   0   0   0
    Phytophthora cambivora       0   0   1   0
    Phytophthora chlamydospora   0   0   1   0
    Phytophthora cinnamomi       0   0   1   0
    Phytophthora gonapodyides    1   0   0   0
    Phytophthora ilicis          1   0   0   0
    Phytophthora kernoviae       1   0   0   0
    Phytophthora lateralis       0   0   1   0
    Phytophthora obscura         1   0   0   0
    Phytophthora plurivora       0   0   1   0
    Phytophthora pseudosyringae  1   0   0   0
    Phytophthora ramorum         1   0   0   0
    Phytophthora syringae        0   0   1   0
    OTHER 213 SPECIES IN DB      0   0   0   213

More usually, you would output to a named file, and look at that:

.. code:: console

    $ thapbi_pict assess -i expected/DNA15MIX.known.tsv DNA15MIX.onebp.tsv \
      -o DNA15MIX.assess.tsv
    Assessed onebp vs known in 1 files (230 species)
    $ cut -f 1-5,9,11 DNA15MIX.assess.tsv
    <SEE TABLE BELOW>

You should be able to open this ``DNA15MIX.assess.tsv`` file in R, Excel, etc,
and focus on the same column selection:

=========================== == == == === ==== ===========
#Species                    TP FP FN TN  F1   Ad-hoc-loss
=========================== == == == === ==== ===========
OVERALL                     8  2  7  213 0.64 0.529
Phytophthora aleatoria      0  1  0  0   0.00 1.000
Phytophthora alpina         0  1  0  0   0.00 1.000
Phytophthora austrocedri    1  0  0  0   1.00 0.000
Phytophthora boehmeriae     0  0  1  0   0.00 1.000
Phytophthora cactorum       1  0  0  0   1.00 0.000
Phytophthora cambivora      0  0  1  0   0.00 1.000
Phytophthora chlamydospora  0  0  1  0   0.00 1.000
Phytophthora cinnamomi      0  0  1  0   0.00 1.000
Phytophthora gonapodyides   1  0  0  0   1.00 0.000
Phytophthora ilicis         1  0  0  0   1.00 0.000
Phytophthora kernoviae      1  0  0  0   1.00 0.000
Phytophthora lateralis      0  0  1  0   0.00 1.000
Phytophthora obscura        1  0  0  0   1.00 0.000
Phytophthora plurivora      0  0  1  0   0.00 1.000
Phytophthora pseudosyringae 1  0  0  0   1.00 0.000
Phytophthora ramorum        1  0  0  0   1.00 0.000
Phytophthora syringae       0  0  1  0   0.00 1.000
OTHER 213 SPECIES IN DB     0  0  0  213 0.00 0.000
=========================== == == == === ==== ===========

The ``OVERALL`` line tells us that there were 8 true positives, 2 false
positives, 7 false negatives, and 213 true negatives. The final number needs a
little explanation. First, 8+2+7+213 = 230, which is the number of species in
the database. With only one sample being considered, 213 is the number of
other species in the database which the tool correctly says are not present.

Following this we get one line per species, considering this species in
isolation (making this a traditional and simpler to interpret classification
problem). Here there is only one sample, so this time TP+FP+FN+TN=1.

You can easily spot the 2 FP in this layout, *Phytophthora alpina* and
*P. aleatoria*, or the 7 FN.

The additional columns (not all shown here) include traditional metrics like
sensitivity, specificity, precision, F1, and Hamming loss. We've shown F1 or
F-measure here (from zero to one for perfect recall), plus our own metric
provisionally called *Ad hoc loss* which is a modification of the Hamming loss
without using the true negative count (which we expect to always be very large
as the database will contain many species, while a community might contain
only ten).

Positive Control - 10 species mix
---------------------------------

The positive controls from the second plate had a different mix of ten
*Phytophthora* species, again listed alphabetically:

- *Phytophthora boehmeriae*
- *Phytophthora cactorum*
- *Phytophthora capsici*
- *Phytophthora castaneae*
- *Phytophthora fallax*
- *Phytophthora foliorum*
- *Phytophthora obscura*
- *Phytophthora plurivora*
- *Phytophthora rubi*
- *Phytophthora siskiyouensis*

Again referring to the sample summary report from running with default settings,
for ``DNA10MIX_undiluted`` and ``DNA10MIX_diluted25x`` we got:

- *Phytophthora agathidicida* (uncertain/ambiguous)
- *Phytophthora capsici*
- *Phytophthora castaneae* (uncertain/ambiguous)
- *Phytophthora fallax*
- *Phytophthora foliorum*
- *Phytophthora glovera* (uncertain/ambiguous)
- *Phytophthora obscura*
- *Phytophthora plurivora*
- *Phytophthora rubi*
- *Phytophthora siskiyouensis*

Plus the results from ``DNA10MIX_bycopynumber`` were almost the same - but this
time there wasn't a sequence only matched to *P. capsici*, giving:

- *Phytophthora agathidicida* (uncertain/ambiguous)
- *Phytophthora capsici* (uncertain/ambiguous)
- *Phytophthora castaneae* (uncertain/ambiguous)
- *Phytophthora fallax*
- *Phytophthora foliorum*
- *Phytophthora glovera* (uncertain/ambiguous)
- *Phytophthora obscura*
- *Phytophthora plurivora*
- *Phytophthora rubi*
- *Phytophthora siskiyouensis*

The exact preparation of the 10 species mixture (with and without dilution,
etc) made little difference.

Leaving aside the ambiguous qualifier, there are ten species predictions, but
only nine are correct (9 TP: *P. capsici*, *P. castaneae*, *P. fallax*,
*P. foliorum*, *P. obscura*, *P. plurivora*, *P. rubi*, *P. siskiyouensis*),
with two wrong guesses (2 FP: *P. agathidicida* and *P. glovera*), and two
missing predictions (2 FN: *P. boehmeriae* and *P. cactorum*).

The uncertain/ambiguous prediction of *Phytophthora agathidicida* is easily
explained, it comes from a sequence present in all three samples with MD5
checksum ``5122dde24762f8e3d6a54e3f79077254``, and this exact sequence is in
the database with entries for both *Phytophthora castaneae* (which was in the
DNA control mixture) and also *Phytophthora agathidicida* (e.g. accession
KP295308).

You can confirm this by looking at the intermediate TSV files, e.g. using
grep to show all lines with this species name:

.. code:: console

    $ grep "Phytophthora agathidicida" summary/thapbi-pict.ITS1.all_reads.onebp.tsv
    29de890989becddc5e0b10ecbbc11b1a_1524  1642459;1642465  Phytophthora agathidicida;Phytophthora castaneae
    $ grep 29de890989becddc5e0b10ecbbc11b1a intermediate/ITS1/*.fasta
    intermediate/ITS1/DNA10MIX_bycopynumber.fasta:>29de890989becddc5e0b10ecbbc11b1a_245
    intermediate/ITS1/DNA10MIX_diluted25x.fasta:>29de890989becddc5e0b10ecbbc11b1a_655
    intermediate/ITS1/DNA10MIX_undiluted.fasta:>29de890989becddc5e0b10ecbbc11b1a_624
    $ thapbi_pict conflicts | grep 29de890989becddc5e0b10ecbbc11b1a
    Loaded taxonomy for 1329 sequences from DB
    29de890989becddc5e0b10ecbbc11b1a  species  Phytophthora agathidicida;Phytophthora castaneae

The same applies to *Phytophthora capsici* and *Phytophthora glovera*.

As noted above, the woody hosts paper concluded the failure to detect
*P. boehmeriae* in either DNA mix was due to inefficient primer annealing
in a species mixture. We have an unexpected FN for *P. cactorum* though.

Assessing multiple samples
--------------------------

When we ran the assess command earlier on a sample sample, we provided pairs
of ``<sample_name>.<method>.tsv`` and ``<sample_name>.known.tsv``, but you
can instead provide the intermediate ``<sample_name>.fasta`` files(s) and the
pooled classifier output.

.. code:: console

    $ thapbi_pict assess -i expected/DNA10MIX_*.known.tsv intermediate/ITS1/DNA10MIX_*.fasta summary/thapbi-pict.ITS1.all_reads.onebp.tsv -o DNA10MIX.assess.tsv
    Assessed onebp vs known in 3 files (230 species)
    $ cut -f 1-5,9,11 DNA10MIX.assess.tsv
    <SEE TABLE BELOW>

New table ``DNA10MIX.assess.tsv`` is similar to what we had before for the
single 15 species sample, but here with three samples the species rows add up
to 3:

========================== == == == === ==== ===========
#Species                   TP FP FN TN  F1   Ad-hoc-loss
========================== == == == === ==== ===========
OVERALL                    24 6  6  654 0.80 0.333
Phytophthora agathidicida  0  3  0  0   0.00 1.000
Phytophthora boehmeriae    0  0  3  0   0.00 1.000
Phytophthora cactorum      0  0  3  0   0.00 1.000
Phytophthora capsici       3  0  0  0   1.00 0.000
Phytophthora castaneae     3  0  0  0   1.00 0.000
Phytophthora fallax        3  0  0  0   1.00 0.000
Phytophthora foliorum      3  0  0  0   1.00 0.000
Phytophthora glovera       0  3  0  0   0.00 1.000
Phytophthora obscura       3  0  0  0   1.00 0.000
Phytophthora plurivora     3  0  0  0   1.00 0.000
Phytophthora rubi          3  0  0  0   1.00 0.000
Phytophthora siskiyouensis 3  0  0  0   1.00 0.000
OTHER 218 SPECIES IN DB    0  0  0  654 0.00 0.000
========================== == == == === ==== ===========

We have 3 FP for *Phytophthora agathidicida* (indistinguishable from
*P. castaneae*), 3 FP for *Phytophthora glovera* (indistinguishable from
*P. capsici*), making 6 FP (two per sample).

We have 3 FN for *Phytophthora boehmeriae* (primer mismatching as before), and
and 3 FN for *Phytophthora cactorum*.

Next, let's run the assess command on all four positive control samples, just
by giving the input directory names (it will work out the common filenames,
just like doing this via the pipeline command - see below):

.. code:: console

    $ thapbi_pict assess -i expected/ intermediate/ITS1/ \
      summary/thapbi-pict.ITS1.all_reads.onebp.tsv -o thabpi-pict.ITS1.assess.tsv
    Assessed onebp vs known in 4 files (230 species)
    $ cut -f 1-5,9,11 thabpi-pict.ITS1.assess.tsv
    <SEE TABLE BELOW>

New table ``thabpi-pict.ITS1.assess.tsv`` is similar:

=========================== == == == === ==== ===========
#Species                    TP FP FN TN  F1   Ad-hoc-loss
=========================== == == == === ==== ===========
OVERALL                     32 8  13 867 0.75 0.396
Phytophthora agathidicida   0  3  0  1   0.00 1.000
Phytophthora aleatoria      0  1  0  3   0.00 1.000
Phytophthora alpina         0  1  0  3   0.00 1.000
Phytophthora austrocedri    1  0  0  3   1.00 0.000
Phytophthora boehmeriae     0  0  4  0   0.00 1.000
Phytophthora cactorum       1  0  3  0   0.40 0.750
Phytophthora cambivora      0  0  1  3   0.00 1.000
Phytophthora capsici        3  0  0  1   1.00 0.000
Phytophthora castaneae      3  0  0  1   1.00 0.000
Phytophthora chlamydospora  0  0  1  3   0.00 1.000
Phytophthora cinnamomi      0  0  1  3   0.00 1.000
Phytophthora fallax         3  0  0  1   1.00 0.000
Phytophthora foliorum       3  0  0  1   1.00 0.000
Phytophthora glovera        0  3  0  1   0.00 1.000
Phytophthora gonapodyides   1  0  0  3   1.00 0.000
Phytophthora ilicis         1  0  0  3   1.00 0.000
Phytophthora kernoviae      1  0  0  3   1.00 0.000
Phytophthora lateralis      0  0  1  3   0.00 1.000
Phytophthora obscura        4  0  0  0   1.00 0.000
Phytophthora plurivora      3  0  1  0   0.86 0.250
Phytophthora pseudosyringae 1  0  0  3   1.00 0.000
Phytophthora ramorum        1  0  0  3   1.00 0.000
Phytophthora rubi           3  0  0  1   1.00 0.000
Phytophthora siskiyouensis  3  0  0  1   1.00 0.000
Phytophthora syringae       0  0  1  3   0.00 1.000
OTHER 205 SPECIES IN DB     0  0  0  820 0.00 0.000
=========================== == == == === ==== ===========

This time the ``OVERALL`` line says we had 32 TP, 8 FP, 13 FN and 867 TN.
Their total, 32+8+13+867 = 920 = 4 * 230, is the number of samples times the
number of species in the database.

This time notice all the per-species lines have TP+FP+FN+TN=4 as there were 4
samples.

Running assessment as part of pipeline
--------------------------------------

Provided they follow the expected naming convention, if you include your
control files ``*.known.tsv`` as one of the pipeline inputs, it will call the
classifier assessment after running the classifier and producing the main
reports:

.. code:: console

    $ thapbi_pict pipeline -i raw_data/ expected/ -s intermediate/ \
      -o summary/with-metadata -n raw_data/NEGATIVE*.fastq.gz \
      -t metadata.tsv -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -x 16
    ...
    $ ls -1 summary/with-metadata.*
    summary/with-metadata.ITS1.all_reads.fasta
    summary/with-metadata.ITS1.all_reads.onebp.tsv
    summary/with-metadata.ITS1.assess.confusion.onebp.tsv
    summary/with-metadata.ITS1.assess.onebp.tsv
    summary/with-metadata.ITS1.assess.tally.onebp.tsv
    summary/with-metadata.ITS1.reads.onebp.tsv
    summary/with-metadata.ITS1.reads.onebp.xlsx
    summary/with-metadata.ITS1.samples.onebp.tsv
    summary/with-metadata.ITS1.samples.onebp.txt
    summary/with-metadata.ITS1.samples.onebp.xlsx
    $ diff summary/with-metadata.ITS1.assess.onebp.tsv \
      thabpi-pict.ITS1.assess.tsv

Output file ``summary/with-metadata.ITS1.assess.onebp.tsv`` will match the
output above.

Interpretation of the mock communities
--------------------------------------

Running our pipeline with the default settings results in a number of false
positives (all unavoidable as they come from conflicting marker sequences in
the database, see the ``thapbi_pict conflicts`` command), and some false
negatives (on top of the explained absence of *Phytophthora boehmeriae*).
Specifically we have 6 unexplained false negatives on the 15 species mix, and
are missing *Phytophthora cactorum* in all three samples of the 10 species mix.

This means that with the default settings THAPBI PICT gives a more cautious
set of predictions than the ``metapy`` tool used in the original data analysis
(see `Riddell et al. (2019) Table 1, Table 2
<https://doi.org/10.7717/peerj.6931/table-1>`_) which appears to consider even
singletons.

Attempting to compare the results in their Table 1 with our own numbers is
complicated since it appears to show just one of the 10 species mixes (so the
TP count is out of 10) while we used all three (for a TP count out of 30).

We can recover most of the missing species (the FN) by dropping the minimum
abundance thresholds (which requires deleting the intermediate FASTA files,
and re-running with lower settings for ``-a`` and ``-f``), at the cost of more
FP.

For instance, we find traces of *P. syringae* with less than 10 reads in
the 15 species mix (consistent with Table 2), and even *P. boehmeriae* with
less than 10 reads in two of the 10 species mix (not reported in Table 2).

Interestingly even excluding only singletons (using ``-a 2 -f 0``), we didn't
find any matches to *Phytophthora cactorum* in the three samples of the 10
species mix. However, there is a sequence perfectly matching database entries
for *P. idaei* present at around 40 to 60 copies, and in light of the original
paper, this is likely what was intended to be in the mixture as *P. cactorum*.

Again even excluding only singletons, we didn't find any matches to
*P. plurivora* in the 15 species mix (Table 2 in the original paper suggests
present with only 2 reads).

Optimising on maximising the F1 score and minimising ad-hoc-loss, the
optimal absolute abundance threshold is in the region of 65 to 80 reads,
giving 34 TP, 9 FP, and 11 FN for an F1 of 0.77 and ad-hoc-loss of 0.370.
