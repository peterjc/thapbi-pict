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

Positive control - 15 species mix
---------------------------------

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

- *Phytophthora aleatoria* (uncertain/ambiguous)
- *Phytophthora austrocedri*
- *Phytophthora cactorum* (uncertain/ambiguous)
- *Phytophthora gonapodyides*
- *Phytophthora ilicis*
- *Phytophthora kernoviae*
- *Phytophthora obscura*
- *Phytophthora pseudosyringae*
- *Phytophthora ramorum*

The good news is that eight are correct classifications (eight true
positives, 8 TP, and one false positive, 1 FP). That false positive
*Phytophthora aleatoria* was indistinguishable from *P. cactorum*
(a similar example is discussed below in more detail).

The bad news is we are missing seven expected species (seven false
negatives, 7 FN):

- *Phytophthora boehmeriae*
- *Phytophthora cambivora*
- *Phytophthora chlamydospora*
- *Phytophthora cinnamomi*
- *Phytophthora lateralis*
- *Phytophthora plurivora*
- *Phytophthora syringae*

The woody hosts paper discusses the failure to detect *P. boehmeriae* in
either this or the ten species DNA mix run on the second plate, and concluded
that differences in the primer binding sites could have led to inefficient
primer annealing and thus limited amplification in a species mix.

That still leaves six unexplained false negatives, meaning with the default
settings THAPBI PICT has given a more cautious set of predictions than the
``metapy`` tool used in the original data analysis (see `Riddel *et al.* (2019)
Table 1 <https://doi.org/10.7717/peerj.6931/table-1>`_).

Dropping the minimum read abundance threshold from the default of 100 to 50
(which requires deleting the intermediate FASTA files, and re-running with
``-a 50``), adds three more true positives (*Phytophthora chlamydospora*,
*P. cinnamomi*, and *P. lateralis*) without any false negatives. That gives us
11 TP, 0 FP, 4 FN, and compares favourably with the original ``metapy``
analysis.

Dropping the minimum abundance further, we also get *P. cambivora* (plus a
false positive sharing the same ITS1 sequence as *Phytophthora x cambivora*,
discussed below). Eventually, with less than 10 reads each, we find sequences
matching *P. syringae* and *P. plurivora* (and thus all the expected species).
However, this also brings in false positives for *P. idaei* and
*P. mississippiae*. We would not trust sequences this poorly supported.


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

Again refering to the sample summary report from running with default settings,
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
with one wrong guess (1 FP: *P. agathidicida*), and one missing prediction
(1 FN: *P. boehmeriae*).

As noted above, the woody hosts paper concluded the failure to detect
*P. boehmeriae* in either DNA mix was due to inefficient primer annealing
in a species mixture.

The uncertain/ambiguous prediction of *Phytophthora agathidicida* is easily
explained, it comes from a sequence present in all three samples with MD5
checksum ``5122dde24762f8e3d6a54e3f79077254``, and this exact sequence is in
the database with entries for both *Phytophthora castaneae* (which was in the
DNA control mixture) and also *Phytophthora agathidicida* (e.g. accession
KP295308).

You can confirm this by looking at the intermediate TSV files, e.g. using
grep to show all lines with this species name:

.. code:: console

    $ grep "Phytophthora agathidicida" intermediate/DNA10MIX_*.onebp.tsv
    intermediate/DNA10MIX_bycopynumber.onebp.tsv:29de890989becddc5e0b10ecbbc11b1a_245  1642459;1642465  Phytophthora agathidicida;Phytophthora castaneae
    intermediate/DNA10MIX_diluted25x.onebp.tsv:29de890989becddc5e0b10ecbbc11b1a_656    1642459;1642465  Phytophthora agathidicida;Phytophthora castaneae
    intermediate/DNA10MIX_undiluted.onebp.tsv:29de890989becddc5e0b10ecbbc11b1a_625     1642459;1642465  Phytophthora agathidicida;Phytophthora castaneae

The same applies to *Phytophthora capsici* and *Phytophthora glovera*,
although in this case both were in the mixture.

Overall, given the uniqueness limitations of the ITS1 marker, the tool has
done a faultless job on these three positive control samples from the ten
species mix.

Running thapbi_pict assess
--------------------------

Comparing a few samples like this by hand is one thing, but doing it at scale
requires automation. For assessing changes to the classifier method and
database, we mainly run ``thapbi_pict assess`` against a set of single isolate
positive controls. This requires a computer readable files listing the
expected species in a particular format.

.. code:: console

    $ thapbi_pict assess -h
    ...

The inputs to this command are pairs of plain text tab separated variable
(TSV) files named ``<sample_name>.known.tsv`` (the expected results) and
``<sample_name>.<method>.tsv`` which is the intermediate TSV file from
running ``thapbi_pict classify`` on ``<sample_name>.fasta``, which in turn
came from running ``thapbi_pict prepare-reads`` on the the pair
``<sample_name>_R1.fastq.gz`` and ``<sample_name>_R2.fastq.gz``.

The "known" file uses the same column based layout as the intermediate TSV
files, but while you can provide the expected species for each unique sequence
in the sample, this can be simiplified to a single wildcard ``*`` line
followed by all the NCBI taxids and species names using semi-colon separators.

Looking at the 15 species mixture, we want to assess the classification in the
file ``intermediate/DNA15MIX.onebp.tsv`` so we will need a file named
``DNA15MIX.known.tsv``. This can be in any folder, but the convention we use
is another folder ``expected/`` for all the ``*.known.tsv`` files.
See :ref:`sample data setup <sample_data>` for where to get this file.

The simplest way to run the assess command is to tell it two input filenames,
and it will default to printing its tabular output to screen - shown here
abridged after piping through the ``cut`` command to pull out just the first
five columns:

.. code:: console

    $ thapbi_pict assess -i expected/DNA15MIX.known.tsv intermediate/DNA15MIX.onebp.tsv | cut -f 1-5
    Assessed onebp vs known in 1 files (1099 species)
    #Species                     TP  FP  FN  TN
    OVERALL                      8   1   7   1083
    Phytophthora aleatoria       0   1   0   0
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
    OTHER 1083 SPECIES IN DB     0   0   0   1083

More usually, you would output to a named file, and look at that:

.. code:: console

    $ thapbi_pict assess -i expected/DNA15MIX.known.tsv intermediate/DNA15MIX.onebp.tsv -o DNA15MIX.assess.tsv
    Assessed onebp vs known in 1 files (1099 species)
    $ cut -f 1-5,11 DNA15MIX.assess.tsv
    <SEE TABLE BELOW>

You should be able to open this ``DNA15MIX.assess.tsv`` file in R, Excel, etc:

=========================== == == == ==== ===========
#Species                    TP FP FN TN   Ad-hoc-loss
=========================== == == == ==== ===========
OVERALL                     8  1  7  1083 0.500
Phytophthora aleatoria      0  1  0  0    1.000
Phytophthora austrocedri    1  0  0  0    0.000
Phytophthora boehmeriae     0  0  1  0    1.000
Phytophthora cactorum       1  0  0  0    0.000
Phytophthora cambivora      0  0  1  0    1.000
Phytophthora chlamydospora  0  0  1  0    1.000
Phytophthora cinnamomi      0  0  1  0    1.000
Phytophthora gonapodyides   1  0  0  0    0.000
Phytophthora ilicis         1  0  0  0    0.000
Phytophthora kernoviae      1  0  0  0    0.000
Phytophthora lateralis      0  0  1  0    1.000
Phytophthora obscura        1  0  0  0    0.000
Phytophthora plurivora      0  0  1  0    1.000
Phytophthora pseudosyringae 1  0  0  0    0.000
Phytophthora ramorum        1  0  0  0    0.000
Phytophthora syringae       0  0  1  0    1.000
OTHER 1083 SPECIES IN DB    0  0  0  1083 0.000
=========================== == == == ==== ===========

The ``OVERALL`` line tells us that there were 8 true positives, 1 false
positives, 7 false negatives, and 1083 true negatives. The final number needs a
little explanation. First, 8+1+7+1083 = 1099, which is the number of species in
the database. With only one sample being considered, 1083 is the number of other
species in the database which the tool correctly says are not present.

The additional columns (not shown) include traditional metrics like
sensitivity, specificity, precision, F1, Hamming loss, plus our own metric
provisionally called *Ad hoc loss* which is a modification of the Hamming loss
without using the true negative count (which we expect to always be very large
as the database will contain many species, while a community might contain
only ten).

Following this we get one line per species, considering this species in
isolation (making this a traditional and simpler to interpret classification
problem). Here there is only one sample, so this time TP+FP+FN+TN=1.

Next, let's run the assess command on all four positive control samples, just
by giving the input directory names (it will work out the common filenames):

.. code:: console

    $ thapbi_pict assess -i expected/ intermediate/ -o thabpi-pict.assess.tsv
    Assessed onebp vs known in 4 files (1099 species)
    $ cut -f 1-5,11 thabpi-pict.assess.tsv
    <SEE TABLE BELOW>

New table ``thabpi-pict.assess.tsv`` is similar:

=========================== == == == ==== ===========
#Species                    TP FP FN TN   Ad-hoc-loss
=========================== == == == ==== ===========
OVERALL                     32 7  13 4344 0.385
Phytophthora agathidicida   0  3  0  1    1.000
Phytophthora aleatoria      0  1  0  3    1.000
Phytophthora austrocedri    1  0  0  3    0.000
Phytophthora boehmeriae     0  0  4  0    1.000
Phytophthora cactorum       1  0  3  0    0.750
Phytophthora cambivora      0  0  1  3    1.000
Phytophthora capsici        3  0  0  1    0.000
Phytophthora castaneae      3  0  0  1    0.000
Phytophthora chlamydospora  0  0  1  3    1.000
Phytophthora cinnamomi      0  0  1  3    1.000
Phytophthora fallax         3  0  0  1    0.000
Phytophthora foliorum       3  0  0  1    0.000
Phytophthora glovera        0  3  0  1    1.000
Phytophthora gonapodyides   1  0  0  3    0.000
Phytophthora ilicis         1  0  0  3    0.000
Phytophthora kernoviae      1  0  0  3    0.000
Phytophthora lateralis      0  0  1  3    1.000
Phytophthora obscura        4  0  0  0    0.000
Phytophthora plurivora      3  0  1  0    0.250
Phytophthora pseudosyringae 1  0  0  3    0.000
Phytophthora ramorum        1  0  0  3    0.000
Phytophthora rubi           3  0  0  1    0.000
Phytophthora siskiyouensis  3  0  0  1    0.000
Phytophthora syringae       0  0  1  3    1.000
OTHER 1075 SPECIES IN DB    0  0  0  4300 0.000
=========================== == == == ==== ===========

This time the ``OVERALL`` line says we had 32 TP, 7 FP, 13 FN and 4344 TN.
That total 32+7+13+4344 = 4396 = 4 * 1099, the number of samples times the
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

    $ thapbi_pict pipeline -i raw_data/ expected/ \
      -s intermediate/ -o summary/ \
      -n raw_data/NEGATIVE*.fastq.gz -r with-metadata \
      -t metadata.tsv -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 -x 16 -f 20
    ...
    $ ls -1 summary/with-metadata.*
    summary/with-metadata.all_reads.fasta
    summary/with-metadata.assess.confusion.onebp.tsv
    summary/with-metadata.assess.onebp.tsv
    summary/with-metadata.assess.tally.onebp.tsv
    summary/with-metadata.edit-graph.onebp.xgmml
    summary/with-metadata.reads.onebp.tsv
    summary/with-metadata.reads.onebp.xlsx
    summary/with-metadata.samples.onebp.tsv
    summary/with-metadata.samples.onebp.txt
    summary/with-metadata.samples.onebp.xlsx

Output file ``summary/with-metadata.assess.onebp.tsv`` will match the output
above.
