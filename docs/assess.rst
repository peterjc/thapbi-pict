.. _classifier_assessment:

Classifier Assessment
=====================
In assessing classification performance, it is the combination of both
classification method (algorithm) and marker database which which matters.
Settings like the abundance threshold are also important, and the tool default
settings partly reflect one of the original project goals being to avoid false
positives.

To objectively assess a metabarcoding classifier we require sequenced samples
of known composition, which generally means single isolates (where a single
marker sequence is typically expected), or mock communities (the bulk of our
:ref:`worked examples <worked_examples>`). Carefully controlled environmental
samples are possible too. We use Muri *et al.* (2020) as a :ref:`worked
example identifying fish species <drained_ponds>` where the lake was drained
to collected and identify all the individual fish, but this is problematic as
the lakes were large enough that DNA from each fish could not be expected at
all the sampling points, giving an inflated false negative count.

Our tool includes an presence/absence based assessment framework based on
supplying expected species lists for control samples, from which the true
positive (TP), false positive (FP), true negative (TN), and false negative
(FN) counts can be computed for each species. These are the basis of standard
metrics like sensitivity (recall), specificity, precision, F-score (F-measure,
or F1), and Hamming Loss. It is simple but not overly helpful to apply metrics
like this to each species, but the overall performance is more informative.

However, some scores like the Hamming Loss are fragile with regards to the TN
count when comparing databases. The Hamming Loss given by the total number of
mis-predicted class entries divided by the number of class-level predictions,
thus (FP + FN) / (TP + FP + FN + TN).
Consider a mock community of ten species, where the classifier made 11
predictions which break down as 9 TP and 2 FP, meaning 10 - 9 = 1 FN.
Suppose the database had a hundred species (including all ten in the mock
community), that leaves 100 - 9 - 1 - 2 = 88 TN, and a Hamming Loss of 3/100
= 0.03. Now suppose the database was extended with additional references not
present in this mock community, perhaps expanding from European *Phytophthora*
species to include distinct entries for tropical species, or a sister group
like *Peronospora*. The denominator would increase, reducing the Hamming Loss,
but intuitively the classifier performance on this mock community has not
changed. To address this, the classifier assessment also includes a modified
*ad-hoc* loss metric calculated as the total number of mis-predicted class
entries divided by the number of class-level predictions ignoring TN, or
(FP + FN) / (TP + FP + FN) which in this example would give 3/12 = 0.25
regardless of the number of species in the database. This is an intuitive
measure weighting FP and FN equally (smaller is better, zero is perfect),
a potential complement to the F-score.

Note that the assessment framework only considers species level predictions,
ignoring genus only predictions, and thus will not distinguish between the
default ``onebp`` classifier and variants like ``1s3g``.
