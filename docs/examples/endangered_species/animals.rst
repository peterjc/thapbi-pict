Universal animal DNA barcodes and mini-barcodes
===============================================

For 16S, COI and cyt-b the paper used two targets, a long barcode and a shorter
mini-barcode. The same names have been used in the ``run.sh`` script provided,
the output of which is referred to below.

16S - long marker
-----------------

The 16S primer set output is disappointing at the default abundance threshold,
with only a single unique sequence observed - I suspect the long product size
is part of the issue, it must be at the upper limit for overlapping MiSeq read
pairs?

.. code:: console

    $ grep -v "^#" summary/16S.tally.tsv | cut -f 1,179
    16S/1f2b15d58f9f40b862486676809d4744_20189  CACCTCCAGCATTCCCAGTATTGGAGGCATTGCCTGCCCAGTGACAACTGTTTAACGGCCGCGGTATCCTGACCGTGCAAAGGTAGCATAATCATTTGTTCTCTAAATAAGGACTTGTATGAATGGCCGCACGAGGGTTTTACTGTCTCTTACTTCCAATCAGTGAAATTGACCTTCCCGTGAAGAGGCGGGAATGCACAAATAAGACGAGAAGACCCTATGGAGCTTTAACTAACCAACCCAAAGAGAATAGATTTAACCATTAAGGAATAACAACAATCTCCATGAGTTGGTAGTTTCGGTTGGGGTGACCTCGGAGAATAAAAAATCCTCCGAGCGATTTTAAAGACTAGACCCACAAGTCAAATCACTCTATCGCTCATTGATCCAAAAACTTGATCAACGGAACAAGTTACCCTAGGGATAACAGCGCAATCCTATTCAAGAGTCCATATCGACAATAGGGTTTACGACCTCGATGTTGGATCAGGACATCCTGATGGTGCAACCGCTATCAAAGGTTCGTTTGTTCAACGATTAAAGTCCT

This perfectly matches *Bos taurus* and was found in most but not all of the
samples expected - perhaps the default abundance threshold is too high?

Mini-16S - short marker
-----------------------

The output from the Mini-16S marker is far more diverse, with 84 unique
sequences:

.. code:: console

    $ grep -c -v "^#" summary/Mini-16S.tally.tsv
    84

The most common is again a perfect match to *Bos taurus*, which this time has
no false negatives (but two false positives?).

We have all the expected *Sus scrofa* matches, and some of *Gallus gallus* and
*Anguilla anguilla* expected in six samples. *Crocodylus niloticus* is also
found but at far lower levels than expected.

We do see *Homo sapiens*, but happily only in the traditional medicine samples
(multiple replicates within ``S3`` and ``S8``). Within those samples, the
laboratory 16 replicates ``S3_Lab_16`` and ``S8_Lab_16`` also had *Rattus
tanezumi* and *Rattus norvegicus* too, respectively.

Overall, again perhaps the default abundance threshold is too high?

COI - long marker
-----------------

Assuming I understood the paper correctly, this used a pool of four left
primers and four right primers. That is not easily handled with THAPBI PICT at
the time of writing.

Mini-COI - short marker
-----------------------

The output from the Mini-COI marker is quite diverse, with 22 unique sequences:

.. code:: console

    $ grep -c -v "^#" summary/Mini-COI.tally.tsv
    22

The species matches are all reasonable, it detects all the *Pieris brassicae*,
most of the *Bos taurus*, *Pleuronectes platessa*, *Sus scrofa*, many of the
*Huso dauricus* and *Gallus gallus*.

We have unexpected *Acipenser schrenckii*, which was also found in the paper
and explained due to sample preparation.

There are also plenty of unclassified sequences from the traditional medicine
samples, based on an NCBI BLAST search many are likely from undescribed fungi.

cyt-b - long marker
-------------------

This gave no sequences at the default abundance threshold, nor at 50. Dropping
to 10 we get a modest number of hits - the only perfect match was unfortunately
to plants in the Asteraceae family.

Mini-cyt-b - short marker
-------------------------

The output from the Mini-COI marker had only 17 unique sequences:

.. code:: console

    $ grep -c -v "^#" summary/Mini-cyt-b.tally.tsv
    17

This found all the expected *Sus scrofa* and *Meleagris gallopavo*, and most
*Bos taurus*, *Crocodylus niloticus*, *Huso dauricus* and some of the
*Anguilla anguilla*.

As above, we have explained false matches for *Acipenser schrenckii*, and
again *Homo sapiens* in the traditional medicine but also in ``EM_8``.
