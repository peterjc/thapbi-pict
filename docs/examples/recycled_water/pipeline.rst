.. _custom_database_pipeline:

Pipeline with custom database
=============================

Running thapbi-pict pipeline
----------------------------

Compared to the original worked example, we must specify the primers and
our custom database built with matching primer trimmed entries:

.. code:: console

    $ mkdir intermediate/ summary/
    $ thapbi_pict pipeline -i raw_data/ -s intermediate/ -o summary/ \
      --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA \
      --right AGCGTTCTTCATCGATGTGC \
      -d Redekar_et_al_2019_sup_table_3.sqlite -m onebp \
      -r recycled-water-custom -t metadata.tsv -x 7 -c 1,2,3,4,5,6
    ...
    onebp classifier assigned species/genus to 3577559 of 9956078 sequences from 384 files
    ...
    Wrote summary/recycled-water-custom.samples.onebp.*
    Wrote summary/recycled-water-custom.reads.onebp.*
    Loaded 3054 unique sequences from 384 FASTA files.
    Matched 81 unique sequences in database
    ...
    $ ls -1 intermediate/SRR*.fasta | wc -l
    384
    $ ls -1 intermediate/SRR*.onebp.tsv | wc -l
    384
    $ ls -1 summary/summary/thapbi-pict.*
    recycled-water-custom.reads.onebp.tsv
    recycled-water-custom.reads.onebp.xlsx
    recycled-water-custom.samples.onebp.tsv
    recycled-water-custom.samples.onebp.txt
    recycled-water-custom.edit-graph.xgmml

Note the classifier method was set explicitly with ``-m`` (or ``--method``),
using the default of ``onebp``. With the narrower set of *Phytophthora*
sequences and comparatively well sampled database, that was a good default.
Recall running with the *Phytophthora* defaults gave a taxonomic assignment
for 1880048 of 2605870 of 2605870 reads - which was 72% of 2.6 million reads.

Here with our relatively sparse database, the ``onebp`` method is perhaps
overly strict - only 36% of the reads were matched (3577559 of 9956078).
However, with the different primer settings, we are examining almost ten million
reads (nearly four times as many), so we're doing about twice as well in terms
of number of raw reads with a classification.

Naturally the more lenient or fuzzy ``blast`` based classifier makes even
more matches:

.. code:: console

    $ thapbi_pict pipeline -i raw_data/ -s intermediate/ -o summary/ \
      --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA \
      --right AGCGTTCTTCATCGATGTGC \
      -d Redekar_et_al_2019_sup_table_3.sqlite -m blast \
      -r recycled-water-custom -t metadata.tsv -x 7 -c 1,2,3,4,5,6
    ...
    blast classifier assigned species/genus to 4281041 of 9956078 sequences from 384 files
    ...
    Wrote summary/recycled-water-custom.samples.blast.*
    Wrote summary/recycled-water-custom.reads.blast.*
    Loaded 3054 unique sequences from 384 FASTA files.
    Matched 81 unique sequences in database
    ...
    $ ls -1 intermediate/SRR*.blast.tsv | wc -l
    384

Better, in that we are up to 43% of the reads with a taxonomic assignment
(4281041 of 9956078 reads). But how many of these are false positives? Sadly,
we don't have any controls for this dataset in order to objectively assess the
classifier performance of the various algorithm and database combinations.

However we can say that this database and indeed the published *Oomycetes*
ITS1 sequences in general is relatively sparse outside *Phytophthora* (and
even there, we as a community have room for improvement).

Results
-------

We will focus on the same four low diversity samples for a brief comparison
of the defaults, custom DB with ``onebp``, and custom DB with ``blast``.
This information is extracted from the sample reports (new files
``summary/recycled-water-custom.samples.onebp.txt`` and
``summary/recycled-water-custom.samples.blast.txt``).

``SRR6303586`` aka ``OSU483``:

- With defaults:

  - *Phytophthora chlamydospora*

- With custom DB and ``onebp`` or ``blast``:

  - *Phytophthora sp.* CAL-2011b (uncertain/ambiguous)
  - *Phytophthora chlamydospora*

``SRR6303588`` aka ``OSU482``:

- With defaults:

  - Unknown
  - *Phytophthora chlamydospora*

- With custom DB and ``onebp`` or ``blast``:

  - *Phytophthora sp.* CAL-2011b (uncertain/ambiguous)
  - *Phytophthora chlamydospora*

``SRR6303596`` aka ``OSU121``:

- With defaults:

  - Unknown

- With custom DB and ``onebp`` or ``blast``:

  - *Phytopythium litorale*
  - *Pythium aff. diclinum* (uncertain/ambiguous)
  - *Pythium aff. dictyosporum* (uncertain/ambiguous)
  - *Pythium aff. dissotocum* (uncertain/ambiguous)
  - *Pythium cf. dictyosporum* (uncertain/ambiguous)
  - *Pythium coloratum* (uncertain/ambiguous)
  - *Pythium diclinum* (uncertain/ambiguous)
  - *Pythium dissotocum* (uncertain/ambiguous)
  - *Pythium lutarium*
  - *Pythium sp.* CAL-2011f (uncertain/ambiguous)
  - *Pythium sp.* group F (uncertain/ambiguous)

``SRR6303948`` aka ``OSU536.s203``:

- With defaults:

  - *Phytophthora ramorum*

- With custom DB and ``onebp`` or ``blast``:

  - Unknown
  - *Phytophthora ramorum*

So, not too dramatic - and on this subset using ``onebp`` versus ``blast``
seems not to matter.

Interestingly the two databases differ on exactly which *Phytophthora* are
present. The main change is with these settings and the new database
``SRR6303596`` aka ``OSU121`` has multiple *Pythium* results (why this
example was selected) plus *Phytopythium litorale* (originally known as
*Pythium litoralis*), and ``SRR6303948`` has some unknown *Oomycete(s)* (as
discussed earlier at the end of the :ref:`primers <custom_database_primers>`
section).
