.. _custom_database_pipeline:

Pipeline with custom database
=============================

Running thapbi-pict pipeline
----------------------------

Compared to the original worked example, we must specify our custom database
(which contains the primer information, and matching primer trimmed entries):

.. code:: console

    $ mkdir -p intermediate_long/ summary/
    $ thapbi_pict pipeline -i raw_data/ -s intermediate_long/ \
      -o summary/recycled-water-custom \
      -d Redekar_et_al_2019_sup_table_3.sqlite -m onebp \
      -t metadata.tsv -x 7 -c 1,2,3,4,5,6
    ...
    onebp classifier assigned species/genus to 3567127 of 10933210 sequences from 1 files
    Wrote summary/recycled-water-custom.ITS1-long.samples.onebp.*
    Wrote summary/recycled-water-custom.ITS1-long.reads.onebp.*
    ...
    $ ls -1 intermediate_long/ITS1-long/SRR*.fasta | wc -l
    384
    $ ls -1 summary/recycled-water-custom.*.onebp.*
    summary/recycled-water-custom.ITS1-long.all_reads.onebp.tsv
    summary/recycled-water-custom.ITS1-long.reads.onebp.tsv
    summary/recycled-water-custom.ITS1-long.reads.onebp.xlsx
    summary/recycled-water-custom.ITS1-long.samples.onebp.tsv
    summary/recycled-water-custom.ITS1-long.samples.onebp.txt
    summary/recycled-water-custom.ITS1-long.samples.onebp.xlsx

Note the classifier method was set explicitly with ``-m`` (or ``--method``),
using the default of ``onebp``. With the narrower set of *Phytophthora*
sequences and comparatively well sampled database, that was a good default.
Recall running with the *Phytophthora* defaults gave a taxonomic assignment
for 2122727 of 2598668 reads - which was 82% of 2.6 million reads.

Here with our relatively sparse database, the ``onebp`` method is perhaps
overly strict - only 33% of the reads were matched (3567127 of 10933210).
However, with the different primer settings, we are examining over ten
million reads (nearly four times as many), so we're doing about twice as well
in terms of number of raw reads with a classification.

Naturally the more lenient or fuzzy ``blast`` based classifier makes even
more matches:

.. code:: console

    $ thapbi_pict pipeline -i raw_data/ -s intermediate_long/ \
      -o summary/recycled-water-custom \
      -d Redekar_et_al_2019_sup_table_3.sqlite -m blast \
      -t metadata.tsv -x 7 -c 1,2,3,4,5,6
    ...
    blast classifier assigned species/genus to 4268503 of 10933210 sequences from 1 files
    Wrote summary/recycled-water-custom.ITS1-long.samples.blast.*
    Wrote summary/recycled-water-custom.ITS1-long.reads.blast.*
    ...
    $ ls -1 summary/recycled-water-custom.*.blast.*
    summary/recycled-water-custom.ITS1-long.all_reads.blast.tsv
    summary/recycled-water-custom.ITS1-long.reads.blast.tsv
    summary/recycled-water-custom.ITS1-long.reads.blast.xlsx
    summary/recycled-water-custom.ITS1-long.samples.blast.tsv
    summary/recycled-water-custom.ITS1-long.samples.blast.txt
    summary/recycled-water-custom.ITS1-long.samples.blast.xlsx

Better, in that we are up to 39% of the reads with a taxonomic assignment
(4268503 of 10933210 reads). But how many of these are false positives? Sadly,
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
