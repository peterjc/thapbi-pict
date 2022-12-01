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
    onebp classifier assigned species/genus to 529 of 3053 unique sequences from 1 files
    Wrote summary/recycled-water-custom.ITS1-long.samples.onebp.*
    Wrote summary/recycled-water-custom.ITS1-long.reads.onebp.*
    ...
    $ ls -1 summary/recycled-water-custom.*.onebp.*
    summary/recycled-water-custom.ITS1-long.all_reads.onebp.tsv
    summary/recycled-water-custom.ITS1-long.reads.onebp.tsv
    summary/recycled-water-custom.ITS1-long.reads.onebp.xlsx
    summary/recycled-water-custom.ITS1-long.samples.onebp.tsv
    summary/recycled-water-custom.ITS1-long.samples.onebp.xlsx

Note the classifier method was set explicitly with ``-m`` (or ``--method``),
using the default of ``onebp``. With the narrower set of *Phytophthora*
sequences and comparatively well sampled database, that was a good default.
Recall running with the *Phytophthora* defaults gave a taxonomic assignment
for 2122757 of 2598566 reads - which was 82% of 2.6 million reads.

Here with our relatively sparse database, the ``onebp`` method is perhaps
overly strict - only 17% of the unique sequences matched (529 of 3053 ASVs),
although it is more like a third if we count the number of reads matched.
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
    blast classifier assigned species/genus to 1036 of 3053 unique sequences from 1 files
    Wrote summary/recycled-water-custom.ITS1-long.samples.blast.*
    Wrote summary/recycled-water-custom.ITS1-long.reads.blast.*
    ...
    $ ls -1 summary/recycled-water-custom.*.blast.*
    summary/recycled-water-custom.ITS1-long.all_reads.blast.tsv
    summary/recycled-water-custom.ITS1-long.reads.blast.tsv
    summary/recycled-water-custom.ITS1-long.reads.blast.xlsx
    summary/recycled-water-custom.ITS1-long.samples.blast.tsv
    summary/recycled-water-custom.ITS1-long.samples.blast.xlsx

Better, in that we are up to 34% of the unique sequences with a taxonomic
assignment (1036 of 3053 ASVs). But how many of these are false positives?
Sadly, we don't have any controls for this dataset in order to objectively
assess the classifier performance of the various algorithm and database
combinations.

However we can say that this database and indeed the published *Oomycetes*
ITS1 sequences in general is relatively sparse outside *Phytophthora* (and
even there, we as a community have room for improvement).

Results
-------

We will focus on the same four low diversity samples for a brief comparison
of the defaults, custom DB with ``onebp``, and custom DB with ``blast``.

Previously with the default DB and default ``onebp`` classifier:

.. code:: console

    $ cut -f 6,7,8 summary/recycled-water-defaults.ITS1.samples.onebp.tsv \
      | grep -E "(SRR6303586|SRR6303586|SRR6303588|SRR6303596|SRR6303948)"
    OSU482       SRR6303588  Phytophthora chlamydospora, Phytophthora x stagnum(*), Unknown
    OSU483       SRR6303586  Phytophthora chlamydospora, Phytophthora x stagnum(*)
    OSU536.s203  SRR6303948  Phytophthora ramorum
    OSU121       SRR6303596  Phytopythium (unknown species)

With the custom DB:

.. code:: console

    $ cut -f 6,7,8 summary/recycled-water-custom.ITS1-long.samples.onebp.tsv \
      | grep -E "(SRR6303586|SRR6303586|SRR6303588|SRR6303596|SRR6303948)"
    OSU482       SRR6303588  Phytophthora chlamydospora, Phytophthora sp. CAL-2011b(*)
    OSU483       SRR6303586  Phytophthora chlamydospora, Phytophthora sp. CAL-2011b(*)
    OSU536.s203  SRR6303948  Phytophthora ramorum, Unknown
    OSU121       SRR6303596  Phytopythium litorale, Pythium aff. diclinum(*), Pythium aff. dictyosporum(*), Pythium aff. dissotocum(*), Pythium cf. dictyosporum(*), Pythium coloratum(*), Pythium diclinum(*), Pythium dissotocum(*), Pythium lutarium, Pythium sp. CAL-2011f(*), Pythium sp. group F(*)

We get the same using the top BLAST hit:

.. code:: console

    $ cut -f 6,7,8 summary/recycled-water-custom.ITS1-long.samples.blast.tsv \
      | grep -E "(SRR6303586|SRR6303586|SRR6303588|SRR6303596|SRR6303948)"
    OSU482       SRR6303588  Phytophthora chlamydospora, Phytophthora sp. CAL-2011b(*)
    OSU483       SRR6303586  Phytophthora chlamydospora, Phytophthora sp. CAL-2011b(*)
    OSU536.s203  SRR6303948  Phytophthora ramorum, Unknown
    OSU121       SRR6303596  Phytopythium litorale, Pythium aff. diclinum(*), Pythium aff. dictyosporum(*), Pythium aff. dissotocum(*), Pythium cf. dictyosporum(*), Pythium coloratum(*), Pythium diclinum(*), Pythium dissotocum(*), Pythium lutarium, Pythium sp. CAL-2011f(*), Pythium sp. group F(*)

On this subset using ``onebp`` versus ``blast`` seems not to matter.
The sample report does not go down to the sequences in each sample,
for that you can use the reads report, or look at the intermediate
FASTA files as discussed in the previous :ref:`primers
<custom_database_primers>` section.

The first two example differ due to the DB curation about exactly
which *Phytophthora* is present. Sample ``OSU121`` aka ``SRR6303596``
went from one *Phytopythium litorale* sequence to being joined
by a much more numerous *Pythium coloratum/dissotocum* sequence
(plus some lower abundance variants of it). Likewise,
``OSU536.s203`` aka ``SRR6303948`` had one sequence for
*Phytophthora ramorum*, but now has multiple unknown sequences.
