.. _custom_database_pipeline

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
      -r recycled-water-custom-onebp -t metadata.tsv -c 1,2,3,4,5,6,7
    ...
    onebp classifier assigned species/genus to 3577559 of 9956078 sequences from 384 files
    ...
    Wrote summary/recycled-water-custom-onebp.samples.*
    Wrote summary/recycled-water-custom-onebp.reads.*
    Loaded 3054 unique sequences from 384 FASTA files.
    Matched 81 unique sequences in database
    ...
    $ ls -1 intermediate/SRR*.fasta | wc -l
    384
    $ ls -1 intermediate/SRR*.onebp.tsv | wc -l
    384
    $ ls -1 summary/summary/thapbi-pict.*
    recycled-water-custom-onebp.reads.tsv
    recycled-water-custom-onebp.reads.xlsx
    recycled-water-custom-onebp.samples.tsv
    recycled-water-custom-onebp.samples.txt
    recycled-water-custom-onebp.edit-graph.xgmml

Note here I have made the classifier method explicit with ``-m`` (or
``--method``), using the default of ``onebp``. With this narrower set
of *Phytophthora* sequences and comparatively well sampled database,
that was a good default. Here with our relatively sparse database, is
is perhaps overly strict - only 36% of the sequences were matched.

The ``blast`` based classifier makes more matches:

.. code:: console

    $ thapbi_pict pipeline -i raw_data/ -s intermediate/ -o summary/ \
      --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA \
      --right AGCGTTCTTCATCGATGTGC \
      -d Redekar_et_al_2019_sup_table_3.sqlite -m blastp \
      -r recycled-water-custom-blast -t metadata.tsv -c 1,2,3,4,5,6,7
    ...
    blast classifier assigned species/genus to 4281041 of 9956078 sequences from 384 files
    ...
    Wrote summary/recycled-water-custom-blast.samples.*
    Wrote summary/recycled-water-custom-blast.reads.*
    Loaded 3054 unique sequences from 384 FASTA files.
    Matched 81 unique sequences in database
    ...
    $ ls -1 intermediate/SRR*.blast.tsv | wc -l
    384

Better, but still only about 43% of the sequences classifiers.
