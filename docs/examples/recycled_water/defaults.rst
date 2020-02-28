.. _custom_database_defaults:

Pipeline with defaults
======================

Running thapbi-pict pipeline
----------------------------

First, we will run the THAPBI PICT pipeline command with largely default
settings, other than including the metadata about the water samples. Note
that this dataset has no blanks or negative controls, so must trust to the
default minimum abundance threshold.

The key values which we will be changing later are the primers and database.

Assuming you have the FASTQ files in ``raw_data/``, run the pipeline command
as follows, and you should get six output report files:

.. code:: console

    $ mkdir intermediate_defaults/ summary_defaults/
    $ thapbi_pict pipeline \
      -i raw_data/ -o summary_defaults/ -s intermediate_defaults/ \
      -r recycled-water-defaults -t metadata.tsv -x 7 -c 1,2,3,4,5,6
    ...
    onebp classifier assigned species/genus to 1880048 of 2605870 sequences from 384 files
    ...
    $ ls -1 intermediate_defaults/SRR*.fasta | wc -l
    384
    $ ls -1 intermediate_defaults/SRR*.onebp.tsv | wc -l
    384
    $ ls -1 summary/summary_defaults/thapbi-pict.*
    recycled-water-defaults.reads.onebp.tsv
    recycled-water-defaults.reads.onebp.xlsx
    recycled-water-defaults.samples.onebp.tsv
    recycled-water-defaults.samples.onebp.xlsx
    recycled-water-defaults.samples.onebp.txt
    recycled-water-defaults.edit-graph.xgmml

Here we used ``-r`` (or ``--report``) to specify a different stem for the
report filenames. The :ref:`sample metadata options <metadata>` were described
earlier -- this is perhaps an idealised example in that ``metadata.tsv`` was
created so that we add the first six columns the table (sorted in that order),
where ``-x 7`` means index to the accession (filename prefix) in column seven.

Notice the output reported a taxonomic assignment for 1880048 of 2605870
reads - that's 72%.

Results
-------

We will compare and contrast the following four samples with the second run
using different primers and a custom database. These were deliberatly picked
from the less diverse samples for clarity.

For now, here is a formatted excerpt from the sample report in file
``recycled-water-defaults.samples.onebp.txt``:

    :Accession: SRR6303586
    :Sample: OSU483
    :Source: Reservoir
    :Site: K
    :Process: Leaf baiting
    :Period: 18
    :Year-Month: 2016-01

    Sequencing sample: SRR6303586

    - *Phytophthora chlamydospora*

    :Accession: SRR6303588
    :Sample: OSU482
    :Source: Reservoir
    :Site: J
    :Process: Leaf baiting
    :Period: 18
    :Year-Month: 2016-01

    Sequencing sample: SRR6303588

    - Unknown
    - *Phytophthora chlamydospora*

    :Accession: SRR6303596
    :Sample: OSU121
    :Source: Runoff
    :Site: H
    :Process: Leaf baiting
    :Period: 2
    :Year-Month: 2015-05

    Sequencing sample: SRR6303596

    - Unknown

    :Accession: SRR6303948
    :Sample: OSU536.s203
    :Source: Runoff
    :Site: H
    :Process: Filtration
    :Period: 22
    :Year-Month: 2016-03

    Sequencing sample: SRR6303948

    - *Phytophthora ramorum*

Three of these four have *Phytophthora*, and two have some unknown(s).
However, this is discarding all the reads which do not match the default
*Phytophthora* centric primers.
