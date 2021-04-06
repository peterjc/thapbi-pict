Classifying sequences
=====================

Running thapbi-pict-classify
----------------------------

.. tip:

   If you don't have the FASTQ files, just the FASTA files, start from here.

The second stage of the pipeline can be run separately as the ``thapbi_pict
classify`` command:

.. code:: console

    $ thapbi_pict classify -h
    ...

There are a number of options here, but for the purpose of this worked example
we will stick with the defaults and tell it to look for FASTA files in the
``intermediate/`` directory.

.. code:: console

    $ thapbi_pict classify -i intermediate/
    ...

Here we have not set the output folder with ``-o`` or ``--output``, which
means the tool will default to writing the TSV output files next to each
input FASTA file. There should now be 122 TSV files, one for each FASTA:

.. code:: console

    $ ls -1 intermediate/*.tsv | wc -l
    122

Intermediate TSV files
----------------------

For each FASTA file named ``<sample_name>.fasta`` a plain text tab separated
variable (TSV) file is generated named ``<sample_name>.<method>.tsv`` where
the default method is ``onebp`` (this looks for perfect matches or up to one
base pair different). The first line is a header comment line (starting with
``#``) labelling the columns, which are:

* Unique sequence name in ``<checksum>_<abundance>`` format.
* NCBI taxid of any predictions (semi-colon separated, same order as species)
* Genus-species of any predictions (semi-colon separated, alphabetical)
* Text note field (arbitrary debug text from the tool)

These files are not really intended for human use, but are readable:

.. code:: console

    $ cat intermediate/Site_1_sample_1.onebp.tsv
    ...

Viewing it like this is not ideal, although there are command line tools which
help. You could open the file in R, Excel, etc. Slightly abridged and
reformatted, we have:

========================================= ============== ==================================================== ====
#sequence-name                            taxid          genus-species                                        note
========================================= ============== ==================================================== ====
``2e4f0ed53888ed39a2aee6d6d8e02206_2271`` 221518         *Phytophthora pseudosyringae*                        ...
``c1a720b2005f101a9858107545726123_716``  78237          *Phytophthora gonapodyides*                          ...
``96e0e2f0475bd1617a4b05e778bb04c9_331``  78237          *Phytophthora gonapodyides*                          ...
``fb30156d7f66c8abf91f9da230f4d19e_208``  164328         *Phytophthora ramorum*                               ...
``dcd6316eb77be50ee344fbeca6e005c7_193``  164328         *Phytophthora ramorum*                               ...
``972db44c016a166de86a2bacab3f4226_182``  53983; 2056922 *Phytophthora cambivora*; *Phytophthora x cambivora* ...
``d9bc3879fdab3b4184c04bfbb5cf6afb_165``  631361         *Phytophthora austrocedri*                           ...
``ed15fefb7a3655147115fc28a8d6d671_113``  78237          *Phytophthora gonapodyides*                          ...
========================================= ============== ==================================================== ====

This says most of the unique sequences here have been assigned a single unique
*Phytophthora* species, except for ``972db44c016a166de86a2bacab3f4226`` (found
in 182 reads for this sample) which has matched *Phytophthora cambivora* (NCBI
taxid 53983) and close relative *Phytophthora x cambivora* (NCBI taxid
2056922).

If you are familiar with the command line search tool ``grep`` and the regular
expression syntax, you should find the format of these intermediate TSV files
lends itself to some simple searches. For example, you could see which samples
had matches to *Phytophthora rubi* using ``grep`` as follows:

.. code:: console

    $ grep "Phytophthora rubi" intermediate/*.tsv
    intermediate/DNA10MIX_bycopynumber.onebp.tsv:d8613e80b8803b13f7ea5d097f8fe46f_279  129364  Phytophthora rubi  Unique taxonomy match
    intermediate/DNA10MIX_diluted25x.onebp.tsv:d8613e80b8803b13f7ea5d097f8fe46f_349    129364  Phytophthora rubi  Unique taxonomy match
    intermediate/DNA10MIX_undiluted.onebp.tsv:d8613e80b8803b13f7ea5d097f8fe46f_271     129364  Phytophthora rubi  Unique taxonomy match

The summary reports would also answer this particular question, but this kind
of search can be useful for exploring specific questions.
