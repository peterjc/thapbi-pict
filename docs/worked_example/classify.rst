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
#sequence-name                            taxid          genus-species:...                                    note
========================================= ============== ==================================================== ====
``a559aa4d00a28f11b83012e762391259_2303`` 221518         *Phytophthora pseudosyringae*                        ...
``140ccd03a87b423a1f06521f08131464_724``  78237          *Phytophthora gonapodyides*                          ...
``868e1ad838c7ec587dfd05b9dd4556ec_339``  78237          *Phytophthora gonapodyides*                          ...
``742f1f7a934f2df075be6f2eea756fc9_210``  164328         *Phytophthora ramorum*                               ...
``7f27d3a8f7150e0ee7ad64073e6da6b5_193``  164328         *Phytophthora ramorum*                               ...
``eaf42569c8b95c8bf4f9bf1b65a96ce4_183``  53983; 2056922 *Phytophthora cambivora*; *Phytophthora x cambivora* ...
``ffb8fbb83fa26a101c2fddf2af13cf95_167``  631361         *Phytophthora austrocedri*                           ...
``af3654932ad7a06c5f4af3c738706c76_114``  78237          *Phytophthora gonapodyides*                          ...
========================================= ============== ==================================================== ====

This says most of the unique sequences here have been assigned a single unique
*Phytophthora* species, except for ``eaf42569c8b95c8bf4f9bf1b65a96ce4`` (found
in 183 reads for this sample) which has matched *Phytophthora cambivora* (NCBI
taxid 53983) and close relative *Phytophthora x cambivora* (NCBI taxid
2056922).

If you are familiar with the command line search tool ``grep`` and the regular
expression syntax, you should find the format of these intermediate TSV files
lends itself to some simple searches. For example, you could see which samples
had matches to *Phytophthora rubi* using ``grep`` twice as follows (exclude
header lines, then find species):

.. code:: console

    $ grep -v "^#" intermediate/*.tsv | grep "Phytophthora rubi"
    intermediate/DNA10MIX_bycopynumber.onebp.tsv:2ba87367bdbb87cc37521bed773ffa37_285  129364  Phytophthora rubi  Unique taxonomy match
    intermediate/DNA10MIX_diluted25x.onebp.tsv:2ba87367bdbb87cc37521bed773ffa37_363    129364  Phytophthora rubi  Unique taxonomy match
    intermediate/DNA10MIX_undiluted.onebp.tsv:2ba87367bdbb87cc37521bed773ffa37_274     129364  Phytophthora rubi  Unique taxonomy match

The summary reports would also answer this paricular question, but this kind
of search can be useful for exploring specific questions.
