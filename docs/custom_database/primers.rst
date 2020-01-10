Pipeline with custom primers
============================

Running thapbi-pict pipeline
----------------------------

We first ran the pipeline command with
`:ref:default settings <recycle_water_defaults>`, but now we will give the
actual primers used. Again we assume you have setup the FASTQ files in
``raw_data/``, now run the pipeline command as follows where we set the
left and right primer sequences specifically:

.. code:: console

    $ mkdir intermediate_primers/ summary_primers/
    $ thapbi_pict pipeline \
      -i raw_data/ -o summary_primers/ -s intermediate_primers/ \
      --left GAAGGTGAAGTCGTAACAAGG --right AGCGTTCTTCATCGATGTGC \
      -r recycled-water-primers -t metadata.tsv -c 1,2,3,4,5,6,7
    ...
    $ ls -1 intermediate_primers/SRR*.fasta | wc -l
    384
    $ ls -1 intermediate_primers/SRR*.onebp.tsv | wc -l
    384
    $ ls -1 summary/summary_primers/thapbi-pict.*
    recycled-water-primers.reads.tsv
    recycled-water-primers.reads.xlsx
    recycled-water-primers.samples.tsv
    recycled-water-primers.samples.txt
    recycled-water-primers.edit-graph.xgmml

Again we used ``-r`` (or ``--report``) to specify a stem for the report
filenames, and the metadata is exactly as before.

Here we said the left primer is ``GAAGGTGAAGTCGTAACAAGG`` (same as the THAPBI
PICT default), but that the right primer is ``AGCGTTCTTCATCGATGTGC``. This has
reverse complement ``GCACATCGATGAAGAACGCT`` and is found about 60bp downstream
of the default right primer in *Phytophthora*, and should also match *Pythium*
and *Phytopythium* species.

i.e. We should now find the *Phytophthora* FASTA sequences extracted are about
60bp longer, and many more non-*Phytophthora* are accepted.


Results
-------

Will pick a couple of samples to compare and contrast with the first run...
