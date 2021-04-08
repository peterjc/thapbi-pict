.. _prepare_reads:

Preparing the sequence data
===========================

Running thapbi-pict prepare-reads
---------------------------------

Calling ``thapbi-pict prepare-reads`` is the first action done by the top
level ``thapbi_pict pipeline`` command.

.. code:: console

    $ thapbi_pict prepare-reads -h
    ...

Assuming you have the FASTQ files in ``raw_data/`` as described above:

.. code:: console

    $ thapbi_pict prepare-reads -i raw_data/ -o intermediate/
    ...

For each input FASTQ file pair ``raw_data/<sample_name>_R1.fastq.gz`` and
``raw_data/<sample_name>_R2.fastq.gz`` you should get a small FASTA file
``intermediate/<sample_name>.fasta``. In this case, there are multiple
replicates from each of 14 sample sites where the file name stem is
``Site_<N>_sample_<X>``, plus the controls.

.. code:: console

    $ ls -1 intermediate/*.fasta | wc -l
    122

You should find 122 small FASTA files in the ``intermediate/`` folder (or you
can get these from the compressed file as described above). Note this is
robust to being interrupted and restarted (e.g. a job might time out on the
cluster).

.. WARNING::

    So far this example omits a key consideration - telling the tool which
    samples are negative controls, and/or manually setting the minimum read
    abundance. See below.

Intermediate FASTA files
------------------------

What the prepare command does can be briefly summarised as follows:

* Quality trim the FASTQ reads (pairs where either read becomes too short are
  discarded).
* Merge the overlapping paired FASTQ reads into single sequences (pairs which
  do not overlap are discarded, for example from unexpectedly long fragments,
  or not enough left after quality trimming).
* Primer trim (reads without both primers are discarded).
* Convert into a non-redundant FASTA file, with the sequence name recording
  the abundance (discarding sequences of low abundance).
* Filter with Hidden Markov Models (HMMs) of ITS1 and our four synthetic
  controls (non-matching sequences are discarded).

For each input ``<sample_name>_R1.fastq.gz`` and ``<sample_name>_R2.fastq.gz``
FASTQ pair we get a single *much* smaller FASTA file ``<sample_name>.fasta``.

.. WARNING::

   The intermediate FASTA files can legitimately be empty when no sequences
   passed the thresholds. This can happen when a PCR failed, and is expected
   to happen on blank negative controls.

.. WARNING::

   The intermediate FASTA files start with an atypical header made up of
   lines starting ``#``. Some tools need this to be removed, but others will
   accept this as valid FASTA format.

For example, here the header tells us this sample started with 6136 reads in
the paired FASTQ files, down to just 4185 after processing.

.. code:: console

    $ head -n 9 intermediate/Site_1_sample_1.fasta
    #left_primer:GAAGGTGAAGTCGTAACAAGG
    #right_primer:GCARRGACTTTCGTCCCYRC
    #raw_fastq:6136
    #flash:5900
    #cutadapt:5892
    #abundance:4185
    #threshold:100
    >2e4f0ed53888ed39a2aee6d6d8e02206_2272
    TTTCCGTAGGTGAACCTGCGGAAGGATCATTACCACACCTAAAAAACTTTCCACGTGAACTGTATCGAACAACTAGTTGG
    GGGTCTTGTTTGGCGTGCGGCTGCTTCGGTAGCTGCTGCTAGGCGAGCCCTATCACGGCGAGCGTTTGGACTTCGGTCTG
    AGCTAGTAGCTATTTTTTAAACCCATTCTTTAATACTGATTATACT

The sequence entries in the FASTA file are named ``<checksum>_<abundance>``.
Here ``<checksum>`` is the `MD5 checksum <https://en.wikipedia.org/wiki/MD5>`_
of the sequence, and this is used as a unique shorthand. It is a 32 character
string of the digits ``0`` to ``9`` and lower cases letters ``a`` to ``f``
inclusive, like ``a559aa4d00a28f11b83012e762391259``. These MD5 checksums are
used later in the pipeline, including in reports. The ``<abundance>`` is just
an integer, the number of paired reads which after processing had this unique
sequence.

The description entry in the FASTA file is currently just the name of any HMM
it matched, allowing us to distinguish biological marker sequences (no match)
from the synthetic controls (HMM match).

Finally, the sequence in the FASTA file is written as a single line in upper
case. With standard FASTA line wrapping at 60 or 80 characters, the ITS1
sequences would need a few lines each. However, they are still short enough
that having them on one line without line breaks is no hardship - and it is
*extremely* helpful for simple tasks like using ``grep`` to look for a
particular sequence fragment at the command line.

Note that for this documentation, the FASTA output has had the sequences line
wrapped at 80 characters.

.. code:: console

    $ grep "^>" intermediate/Site_1_sample_1.fasta
    >2e4f0ed53888ed39a2aee6d6d8e02206_2272
    >c1a720b2005f101a9858107545726123_716
    >96e0e2f0475bd1617a4b05e778bb04c9_331
    >fb30156d7f66c8abf91f9da230f4d19e_212
    >dcd6316eb77be50ee344fbeca6e005c7_194
    >972db44c016a166de86a2bacab3f4226_182
    >d9bc3879fdab3b4184c04bfbb5cf6afb_165
    >ed15fefb7a3655147115fc28a8d6d671_113

The final output has just eight unique sequences accepted, happily none of
which match the synthetic controls. The most common is listed first, and had
MD5 checksum ``2e4f0ed53888ed39a2aee6d6d8e02206`` and was seen in 2272 reads.

You could easily find out which other samples had this unique sequence using
the command line search tool ``grep`` as follows:

.. code:: console

    $ grep 2e4f0ed53888ed39a2aee6d6d8e02206 intermediate/*.fasta
    ...

Or, since we deliberately record the sequences without line wrapping, you
could use ``grep`` with the actual sequence instead (which might spot some
slightly longer entries as well).

You can also answer this example question from the read report produced later.

Abundance thresholds
--------------------

As you might gather from reading the command line help, there are two settings
to do with the minimum read abundance threshold, ``-a`` or ``--abundance``
(default 100), and ``-n`` or ``--negctrls`` for specifying negative controls
(default none).

If any negative controls are specified, those paired FASTQ files are processed
*first*, using the specified minimum abundance (default 100). If any of these
contained ITS1 sequences above the threshold, that higher number is used as
the minimum abundance threshold for the non-control samples. For example, say
one control had several ITS1 sequences with a maximum abundance of 124, and
another control had a maximum ITS1 abundance of 217, while the remaining
controls had no ITS1 sequence above the default level. In that case, the tool
would take maximum 217 as the abundance threshold for the non-control samples.

If you wished to lower the threshold from the default to 50, you could use:

.. code:: console

    $ rm -rf intermediate/*.fasta  # Are you sure?
    $ thapbi_pict prepare-reads -i raw_data/ -o intermediate/ -a 50
    ...

.. WARNING::

   By default ``thapbi_pict prepare-reads`` and ``thapbi_pict pipeline`` will
   reuse existing intermediate FASTA files, so you must explicitly delete any
   old FASTA files before the new abundance threshold will have any effect.

.. WARNING::

    Setting the abundance threshold low (say under 50) risks background
    contamination coming through into the results. Do not do this without
    strong justification (e.g. look at suitable controls over multiple plates
    from your own laboratory procedure).

.. WARNING::

    Setting the abundance threshold *very* low (under 10) has the additional
    problem that the number of unique sequences accepted will increase many
    times over. This will *dramatically* slow down the rest of the analysis.
    This is only advised for investigating single samples.

For the woody host data, each plate had a negative control sample which should
contain no ITS1 sequences. We can specify the negative controls with ``-n`` or
``--negctrls`` by entering the four FASTQ filenames in full, but since they
have a common prefix we can use a simple wildcard:

.. code:: console

    $ thapbi_pict prepare-reads -i raw_data/ -o intermediate/ -n raw_data/NEGATIVE*.fastq.gz
    ...

For this sample data, happily neither of the negative controls have any ITS1
present above the default threshold, so this would have no effect.

For the THAPBI Phyto-Threats project we now run each 96-well PCR plate with
multiple negative controls. Rather than a simple blank, these include a known
mixture of synthetic sequences of the same length, same nucelotide
composition, and also same di-nucleotide composition as real *Phytophthora*
ITS1. This means we might have say 90 biological samples which should contain
ITS1 but not the synthetics controls, and 6 negative controls which should
contain synthetic controls but not ITS1.

We therefore run ``thapbi_pict prepare-reads`` separately for each plate,
where any ITS1 contamination in the synthetic controls is used to set a plate
specific minimum abundance. This means we cannot run ``thapbi_pict pipeline``
on multiple plates at once (although we could run it on each plate, we
generally want to produce reports over multiple plates).
