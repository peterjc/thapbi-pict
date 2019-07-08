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

You should find 122 small FASTQ files in the ``intermediate/`` folder (or you
can get these from the compressed file as described above). Note this is
robust to being interupted and restarted (e.g. a job might time out on the
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
  do not overlap are discarded, for example from unexpectedly long fragements,
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

The sequence entries in the FASTA file are named ``<checksum>_<abundance>``.
Here ``<checksum>`` is the `MD5 checksum <https://en.wikipedia.org/wiki/MD5>`_
of the sequence, and this is used as a unique shorthand. It is a 32 character
string of the digits ``0`` to ``9`` and lower cases letters ``a`` to ``f``
inclusive, like ``a559aa4d00a28f11b83012e762391259``. These MD5 checksums are
used later in the pipeline, including in reports. The ``<abundance>`` is just
an integer, the number of paired reads which after processing had this unique
sequence.

The description entry in the FASTA file is currently just the name of the HMM
it matched, allowing us to distinguish the biological ITS1 sequences from the
synthetic controls.

Finally, the sequence in the FASTA file is written as a single line in upper
case. With standard FASTA line wrapping at 60 or 80 characters, the ITS1
sequences would need a few lines each. However, they are still short enough
that having them on one line without line breaks is no hardship - and it is
*extremely* helpful for simple tasks like using ``grep`` to look for a
particular sequence fragment at the command line.

For example,

.. code:: console

    $ cat intermediate/Site_1_sample_1.fasta
    >a559aa4d00a28f11b83012e762391259_2303 phytophthora_its1
    CCACACCTAAAAAACTTTCCACGTGAACTGTATCGAACAACTAGTTGGGGGTCTTGTTTGGCGTGCGGCTGCTTCGGTAGCTGCTGCTAGGCGAGCCCTATCACGGCGAGCGTTTGGACTTCGGTCTGAGCTAGTAGCTATTTTTTAAACCCATTCTTTAATACTGATTATACT
    >140ccd03a87b423a1f06521f08131464_724 phytophthora_its1
    CCACACCTAAAAAAACTTTCCACGTGAACCGTATCAACCCCTATAATTTGGGGGCTTGCTCGGCGGCGTGTGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCGGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAGCTAGTAGCTATCAATTTTAAACCCTTTCTTAAATACTGAACATACT
    >868e1ad838c7ec587dfd05b9dd4556ec_339 phytophthora_its1
    CCACACCTAAAAAAAACTTTCCACGTGAACCGTATCAACCCCTATAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCGGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAGCTAGTAGCTATCAATTTTAAACCCTTTCTTAAATACTGAACATACT
    >742f1f7a934f2df075be6f2eea756fc9_210 phytophthora_its1
    CCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAACCGTTAGTTGGGGGCTTCTGTTCGGCTGGCTTCGGCTGGCTGGGCGGCGGCTCTATCATGGCGAGCGCTTGAGCCTTCGGGTCTGAGCTAGTAGCCCACTTTTTAAACCCATTCCTAAATACTGAATATACT
    >7f27d3a8f7150e0ee7ad64073e6da6b5_193 phytophthora_its1
    CCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAACCCTTAGTTGGGGGCTTCTGTTCGGCTGGCTTCGGCTGGCTGGGCGGCGGCTCTATCATGGCGAGCGCTTGAGCCTTCGGGTCTGAGCTAGTAGCCCACTTTTTAAACCCATTCCTAAATACTGAATATACT
    >eaf42569c8b95c8bf4f9bf1b65a96ce4_183 phytophthora_its1
    CCACACCTAAAAAACTTTCCACGTGAACCGTATCAACCCACTTAGTTGGGGGCTAGTCCCGGCGGCTGGCTGTCGATGTCAAAGTTGACGGCTGCTGCTGTGTGTCGGGCCCTATCATGGCGAGCGTTTGGGTCCCTCTCGGGGGAACTGAGCCAGTAGCCCTTATTTTTTAAACCCATTCTTGAATACTGAATATACT
    >ffb8fbb83fa26a101c2fddf2af13cf95_167 phytophthora_its1
    CCACACCTAAAAAACTTTCCACGTGAACCGTATCAAAATCCTTTTATTGGGGGCTTCTGTCTGGTCTGGCTTCGGCTGGTCTGGGTGGCGGCTCTATCATGGTGACCGCTCTGGGCTTCGGCTTGGAGTTAGTAGCCCACTTTTTAAACCCATTCTTAATTACTGAACATACT
    >af3654932ad7a06c5f4af3c738706c76_114 phytophthora_its1
    CCACACCTAAAAAAACTTTCCACGTGAACCGTATCAACCCCTATAATTTGGGGGCTTGCTCGGCGGCGTGCGTGCTGGCCTGTAATGGGTCGGCGTGCTGCTGCTGGGCGGGCTCTATCATGGGCGAGCGTTTGGGCTTCGGCTCGAGCTAGTAGCTATCAATTTTAAACCCTTTCTTAAATACTGAACATACT

We see this sample had eight unique sequences accepted, all matched the ITS1
HMM (happily none match the synthetic controls). The most common had MD5
checksum ``a559aa4d00a28f11b83012e762391259`` and was seen in 2303 reads.

You could easily find out which other samples had this unique sequence using
the command line search tool ``grep`` as follows:

.. code:: console

    $ grep a559aa4d00a28f11b83012e762391259 intermediate/*.fasta
    ...

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

For example, to lower the threshold from the default to 50, you could use:

.. code:: console

    $ rm -rf intermediate/*.fasta
    $ thapbi_pict prepare-reads -i raw_data/ -o intermediate/ -a 50
    ...

.. WARNING::

   By default ``thapbi_pict prepare-reads`` and ``thapbi_pict pipeline`` will
   reuse existing intermediate FASTA files, so you must explicitly delete any
   old FASTA files before the new abundance threshold will have any effect.

.. WARNING::

    Setting the abundance threhold low (say under 50) risks letting background
    contamination through into the results. Do not do this without strong
    justification (e.g. look at suitable controls over multiple plates from
    your own laboratory procedure).

.. WARNING::

    Setting the abundance threshold *very* low (under 10) has the additional
    problem that the number of unique sequences accepted will increase many
    times over. This will *dramatically* slow down the rest of the analysis.
    This is only advised for investigating single samples.

For the woody host data, each plate had a negative control sample which should
contain no ITS1 sequences. We can specify the negative controls with ``-n`` or
``--negctrls`` by entering the four FASTQ filenames in full, but since they
have a common prefix we can use a simple wild card:

.. code:: console

    $ thapbi_pict prepare-reads -i raw_data/ -o intermediate/ -n raw_data/NEGATIVE*.fastq.gz
    ...

For this sample data, happily neither of the negative controls have any ITS1
present above the default threshold, so this would have no effect.

For the THAPBI project we now run each 96-well PCR plate with multiple
negative controls. Rather than a simple blank, these include a known mixture
of synthetic sequences of the same length, same nucelotide composition, and
also same di-nucleotide composition as real *Phytophthora* ITS1. This means we
might have say 90 biological samples which should contain ITS1 but not the
synthetics controls, and 6 negative controls which should contain synthetic
controls but not ITS1. We then run ``thapbi_pict prepare-reads`` separately
for each plate, where any ITS1 contamination in the synthetic controls is
used to set a plate specific minimum abundance. This means we cannot run
``thapbi_pict pipeline`` on multiple plates at once (although we could run it
on each plate, we generally want to produce reports over multiple plates).
