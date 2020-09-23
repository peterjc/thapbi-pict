.. _negative_controls:

Negative Controls
=================

Any negative control sample is not expected to contain any of the target
sequences.

On a typical 96-well plate of PCR products which will go on to be multiplexed
for Illumina MiSeq sequencing, most of the samples are biological - but some
should be negative controls (e.g. PCR blanks, or synthetic sequences).
The presence of biological sequence reads in the negative control samples is
indicative of some kind of cross contamination.

Minimum Abundance Threshold
---------------------------

During the read preparation step, the ``-n`` or ``--negctrls`` argument gives
the filename stems of any negative controls. For example, with a consistent
prefix you might use something like ``-n CTRL*.fastq.gz`` for this. The
negative control samples are processed first, and if any biological sequences
are found above the minimum abundance threshold, the threshold is raised to
that level for the other samples on that plate. This assumes if you have
multiple 96-well plates, or other logical groups, their raw FASTQ files are
separated into a separate sub-folder per plate.

For example, if running with the default minimum abundance threshold of 100,
and a negative control contains a biological sequence at abundance 136, then
the threshold for the non-control samples in that folder is raised to 136.

Synthetic Controls
------------------

Four synthetic sequences were designed for Phyto-Threats project which funded
THAPBI PICT. These were of the typical expected ITS1 fragment length and base
content, had the typical fixed 32bp header, but were otherwise shuffled with
no biological meaning (avoiding any secondary structure forming). They were
synthesised using `Integrated DNA Technologies gBlocks Gene Fragments
<https://www.idtdna.com/pages/products/genes-and-gene-fragments/double-stranded-dna-fragments/gblocks-gene-fragments>`_.

Our 96-well PCR plates included multiple samples which were known ratios of
these synthetic sequences, rather than environmental DNA.

The tool needed a way to distinguish biological marker sequences (for which
we wanted to make as few assumptions as possible) from the synthetic ones
(where the template sequences were known, subject only to PCR noise).
This is implemented by providing a HMMER3 format Hidden Markov Model of the
synthetic reads during the read preparation step via the ``--hmm`` argument.
Any sequence matching the specified HMM is treated as a synthetic control,
anything not matching is assumed to be biological sequence. If you are using
blank PCR controls (or no negative controls at all), the default HMM will not
match anything, but using ``--hmm ""`` to explicitly skip this step would be
slightly faster.

Conversely, the presence of the synthetic controls in any of the biological
samples is also problematic. To that end, the synthetic control sequences are
also included in the default database and can be matched by the classifier,
and appear in the reports.
