.. _negative_controls:

Abundance & Negative Controls
=============================

Any negative control sample is not expected to contain any of the target
sequences (although it may contain spike-in synthetic control sequences).

On a typical 96-well plate of PCR products which will go on to be multiplexed
for Illumina MiSeq sequencing, most of the samples are biological - but some
should be negative controls (e.g. PCR blanks, or synthetic sequences).
The presence of biological sequence reads in the negative control samples is
indicative of some kind of cross contamination. Likewise, synthetic sequences
in the biological samples are a warning sign.

The tool implements both absolute and fractional abundance thresholds, which
can be specified at the command line.  Moreover, control samples can be used
to automatically raise the threshold for batches of samples. Simple negative
controls can be used to set an absolute abundance threshold, but to set the
fractional abundance threshold we need to be able to distinguish expected
sequences from unwanted ones. For this we require known spike-in control
sequences, which are clearly distinct from the biological markers.

Spike-in Controls
-----------------

Four synthetic sequences were designed for Phyto-Threats project which funded
THAPBI PICT. These were of the typical expected ITS1 fragment length and base
content, had the typical fixed 32bp header, but were otherwise shuffled with
no biological meaning (avoiding any secondary structure forming). They were
synthesised using `Integrated DNA Technologies gBlocks Gene Fragments
<https://www.idtdna.com/pages/products/genes-and-gene-fragments/double-stranded-dna-fragments/gblocks-gene-fragments>`_.

Our 96-well PCR plates included multiple control samples which were known
ratios of these synthetic sequences, rather than environmental DNA.

The tool needs a way to distinguish biological marker sequences (for which
we wanted to make as few assumptions as possible) from the synthetic ones
(where the template sequences were known, subject only to PCR noise).

The spike-in controls are assumed to be in the database, by default under
the synthetic "genus" but that is configurable. Similar sequences in the
samples are considered to be spike-ins. While the PCR noise is typically just
a few base pair changes, we also found large deletions relatively common. The
matching is therefore relaxed, currently based on *k*-mer content.

Conversely, the presence of the synthetic controls in any of the biological
samples is also problematic. Since our synthetic control sequences are in
the default database, they can be matched by the chosen classifier, and
appear in the reports.


Minimum Absolute Abundance Threshold
------------------------------------

The initial absolute abundance threshold is set at the command line with
``-a`` or ``--abundance`` giving an integer value. If your samples have
dramatically different read coverage, then the fractional abundance threshold
may be more appropriate (see below).

During the sample tally step, the ``-n`` or ``--negctrls`` argument gives
the sample filenames of any negative controls to use to potentially increase
the absolute abundance threshold (see below). If you have no spike-in
controls, then any sequences in these negative controls can raise the
threshold - regardless of what they may or may not match in the reference
database.

Minimum Fractional Abundance Threshold
--------------------------------------

The initial fractional abundance threshold is set at the command line with
``-f`` or ``--abundance-fraction`` to an floating point number between zero
and one, thus ``-f 0.001`` means 0.1%. This is a percentage of the reads
identified for each marker after merging the overlapping pairs and primer
matching.

During the read preparation step, the ``-y`` or ``--synctrls`` argument gives
the sample filenames of any synthetic controls to use to potentially increase
the absolute abundance threshold. This setting works in conjunction with the
database which must include the spike-in sequences under the genera specified
at the command lines with ``--synthetic`` (by default "synthetic").

Automatic thresholds
--------------------

Any control samples are processed first, before the biological samples, and
high read counts can raise the threshold to that level for the other samples
in that folder. This assumes if you have multiple 96-well plates, or other
logical groups, their raw FASTQ files are separated into a sub-folder per
plate.

Control samples given via ``-n`` can raise the absolute abundance threshold
(any synthetic spike-in reads are ignored for this), while controls given via
``-y`` can raise the fractional abundance threshold (but must have synthetic
spike-in reads in order to give a meaningful fraction).

For example, if running with the default minimum abundance threshold of 100
(set via ``-a 100``), and a negative control (set via
``-n raw_data/CTRL*.fastq.gz``) contains a non-spike-in (and thus presumably
biological) sequence at abundance 136, then the threshold for the non-control
samples in that folder is raised to 136.

Alternatively, you might have synthetic spike-in controls listed with
``-y raw_data/SPIKES*.fastq.gz`` and use ``-f 0.001`` to set a default
fractional abundance threshold of 0.1%. Suppose a control had 100,000 reads
for a marker passing the overlap merging and primer matching, of which 99,800
matched the spike-ins leaving 200 unwanted presumably biological reads, of
which the most abundant was at 176 copies. Then the fractional abundance
threshold would be raised slightly to 176 / 100000 = 0.00176 or 0.176%.

Note that a control sample can be used with *both* ``-n`` and ``-y``, so in
this second example that would *also* raise the absolute abundance threshold
to 176 reads.

Potentially a synthetic control sample can have unusually low read coverage,
meaning even a low absolute number of non-spike-in reads (at noise level)
would give a spuriously high inferred fractional abundance threshold. To guard
against this corner case, as a heuristic half the absolute abundance threshold
is applied to the synthetic control samples. Likewise, half of any fractional
abundance threshold is applied to the negative control samples, which guards
against spurious raising of the absolute threshold.

A similar problem would occur if you accidentally use ``-y`` on a sample
without any expected spike-in controls. This would suggest result in an overly
high fractional threshold, treated as an error.
