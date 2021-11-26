.. _negative_controls:

Negative Controls
=================

Any negative control sample is not expected to contain any of the target
sequences (although it may contain spike-in synthetic control sequences).

On a typical 96-well plate of PCR products which will go on to be multiplexed
for Illumina MiSeq sequencing, most of the samples are biological - but some
should be negative controls (e.g. PCR blanks, or synthetic sequences).
The presence of biological sequence reads in the negative control samples is
indicative of some kind of cross contamination. Likewise, any reads from your
synthetic sequences in the biological samples are a warning sign.

Minimum Abundance Threshold
---------------------------

During the read preparation step, the ``-n`` or ``--negctrls`` argument gives
the filename stems of any negative controls. For example, with a consistent
prefix you might use something like ``-n CTRL*.fastq.gz`` for this. Negative
control samples are processed first, and if any non-spike-in sequences are
found above the minimum abundance threshold, the threshold is raised to that
level for the other samples on that plate. This assumes if you have multiple
96-well plates, or other logical groups, their raw FASTQ files are separated
into a separate sub-folder per plate.

For example, if running with the default minimum abundance threshold of 100,
and a negative control contains a biological sequence at abundance 136, then
the threshold for the non-control samples in that folder is raised to 136.

If you have no spike-in controls, then any sequences in the negative controls
can raise the threshold - regardless of what they may or may not match in the
reference database.

Currently the controls are not used to adjust the percentage abundance
threshold, only the absolite abundance threshold.

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
samples are considered to be spike-ins and ignored in the threshold
calculations above. While the PCR noise is typically just a few base pair
changes, we also found large deletions relatively common. The matching is
therefore relaxed, currently based on *k*-mer content.

Conversely, the presence of the synthetic controls in any of the biological
samples is also problematic. Since our synthetic control sequences are in
the default database, they and can be matched by the chosen classifier, and
appear in the reports.
