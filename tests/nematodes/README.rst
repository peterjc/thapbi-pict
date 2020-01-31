Nematode test dataset
=====================

Real world test data for the ``prepare-reads`` command with the ``--flip``
option.

This is a tiny sample from some currently unpublished nematode data from East
Africa, from a collaboration including John Jones, Danny Coyne, Laura Cortada,
Harrison Mburu, Jamie Orr, James Price and Kyriakos Varypatakis.

Here custom barcodes were added to the start and end of the PCR amplified
markers in a non-strand specific manner. This requires checking both strands
for the marker, as about half the reads are reversed in this dataset.

This sample is about 6% of a single barcode combination, extracted from the
original FASTQ using cutadapt for the demultiplexing, and then my own FASTQ
sampling script taking just 50 paired reads.

The test uses a minimum abundance threshold of 10 (as it is only a subsample),
which with the defaults gives only 9 copies of the most abundant sequence, and
thus no output. Adding --flip to the command brings this to 12 copies, and a
single accepted sequence.

From an NCBI BLAST online, the dominant sequence and handful of the lower
abundance variants are perfect or 99% identical matches to published Globodera
rostochiensis ITS1.
