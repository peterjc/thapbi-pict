Fungal mock community example
=============================

Based on the two amplicon sequencing libraries from this paper:

    Bakker (2018) A fungal mock community control for amplicon sequencing experiments
    https://doi.org/10.1111/1755-0998.12760
    https://www.ebi.ac.uk/ena/data/view/PRJNA377530

Both amplicon sequencing libraries were sequenced on an Illumina MiSeq
(multiplexed with other unrelated samples), giving 61 FASTQ paired files.

Amplicon library one amplified a small region of ITS1 using the BITS/B58S3
primer pair (ACCTGCGGARGGATC and GAGATCCRTTGYTRAAAGTT), as shown in the
paper's supplementary Table S4.

Amplicon library two amplified a larger region of ITS1 using the ITS1f/ITS2
primer pair (CTTGGTCATTTAGAGGAAGTAA GCTGCGTTCTTCATCGATGC), which includes the
first library’s target region entirely. Similar yields as per supplementary
Table S4 vs S5.

Additionally, amplicon library two amplified ITS2 using the ITS3‐KYO2 and
ITS4‐KYO3 primer pair (GATGAAGAACGYAGYRAA and CTBTTVCCKCTTCACTCG), with lower
yields as per supplementary Table S5 vs S6.

The ``metadata.tsv`` files for each library were based on the metadata
downloaded from the ENA, with some reformatting. The split into amplicon one
and two was based on supplementary Tables S4, S5 and S6 (for the mock
community samples) and reading the paper (for placing the negative controls).

Files ``ITS1.fasta`` and ``ITS2.fasta`` were extracted from supplementary
materials appendix S2, with the species name alone added to the FASTA titles
(for use with ``thapbi_pict curated-import``).

The provided ``negative_control.known.tsv`` and ``mock_community.known.tsv``
files lists the expected species in the negative controls (none) and the mock
community samples (the same 19 species, although not always in equal ratios).

Script ``setup.py`` will download the raw gzip compressed FASTQ files from
the ENA (124 files, 61 pairs, under 400MB in total), and setup appropriate
per-sample symlinks to the expected output in the ``expected/``
sub-directories for use with classifier assessment.

Script ``run.py`` will (assuming ``setup.py`` worked as expected) run the
THAPBI PICT pipeline three times (once on amplicon library one with a single
primer set; twice on amplicon library two with the alternative primers), each
time creating a classifier database from ``ITS1.fasta`` or ``ITS2.fasta`` as
appropriate.
