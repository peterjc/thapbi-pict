#!/bin/bash
set -euo pipefail

echo "NOTE: Expected first time run time is just over 20 minutes"
echo "(or about 8 minutes from the merged reads cache),"
echo "repeat runs about a minute to regenerate reports."
echo

mkdir -p tmp_merged/ intermediate/ summary/

# Illumina primers via Supplemental_primers.xlsx
#
# Adapters
# ILL_LEFT=ACACTCTTTCCCTACACGACGCTCTTCCGATCT
# ILL_RIGHT=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT

# Illumina 5'  fITS7
LEFT=GTGARTCATCGAATCTTTG
# Illumina 3'  ITS4
RIGHT=TCCTCCGCTTATTGATATGC
#Gives RC GCATATCAATAAGCGGAGGA

if [ ! -f references.sqlite ]; then
    echo =================
    echo Creating database
    echo =================

    thapbi_pict import -d references.sqlite \
        -i references.fasta environment.fasta -s ";" -x \
        -k ITS2 --left $LEFT --right $RIGHT
fi

echo ==========================================
echo Running analysis excluding only singletons
echo ==========================================

# Dataset has a single synthetic negative control SRR7109420, on the m6 run.
# This has 1199806 reads after cutadapt, and highest non-synthetic-spike-in
# sequence is seen 187 times, suggesting using an absolute threshold -a 187
# (via -n raw_data/SRR7109420_*.fastq.gz), or a fractional threshold (via
# -y raw_data/SRR7109420_*.fastq.gz) of 0.0115%, or -f 0.000115019, which
# is 187/1199806.
#
# This means while default absolute threshold -a 100 is not high enough (since
# this dataset is very high coverage), the default fractional threshold 0.1%
# via -f 0.001 would be 1200 reads on this sample, neadly x10 higher than it
# needs to be to exclude all the apparent Illumina tag switching to/from
# biological to synthetic samples.
#
# Even using 0.01% as the fractional threshold would give a threshold of 120
# reads on the control, which is not quite enough to discard the top to
# non-synthetic entries which are consistent with but about 10x higher than
# the worst Illumina tag-swapping in the other direction (synthetics appearing
# in bioligical samples; see Figure 6).

echo "Very low threshold (excluding only singletons) to compare with Figure 6"
mkdir -p intermediate_a2/
# NOT using negative controls YET, -n raw_data/SRR7109420_*.fastq.gz
# or -y raw_data/m6/SRR7109420_*.fastq.gz
thapbi_pict pipeline -d references.sqlite \
    -i raw_data/ expected/ --merged-cache tmp_merged/ \
    -s intermediate_a2/ -o summary/a2 -a 2 -f 0 \
    -t metadata.tsv -x 1 -c 3,4

echo ========================================
echo Running analysis using synthetic control
echo ========================================

echo "Running m6 plate using synthetic control for percentage abundance threshold"
mkdir -p intermediate_ctrl/
thapbi_pict pipeline -d references.sqlite \
    -i raw_data/ expected/ --merged-cache tmp_merged/ \
    -y raw_data/SRR7109420_*.fastq.gz \
    -s intermediate_ctrl/ -o summary/ctrl -a 100 -f 0 \
    -t metadata.tsv -x 1 -c 3,4

echo ====
echo Done
echo ====
