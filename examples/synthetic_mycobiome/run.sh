#!/bin/bash
set -euo pipefail

echo "NOTE: Expected first time run time is XXX"
echo "(or about XXX minutes from the merged reads cache),"
echo "repeat runs take seconds just to regenerate reports."
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

    NAME=ITS2
    thapbi_pict import -d references.sqlite \
        -i references.fasta environment.fasta -s ";" -x \
        -k $NAME --left $LEFT --right $RIGHT -v
fi

echo ================
echo Running analysis
echo ================

# Using defaults -a 100 -f 0.001 (defaults) with sample counts after cutadapt
# from 221k to 1199k reads, gives thresholds 221 to 1199 (all above -a value).
#mkdir -p intermediate
#thapbi_pict pipeline -d references.sqlite \
#    -i raw_data/ expected/ \
#    --merged-cache tmp_merged/ \
#    -s intermediate/ -o summary/f4k \
#    -t metadata.tsv -x 1 -c 3,4,5 -v

echo "Running m6 plate with -f 0 -a 5 to compare with Figure 6"
mkdir -p intermediate_a5/
# NOT using negative controls, -n raw_data/SRR7109420_*.fastq.gz
# or -y raw_data/m6/SRR7109420_*.fastq.gz
thapbi_pict pipeline -d references.sqlite \
    -i raw_data/m6/ expected/ --merged-cache tmp_merged/ \
    -s intermediate_a5/ -o summary/a5 -a 5 -f 0 -m identity \
    -t metadata.tsv -x 1 -c 3,4,5

exit

# Using -a 50 -f 0.00025 (1/4000) with sample counts after cutadapt
# from 221k to 1199k reads, gives thresholds 55 to 300 (all above -a value).
#mkdir -p intermediate_f4k
#thapbi_pict pipeline -d references.sqlite \
#    -i raw_data/ expected/ -a 50 -f 0.00025 \
#    --merged-cache tmp_merged/ \
#    -s intermediate_f4k/ -o summary/f4k \
#    -t metadata.tsv -x 1 -c 3,4,5 -v

#mkdir -p intermediate_nc/
#
# Using defaults with negative controls, not -a 5 -f 0
#thapbi_pict pipeline -d references.sqlite \
#    -i raw_data/ expected/ -n raw_data/SRR7109420_*.fastq.gz \
#    --merged-cache tmp_merged/ \
#    -s intermediate_nc/ -o summary/negctrl \
#    -t metadata.tsv -x 1 -c 3,4,5 -v

mkdir -p intermediate_a5/ intermediate_p01/

# NOT using negative controls, -n raw_data/SRR7109420_*.fastq.gz
# or -y raw_data/m6/SRR7109420_*.fastq.gz
thapbi_pict pipeline -d references.sqlite \
    -i raw_data/ expected/ --merged-cache tmp_merged/ \
    -s intermediate_a5/ -o summary/a5 -a 5 -f 0 -m identity \
    -t metadata.tsv -x 1 -c 3,4,5

# Now using the control for m6 plate to set fractional abundance threshold,
# and sensible initial value -f 0.0001 (i.e. 0.01%) so m4A also clean:
thapbi_pict pipeline -d references.sqlite \
    -i raw_data/ expected/ --merged-cache tmp_merged/ \
    -y raw_data/m6/SRR7109420_*.fastq.gz \
    -s intermediate_p01/ -o summary/p01 -a 100 -f 0.0001 -m identity \
    -t metadata.tsv -x 1 -c 3,4,5


# Dataset has a single synthetic negative control SRR7109420, on the m6 run.
# This has 1199806 reads after cutadapt, and highest non-synthetic-spike-in
# sequence is seen 187 times, suggesting using an absolute threshold -a 187
# (via -n raw_data/m6/SRR7109420_*.fastq.gz), or a fractional threshold (via
# -y raw_data/m6/SRR7109420_*.fastq.gz) of 0.0115%, or -f 0.000115019, which
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
#
# Negative control highest non-synthetics are:
#
# >2ead414a3fc878d2ba5f9f4d3ca2ad80_187
# AACGCACCTTGCGCCCCTCGGTCATCCGAGGAGCACGCCCGTTCGAGTATCACGTTAACCTCTTGCTTACTGTAGTAATACGGTCTGCGAGGACTTGGGTGCTGCCGGCCCTGCGTCGGCTCGCCTCGAAATGCATTAGTGGGGACCTGCCCTGCGCGGTGTTGATAATTGTCTACGTCGTGTCGGTGTTGGCTCTTTGCTTCGAACTTCCCTAGGGAAACTTTTCAAAGTTTGATCTCGAATCGGGTGGGACTACCCGCTGAACTTAA
#
# This is 1bp away from OU498434.1 uncultured fungus.
#
# It is most common in m6-755-1 aka SRR7109412, and m6-736-2 aka SRR7109418
#
# Next is something which is perhaps synthetic in origin(?):
#
# >4271381c66214c66c3f2c86e04c09566_139
# AACGCACCTTGCGCTCCTTGGTATTCCGAGGAGCATGCCTGTTTGAGTGTCGTGAAATTATCAACTCTCTTACTTTATTGTAAGATGAGCTTGGACTTGGGGATTGCTGGTGTAAATCAGCTTCTCTTGAATGCATTAGCTGGAATTTAGTTCGCAGCATATATGCGGTGTGATAATGTCGTCACTGTGTAGCTCGAATTGTCTGGCTTCTAATCGTCCTTCACGGGACAATTTGATCATTTTGACCTCAAATCAGGTAGGACTACCCGCTGAACTTAA
#
# This is most common in m6-736-1 aka SRR7109414
#
# Third is another biological sample (fungi):
#
# >bd8a1a5a44bd54b26fa22eb962638777_34
# AACGCACCTTGCGCTCCTTGGTATTCCGAGGAGCATGCCTGTTTGAGTGTCGTGAAATTATCAACTCCCATTCTTTGTTGATTGGTGAGCTTGGACTTGGGGACTGCTGGTGCAAATCAGCTTCCTTTGAATGAATTAGCTGGAATTTGATTCGCAGCATATATGCGGTGTGATAATGTCGTCACTGTGTAGCTCGGATTGTCTGGCTTCTAATCGTCCTTCACGGACAACTTGATCTTTTGACCTCAAATCAGGTAGGACTACCCGCTGAACTTAA
#
# This is 1bp away from MT595066.1 Hyphodontia tropica.
#
# This is most common in m6-736-1 aka SRR7109414

echo ====
echo Done
echo ====
