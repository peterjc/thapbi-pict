#!/bin/bash
IFS=$'\n\t'
set -eux

# Note not using "set -o pipefile" as want to use that
# with grep to check error messages

export TMP=${TMP:-/tmp}

echo "Checking prepare-reads"
thapbi_pict prepare-reads 2>&1 | grep "the following arguments are required"

# Try a real example
rm -rf $TMP/DNAMIX_S95_L001.prepared.fasta
thapbi_pict prepare-reads -o $TMP tests/reads/DNAMIX_S95_L001_*.fastq.gz -a 0
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.prepared.fasta` -ne "921" ]; then echo "Wrong FASTA output count"; false; fi

rm -rf $TMP/DNAMIX_S95_L001.prepared.fasta
thapbi_pict prepare-reads -o $TMP tests/reads/DNAMIX_S95_L001_*.fastq.gz -a 5
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.prepared.fasta` -ne "27" ]; then echo "Wrong FASTA output count"; false; fi

echo "Generating mock control file"
# Using just 50 real reads (50 * 4 = 200 lines)
zcat tests/reads/DNAMIX_S95_L001_R1_001.fastq.gz | head -n 200 > $TMP/MOCK_CONTROL_R1.fastq
zcat tests/reads/DNAMIX_S95_L001_R2_001.fastq.gz | head -n 200 > $TMP/MOCK_CONTROL_R2.fastq

rm -rf $TMP/DNAMIX_S95_L001.prepared.fasta
rm -rf $TMP/MOCK_CONTROL.prepared.fasta
# Starting low threshold, should be increased to 19, so get new output count...
thapbi_pict prepare-reads -o $TMP tests/reads/DNAMIX_S95_L001_*.fastq.gz -a 5 -c $TMP/MOCK_CONTROL_R?.fastq
if [ `grep -c "^>" $TMP/MOCK_CONTROL.prepared.fasta` -ne "16" ]; then echo "Wrong FASTA control output count"; false; fi
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.prepared.fasta` -ne "9" ]; then echo "Wrong FASTA output count"; false; fi

rm -rf $TMP/DNAMIX_S95_L001.prepared.fasta
thapbi_pict prepare-reads -o $TMP tests/reads/DNAMIX_S95_L001_*.fastq.gz -a 100
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.prepared.fasta` -ne "7" ]; then echo "Wrong FASTA output count"; false; fi

echo "$0 passed"
