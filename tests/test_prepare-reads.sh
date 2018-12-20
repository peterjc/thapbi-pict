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
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.prepared.fasta` -ne "892" ]; then echo "Wrong FASTA output count"; false; fi

rm -rf $TMP/DNAMIX_S95_L001.prepared.fasta
thapbi_pict prepare-reads -o $TMP tests/reads/DNAMIX_S95_L001_*.fastq.gz -a 100
if [ `grep -c "^>" $TMP/DNAMIX_S95_L001.prepared.fasta` -ne "7" ]; then echo "Wrong FASTA output count"; false; fi

echo "$0 passed"
