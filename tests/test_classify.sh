#!/bin/bash
IFS=$'\n\t'
set -eux
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "Checking classify"
thapbi_pict classify 2>&1 | grep "the following arguments are required"
thapbi_pict classify -d "sqlite:///:memory:" hypothetical_example.fasta 2>&1 | grep "cannot classify anything"
set -o pipefail

export DB=$TMP/legacy_004_and_005_validated.sqlite
if [ ! -f $DB ]; then echo "Run test_legacy-import.sh to setup test DB"; false; fi

rm -rf database/legacy/*.identity.tsv

# Passing one filename; default output dir:
thapbi_pict classify -m identity -d $DB database/legacy/database.fasta
if [ "`grep -c -v '^#' database/legacy/database.identity.tsv`" -ne "`grep -c '^>' database/legacy/database.fasta`" ]; then echo "Expected one line per input seq"; false; fi

rm -rf $TMP/DNAMIX_S95_L001.identity.tsv
rm -rf $TMP/thapbi_onebp
rm -rf $TMP/thapbi_swarm
rm -rf $TMP/thapbi_blast
mkdir -p $TMP/thapbi_onebp
mkdir -p $TMP/thapbi_swarm
mkdir -p $TMP/thapbi_blast

# Explicitly setting output directory, would be here anyway:
thapbi_pict classify -m identity -d $DB tests/prepare-reads/DNAMIX_S95_L001.fasta -o $TMP/
thapbi_pict classify -m onebp -d $DB tests/prepare-reads/DNAMIX_S95_L001.fasta -o $TMP/thapbi_onebp
thapbi_pict classify -m blast -d $DB tests/prepare-reads/DNAMIX_S95_L001.fasta -o $TMP/thapbi_blast
thapbi_pict classify -m swarm -d $DB tests/prepare-reads/DNAMIX_S95_L001.fasta -o $TMP/thapbi_swarm
thapbi_pict classify -m swarmid -d $DB tests/prepare-reads/DNAMIX_S95_L001.fasta

# Passing one directory name (should get all three FASTA files):
rm -rf $TMP/legacy
mkdir -p $TMP/legacy
thapbi_pict classify -m identity -d $DB database/legacy/ -o $TMP/legacy
if [ "`ls -1 $TMP/legacy/*.identity.tsv | wc -l`" -ne "3" ]; then echo "Expected 3 files"; false; fi

# Test using sequences from a single isolate control,
rm -rf $TMP/P-infestans-T30-4.*.tsv
for M in identity onebp blast swarm swarmid; do
    echo "Checking single isolate control with $M"
    thapbi_pict classify -d $DB tests/classify/P-infestans-T30-4.fasta -o $TMP/ -m $M
    diff $TMP/P-infestans-T30-4.$M.tsv tests/classify/P-infestans-T30-4.$M.tsv
done

echo "$0 - test_classify.sh passed"
