#!/bin/bash
IFS=$'\n\t'
set -eux

# Note not using "set -o pipefile" as want to use that
# with grep to check error messages

export TMP=${TMP:-/tmp}

echo "Checking classify-reads"
thapbi_pict classify-reads 2>&1 | grep "the following arguments are required"
thapbi_pict classify-reads -d "sqlite:///:memory:" hypothetical_example.fasta 2>&1 | grep "cannot classify anything"

export DB=$TMP/legacy_004_and_005_validated.sqlite
if [ ! -f $DB ]; then echo "Run test_legacy-import.sh to setup test DB"; false; fi

rm -rf database/legacy/*.identity-reads.tsv
rm -rf database/legacy/*.identity-tax.tsv

# Passing one filename; default output dir:
thapbi_pict classify-reads -m identity -d $DB database/legacy/database.fasta
# Now fails as we expect reads to have been prepared and trimmed with HMM
if [ `wc -l database/legacy/database.identity-reads.tsv` -ne "0" ]; then echo "Expected no matches"; false; fi

if [ ! -f $TMP/DNAMIX_S95_L001.prepared.fasta]; then echo "Run test_prepare-reads.sh to setup test input"; false; fi
rm -rf $TMP/DNAMIX_S95_L001.prepared.swarm-reads.tsv
rm -rf $TMP/DNAMIX_S95_L001.prepared.swarm-tax.tsv

thapbi_pict classify-reads -m swarm -d $DB $TMP/DNAMIX_S95_L001.prepared.fasta
cut -f 5 $TMP/DNAMIX_S95_L001.prepared.swarm-reads.tsv | sort | uniq -c
# grep -c Phytophthora $TMP/DNAMIX_S95_L001.prepared.swarm-tax.tsv

# Passing one directory name (should get all three FASTA files):
rm -rf $TMP/legacy/*.identity-*.tsv
mkdir -p $TMP/legacy
thapbi_pict classify-reads -m identity -d $DB database/legacy/ -o $TMP/legacy
ls -1 $TMP/legacy/*.identity-*.tsv  # Should be 6 pairs

echo "$0 passed"
