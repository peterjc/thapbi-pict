#!/bin/bash
IFS=$'\n\t'
set -eux

# Note not using "set -o pipefile" as want to use that
# with grep to check error messages

export TMP=${TMP:-/tmp}

echo "Checking classify"
thapbi_pict classify 2>&1 | grep "the following arguments are required"
thapbi_pict classify -d "sqlite:///:memory:" hypothetical_example.fasta 2>&1 | grep "cannot classify anything"

export DB=$TMP/legacy_004_and_005_validated.sqlite
if [ ! -f $DB ]; then echo "Run test_legacy-import.sh to setup test DB"; false; fi

rm -rf database/legacy/*.identity.tsv

# Passing one filename; default output dir:
thapbi_pict classify -m identity -d $DB database/legacy/database.fasta
# Now fails as we expect reads to have been prepared and trimmed with HMM
if [ `wc -l database/legacy/database.identity.tsv` -ne "0" ]; then echo "Expected no matches"; false; fi

if [ ! -f $TMP/DNAMIX_S95_L001.fasta ]; then echo "Run test_prepare.sh to setup test input"; false; fi
rm -rf $TMP/DNAMIX_S95_L001.swarm.tsv
rm -rf $TMP/DNAMIX_S95_L001.identity.tsv

# Explicitly setting output directory, would be here anyway:
thapbi_pict classify -m identity -d $DB $TMP/DNAMIX_S95_L001.fasta -o $TMP/

thapbi_pict classify -m swarm -d $DB $TMP/DNAMIX_S95_L001.fasta
cut -f 5 $TMP/DNAMIX_S95_L001.swarm.tsv | sort | uniq -c

# Passing one directory name (should get all three FASTA files):
rm -rf $TMP/legacy/*.identity.tsv
mkdir -p $TMP/legacy
thapbi_pict classify -m identity -d $DB database/legacy/ -o $TMP/legacy
if [ `ls -1 $TMP/legacy/*.identity.tsv | wc -l` -ne `3` ]; then echo "Expected 3 files;" false; fi

echo "$0 passed"
