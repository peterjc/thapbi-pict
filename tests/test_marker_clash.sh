#!/bin/bash

# Copyright 2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -euo pipefail

export TMP=${TMP:-/tmp/thapbi_pict}/marker_clash
rm -rf $TMP
mkdir -p $TMP

echo "======================="
echo "Checking marker clashes"
echo "======================="

export DB=$TMP/clashes.sqlite

set -x
# This also defines the marker name, primers, and length range for preparing reads
thapbi_pict import -d $DB -i database/controls.fasta -x \
    -k ITS1 -l GAAGGTGAAGTCGTAACAAGG -r GCARRGACTTTCGTCCCYRC \
    --minlen 100 --maxlen 1000

# Add a biological ITS1 sequence too
thapbi_pict import -d $DB -k ITS1 -x -i tests/marker_clash/Phytophthora_cinnamomi.fasta

# Now setup the clash, import the controls AGAIN under a different marker name
# File needs to have different MD5 checksum due to current import assumptions
head -n 6 database/controls.fasta | sed "s/synthetic/unreal/g" > $TMP/three_controls.fasta
thapbi_pict import -d $DB -i $TMP/three_controls.fasta -x \
    -k NotITS1 -l ATGCGATACTTGGTGTGAAT -r GACGCTTCTCCAGACTACAAT \
    --minlen 50 --maxlen 250

if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_definition;"` -ne "2" ]; then echo "Wrong marker count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "3" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM sequence_source;"` -ne "8" ]; then echo "Wrong sequence_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_sequence;"` -ne "5" ]; then echo "Wrong marker_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "8" ]; then echo "Wrong taxonomy count"; false; fi
if [ `sqlite3 $DB "SELECT DISTINCT genus, species FROM taxonomy;" | wc -l` -ne 8 ]; then echo "Wrong species count"; false; fi

# This should fail - it should spot the 3 marker level conflicts
set +o pipefail
thapbi_pict conflicts -d $DB -o $TMP/conflicts.tsv | true
set -o pipefail
diff $TMP/conflicts.tsv tests/marker_clash/conflicts.tsv

for M in identity onebp substr 1s2g 1s3g 1s4g 1s5g blast; do
    # This should work, making it explicit to use ITS1:
    rm -rf $TMP/*.tsv
    thapbi_pict classify -m $M  -d $DB -i tests/marker_clash/GBLOCK-example.fasta -o $TMP/ -k ITS1
    diff $TMP/GBLOCK-example.${M}.tsv tests/marker_clash/GBLOCK-example.ITS1.${M}.tsv

    # Doesn't yet, but this should ideally fail (mismatch between DB and FASTA header information):
    # For now it should check for synthetics 1 (not present here), 2 and 3 (but not 4 and not P. cinnamomi)
    rm -rf $TMP/*.tsv
    thapbi_pict classify -m $M -d $DB -i tests/marker_clash/GBLOCK-example.fasta -o $TMP/ -k NotITS1
    diff $TMP/GBLOCK-example.${M}.tsv tests/marker_clash/GBLOCK-example.NotITS1.tsv

    # This should fail (multiple markers and not told which one to use):
    set +o pipefail
    thapbi_pict classify -m $M -d $DB  -i tests/marker_clash/GBLOCK-example.fasta -o - 2>&1 | grep "DB has multiple amplicon markers"
    set -o pipefail
done

set +x

echo "===="
echo "Done"
echo "===="

echo "$0 - test_marker_clash.sh passed"
