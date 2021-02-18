#!/bin/bash

# Copyright 2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -euo pipefail

export TMP=${TMP:-/tmp}
export `grep ^TAX= database/build_ITS1_DB.sh`

if [ ! -f "${TAX}.zip" ]; then curl -L -O "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/${TAX}.zip"; fi
if [ ! -d "database/${TAX}" ]; then unzip ${TAX}.zip nodes.dmp names.dmp -d database/${TAX}; fi

echo "===================================================="
echo "Loading expected DB and exporting as table and FASTA"
echo "===================================================="
set -x

export DB=$TMP/expected_ITS1_DB
rm -rf $DB.*
sqlite3 $DB.sqlite < database/ITS1_DB.sql

thapbi_pict dump -m -f fasta -d "$DB.sqlite" -o "$DB.fasta"
thapbi_pict dump -m -d "$DB.sqlite" -o "$DB.txt"

set +x
echo "========================="
echo "Rebuilding the default DB"
echo "========================="
set -x

cd database/
rm -rf ITS1_DB.sqlite ITS1_DB.sql ITS1_DB.txt ITS1_DB.fasta
./build_ITS1_DB.sh
cd ..

set +x
echo "========================================"
echo "Checking rebuilt DB matches expectations"
echo "========================================"
set -x

diff "$DB.fasta" database/ITS1_DB.fasta
diff "$DB.txt" database/ITS1_DB.txt
# git diff database/ITS1_DB.sql

echo "$0 - test_build_db.sh passed"
