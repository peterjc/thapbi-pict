#!/bin/bash

# Copyright 2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "======================"
echo "Checking curated-import"
echo "======================"
set -x
thapbi_pict curated-import 2>&1 | grep "the following arguments are required"
thapbi_pict curated-import -d "sqlite:///:memory:" -i database/controls.fasta 2>&1 | grep "Taxonomy table empty"
set -o pipefail

# Without the synthetic controls in the taxonomy
export DB=$TMP/contols_lax.sqlite
rm -rf $DB
thapbi_pict curated-import -d $DB -i database/controls.fasta -x
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_source;"` -ne "4" ]; then echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "4" ]; then echo "Wrong its1_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "4" ]; then echo "Wrong taxonomy count"; false; fi
if [ `sqlite3 $DB "SELECT DISTINCT genus, species FROM taxonomy;" | wc -l` -ne 4 ]; then echo "Wrong species count"; false; fi

export DB=$TMP/contols_strict.sqlite
rm -rf $DB
thapbi_pict load-tax -d $DB -t new_taxdump_2019-09-01
# Add the G-BLOCK synthetic controls
sqlite3 $DB "INSERT INTO taxonomy (ncbi_taxid, genus, species) VALUES (32630, 'synthetic', 'construct C1');"
sqlite3 $DB "INSERT INTO taxonomy (ncbi_taxid, genus, species) VALUES (32630, 'synthetic', 'construct C2');"
sqlite3 $DB "INSERT INTO taxonomy (ncbi_taxid, genus, species) VALUES (32630, 'synthetic', 'construct C3');"
sqlite3 $DB "INSERT INTO taxonomy (ncbi_taxid, genus, species) VALUES (32630, 'synthetic', 'construct C4');"
thapbi_pict curated-import -d $DB -i database/controls.fasta -v
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_source;"` -ne "4" ]; then echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "4" ]; then echo "Wrong its1_sequence count"; false; fi

echo "$0 - test_curated-import.sh passed"
