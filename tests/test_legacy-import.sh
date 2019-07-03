#!/bin/bash

# Copyright 2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eux
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "Checking legacy-import"
thapbi_pict legacy-import 2>&1 | grep "the following arguments are required"
thapbi_pict legacy-import -d "sqlite:///:memory:" -i database/legacy/Phytophthora_ITS_database_v0.004.fasta 2>&1 | grep "Taxonomy table empty"
set -o pipefail

export DB=$TMP/dup_seqs.sqlite
rm -rf $DB
thapbi_pict legacy-import -x -d $DB -i tests/legacy-import/dup_seqs.fasta
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_source;"` -ne "8" ]; then echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "2" ]; then echo "Wrong its1_sequence count"; false; fi

export DB=database/legacy/database_lax.sqlite
rm -rf $DB
thapbi_pict legacy-import -x -d $DB -i database/legacy/database.fasta

thapbi_pict legacy-import -x -d "sqlite:///:memory:" -i database/legacy/Phytophthora_ITS_database_v0.004.fasta

rm -rf database/legacy/Phytophthora_ITS_database_v0.005.sqlite
thapbi_pict legacy-import -x -d "database/legacy/Phytophthora_ITS_database_v0.005.sqlite" -i database/legacy/Phytophthora_ITS_database_v0.005.fasta
# We use this DB later...

export DB=$TMP/legacy_004_and_005_lax.sqlite
rm -rf $DB
thapbi_pict legacy-import -x -d $DB -i database/legacy/Phytophthora_ITS_database_v0.004.fasta -n "Legacy DB v0.004"
thapbi_pict legacy-import -x -d $DB -i database/legacy/Phytophthora_ITS_database_v0.005.fasta -n "Legacy DB v0.005"
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "2" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_source;"` -ne "382" ]; then echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "176" ]; then echo "Wrong its1_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "165" ]; then echo "Wrong taxonomy count"; false; fi
if [ `thapbi_pict dump -d $DB -f fasta | grep -c "^>"` -ne "382" ]; then echo "Wrong FASTA record count"; false; fi
if [ `thapbi_pict dump -d $DB | grep "synthetic" -c` -ne 4 ]; then echo "Missing four synthetic controls"; false; fi


if [ ! -f "new_taxdump_2018-12-01.zip" ]; then curl -L -O "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/new_taxdump_2018-12-01.zip"; fi
if [ ! -d "new_taxdump_2018-12-01" ]; then unzip new_taxdump_2018-12-01.zip -d new_taxdump_2018-12-01; fi

# Now test with species name validation, load with Phytophthora
# Then add 172 entries from v5, then 170 entries from v4
export DB=$TMP/legacy_004_and_005_validated.sqlite
rm -rf $DB
thapbi_pict load-tax -d $DB -t new_taxdump_2018-12-01 -a 4783
if [ `sqlite3 $DB "SELECT COUNT(DISTINCT genus) FROM taxonomy;"` -ne "1" ]; then echo "Wrong genus count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(DISTINCT species) FROM taxonomy;"` -ne "252" ]; then echo "Wrong species count"; false; fi
thapbi_pict legacy-import -d $DB -i database/legacy/Phytophthora_ITS_database_v0.005.fasta -n "Legacy DB v0.005"
if [ `sqlite3 $DB "SELECT COUNT(DISTINCT species) FROM taxonomy;"` -ne "252" ]; then echo "Wrong species count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_source;"` -ne "175" ]; then echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "159" ]; then echo "Wrong its1_sequence count"; false; fi
thapbi_pict legacy-import -d $DB -i database/legacy/Phytophthora_ITS_database_v0.004.fasta -n "Legacy DB v0.004"
if [ `sqlite3 $DB "SELECT COUNT(DISTINCT species) FROM taxonomy;"` -ne "252" ]; then echo "Wrong species count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "2" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_source;"` -ne "348" ]; then echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "159" ]; then echo "Wrong its1_sequence count"; false; fi
if [ `thapbi_pict dump -d $DB -f fasta | grep -c "^>"` -ne "348" ]; then echo "Wrong FASTA record count"; false; fi

#thapbi_pict dump 2>&1 | grep "the following arguments are required"
thapbi_pict dump -d database/legacy/Phytophthora_ITS_database_v0.005.sqlite -o /dev/null
thapbi_pict dump -d "sqlite:///database/legacy/Phytophthora_ITS_database_v0.005.sqlite" -o /dev/null -c 8a,8b

export DB=database/legacy/database_lax.sqlite
thapbi_pict dump -d $DB -o /dev/null -c -
thapbi_pict dump -d $DB -o /dev/null -g Phytophthora
thapbi_pict dump -d $DB -o /dev/null -g Phytophthora -s "ilicis, sp. aff. meadii"

set +o pipefail
thapbi_pict dump -d $DB -o /dev/null -s "ambiguous" 2>&1 | grep "requires a single genus"
thapbi_pict dump -d $DB -o /dev/null -g "Phytophthora" -s "ambiguous" 2>&1 | grep "not in database"
thapbi_pict dump -d $DB -o /dev/null -g "Phytopthora" 2>&1 | grep "not in database"
thapbi_pict dump -d $DB -o /dev/null -c "123" 2>&1 | grep "not in database"
set -o pipefail

echo "$0 - test_legacy-import.sh passed"
