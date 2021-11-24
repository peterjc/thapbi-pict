#!/bin/bash

# Copyright 2019-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp/thapbi_pict}/curated_import
rm -rf $TMP
mkdir -p $TMP

export LEFT=GAAGGTGAAGTCGTAACAAGG
export RIGHT=GCARRGACTTTCGTCCCYRC

echo "======================="
echo "Checking curated-import"
echo "======================="
set -x
thapbi_pict import 2>&1 | grep "the following arguments are required"
thapbi_pict import -d "sqlite:///:memory:" -i database/controls.fasta 2>&1 | grep "Taxonomy table empty"
set -o pipefail

echo "Controls (lax mode, without the synthetic controls in the taxonomy)"
export DB=$TMP/contols_lax.sqlite
rm -rf $DB
thapbi_pict import -d $DB -k ITS1 -l $LEFT -r $RIGHT \
    -i database/controls.fasta -x --minlen 150 --maxlen 250
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM sequence_source;"` -ne "4" ]; then echo "Wrong sequence_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_sequence;"` -ne "4" ]; then echo "Wrong marker_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "4" ]; then echo "Wrong taxonomy count"; false; fi
if [ `sqlite3 $DB "SELECT DISTINCT genus, species FROM taxonomy;" | wc -l` -ne 4 ]; then echo "Wrong species count"; false; fi

echo "Curated ITS1 with taxdump"
# See also database/build_CURATED.sh
export DB=$TMP/curated.sqlite
rm -rf $DB
thapbi_pict load-tax -d $DB -t new_taxdump_2019-09-01
thapbi_pict import -d $DB -k ITS1 -l $LEFT -r $RIGHT --minlen 100 --maxlen 250 \
    -i database/Phytophthora_ITS1_curated.fasta -s ";"
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM sequence_source;"` -ne "260" ]; then echo "Wrong sequence_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_sequence;"` -ne "225" ]; then echo "Wrong marker_sequence count"; false; fi

echo "With duplicated sequences (and multi-entry records using Ctrl+A)"
# Note this file uses Ctrl+A separator:
#
# $ grep $'\001' tests/curated-import/dup_seqs.fasta
# >DDD Phytophthora fourEEE Phytophthora five
# >FFF Phytophthora alphaGGG Phytophthora beta
#
# Only 6 FASTA records, but two are double entries so want 8 here
export DB=$TMP/dup_seqs.sqlite
rm -rf $DB
thapbi_pict import -d $DB -k ITS1 -l $LEFT -r $RIGHT \
    -i tests/curated-import/dup_seqs.fasta -x -c ncbi -s $'\001'
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM sequence_source;"` -ne "8" ]; then echo "Wrong sequence_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_sequence;"` -ne "2" ]; then echo "Wrong marker_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "8" ]; then echo "Wrong taxonomy count"; false; fi
if [ `sqlite3 $DB "SELECT DISTINCT genus, species FROM taxonomy;" | wc -l` -ne 8 ]; then echo "Wrong species count"; false; fi


echo "With duplicated sequences (ignoring the multi-entry record naming)"
# Now try using wrong separator, semi-colon.
# Expect 6 entries, and two very silly species names!
export DB=$TMP/dup_seqs_bad.sqlite
rm -rf $DB
thapbi_pict import -d $DB -k ITS1 -l $LEFT -r $RIGHT \
    -i tests/curated-import/dup_seqs.fasta -x -c simple -s ";"
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_definition;"` -ne "1" ]; then echo "Wrong marker_definition count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM sequence_source;"` -ne "6" ]; then echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_sequence;"` -ne "2" ]; then echo "Wrong its1_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "6" ]; then echo "Wrong taxonomy count"; false; fi
if [ `sqlite3 $DB "SELECT DISTINCT genus, species FROM taxonomy;" | wc -l` -ne 6 ]; then echo "Wrong species count"; false; fi

echo "Redekar supplementary table 3"
export DB=$TMP/Redekar_et_al_2019_sup_table_3.sqlite
rm -rf $DB
# This uses the semi-colon separator
thapbi_pict import -d $DB -k ITS1 -l $LEFT -r $RIGHT \
    -i examples/recycled_water/Redekar_et_al_2019_sup_table_3.fasta -x -s ";"
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_definition;"` -ne "1" ]; then echo "Wrong marker_definition count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM sequence_source;"` -ne "1451" ]; then echo "Wrong sequence_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_sequence;"` -ne "838" ]; then echo "Wrong marker_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "423" ]; then echo "Wrong taxonomy count"; false; fi
if [ `sqlite3 $DB "SELECT DISTINCT genus, species FROM taxonomy;" | wc -l` -ne "423" ]; then echo "Wrong species count"; false; fi

echo "Example based on recycled_water worked example"
export DB=$TMP/recycled_water.sqlite
rm -rf $DB
# Loading NCBI taxonomy for handling synonyms,
# using same version as rest of tests - not as per the worked example:
thapbi_pict load-tax -d $DB -t new_taxdump_2019-09-01/
# Using -x / --lax (does not insist on taxonomy match)
# Adding 32bp conserved TTTCCGTAGGTGAACCTGCGGAAGGATCATTA to left primer
# Not giving primers, sequences are already trimmed
thapbi_pict import -x -s ";" -d $DB \
    -i examples/recycled_water/Redekar_et_al_2019_sup_table_3.fasta \
    --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA \
    --right AGCGTTCTTCATCGATGTGC --marker ITS1-long
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_definition;"` -ne "1" ]; then echo "Wrong marker_definition count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM sequence_source;"` -ne "1451" ]; then echo "Wrong sequence_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_sequence;"` -ne "838" ]; then echo "Wrong marker_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "3008" ]; then echo "Wrong taxonomy count"; false; fi
if [ `sqlite3 $DB "SELECT DISTINCT genus, species FROM taxonomy;" | wc -l` -ne "3008" ]; then echo "Wrong species count"; false; fi

echo "$0 - test_curated-import.sh passed"
