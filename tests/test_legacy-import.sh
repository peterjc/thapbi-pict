#!/bin/bash
IFS=$'\n\t'
set -eux

# Note not using "set -o pipefile" as want to use that
# with grep to check error messages

export TMP=${TMP:-/tmp}

echo "Checking legacy-import"
thapbi_pict legacy-import 2>&1 | grep "the following arguments are required"

rm -rf $TMP/dup_seqs.sqlite
thapbi_pict legacy-import -d $TMP/dup_seqs.sqlite tests/legacy-import/dup_seqs.fasta
if [ `sqlite3 $TMP/dup_seqs.sqlite "SELECT COUNT(id) FROM data_source;"` == "1" ]; then true; else echo "Wrong data_source count"; false; fi
if [ `sqlite3 $TMP/dup_seqs.sqlite "SELECT COUNT(id) FROM its1_source;"` == "8" ]; then true; else echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $TMP/dup_seqs.sqlite "SELECT COUNT(id) FROM its1_sequence;"` == "2" ]; then true; else echo "Wrong its1_sequence count"; false; fi

rm -rf database/legacy/database.sqlite
thapbi_pict legacy-import -d database/legacy/database.sqlite database/legacy/database.fasta

thapbi_pict legacy-import -d "sqlite:///:memory:" database/legacy/Phytophthora_ITS_database_v0.004.fasta

rm -rf database/legacy/Phytophthora_ITS_database_v0.005.sqlite
thapbi_pict legacy-import -d "database/legacy/Phytophthora_ITS_database_v0.005.sqlite" database/legacy/Phytophthora_ITS_database_v0.005.fasta
# We use this DB later...

rm -rf $TMP/legacy_004_and_005.sqlite
thapbi_pict legacy-import -d $TMP/legacy_004_and_005.sqlite database/legacy/Phytophthora_ITS_database_v0.004.fasta -n "Legacy DB v0.004"
thapbi_pict legacy-import -d $TMP/legacy_004_and_005.sqlite database/legacy/Phytophthora_ITS_database_v0.005.fasta -n "Legacy DB v0.005"
if [ `sqlite3 $TMP/legacy_004_and_005.sqlite "SELECT COUNT(id) FROM data_source;"` == "2" ]; then true; else echo "Wrong data_source count"; false; fi
if [ `sqlite3 $TMP/legacy_004_and_005.sqlite "SELECT COUNT(id) FROM its1_source;"` == "378" ]; then true; else echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $TMP/legacy_004_and_005.sqlite "SELECT COUNT(id) FROM its1_sequence;"` == "172" ]; then true; else echo "Wrong its1_sequence count"; false; fi

thapbi_pict dump 2>&1 | grep "the following arguments are required"
thapbi_pict dump -d database/legacy/Phytophthora_ITS_database_v0.005.sqlite -o /dev/null
thapbi_pict dump -d "sqlite:///database/legacy/Phytophthora_ITS_database_v0.005.sqlite" -o /dev/null -c 8a,8b
thapbi_pict dump -d database/legacy/database.sqlite -o /dev/null -c -

echo "$0 passed"
