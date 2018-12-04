#!/bin/bash
IFS=$'\n\t'
set -eux

# Note not using "set -o pipefile" as want to use that
# with grep to check error messages

export TMP=${TMP:-/tmp}

echo "Checking ncbi-import"
thapbi_pict ncbi-import 2>&1 | grep "the following arguments are required"

if [ ! -f $TMP/ncbi_sample.fasta ]; then esearch -db nucleotide -query "its1 AND Phytophthora[Organism] AND 150:800[Sequence Length]" | efetch -format fasta > $TMP/ncbi_sample.fasta; fi

grep -c "^>" $TMP/ncbi_sample.fasta

rm -rf $TMP/ncbi_sample.sqlite
thapbi_pict ncbi-import -d $TMP/ncbi_sample.sqlite $TMP/ncbi_sample.fasta

if [ `sqlite3 $TMP/ncbi_sample.sqlite "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $TMP/ncbi_sample.sqlite "SELECT COUNT(id) FROM taxonomy;"` -lt "100" ]; then echo "Taxonomy count too low"; false; fi

thapbi_pict dump 2>&1 | grep "the following arguments are required"
thapbi_pict dump -d $TMP/ncbi_sample.sqlite -o /dev/null

echo "$0 passed"
