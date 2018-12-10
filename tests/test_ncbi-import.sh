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

# Cannot use validation without having some taxonomy entries
thapbi_pict ncbi-import -d sqlite:///:memory: $TMP/ncbi_sample.fasta -s 2>&1 | grep "Taxonomy table empty"

export DB=$TMP/ncbi_sample.sqlite
rm -rf $DB
thapbi_pict ncbi-import -d $DB $TMP/ncbi_sample.fasta

if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -lt "100" ]; then echo "Taxonomy count too low"; false; fi
# Other values subject to change

thapbi_pict dump 2>&1 | grep "the following arguments are required"
thapbi_pict dump -d $DB -o /dev/null

if [ ! -f "new_taxdump_2018-12-01.zip" ]; then curl -L -O "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/new_taxdump_2018-12-01.zip"; fi
if [ ! -d "new_taxdump_2018-12-01" ]; then unzip new_taxdump_2018-12-01.zip -d new_taxdump_2018-12-01; fi

export DB=$TMP/ncbi_sample_validated.sqlite
rm -rf $DB
thapbi_pict load-tax -d $DB -t new_taxdump_2018-12-01 -a 4783
if [ `sqlite3 $DB "SELECT COUNT(DISTINCT genus) FROM taxonomy;"` -ne "1" ]; then echo "Wrong genus count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(DISTINCT species) FROM taxonomy;"` -ne "251" ]; then echo "Wrong species count"; false; fi
thapbi_pict ncbi-import -d $DB $TMP/ncbi_sample.fasta -s
if [ `sqlite3 $DB "SELECT COUNT(DISTINCT species) FROM taxonomy;"` -ne "251" ]; then echo "Wrong species count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -lt "100" ]; then echo "Taxonomy count too low"; false; fi
# Other values subject to change

thapbi_pict dump -d $DB -o /dev/null
thapbi_pict dump -d $DB -o /dev/null -g Phytophthora
thapbi_pict dump -d $DB -o /dev/null -g Phytophthora -s "ilicis, sp. aff. meadii"

echo "$0 passed"
