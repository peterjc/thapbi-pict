#!/bin/bash
IFS=$'\n\t'
set -eux
# Note not using "set -o pipefile" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "Checking ncbi-import"
thapbi_pict ncbi-import 2>&1 | grep "the following arguments are required"
thapbi_pict dump 2>&1 | grep "the following arguments are required"
set -o pipefail

# Couldn't see how to set a limit on the number of records via
# esearch/efetch command line. Instead to avoid timeouts setting
# another limit on the publication date - using last century so
# the search results should now be stable.
if [ ! -f $TMP/20th_Century_ITS1.fasta ]; then esearch -db nucleotide -query "its1 AND Phytophthora[Organism] AND 150:800[Sequence Length] AND 1900:2000[Publication Date]" | efetch -format fasta > $TMP/20th_Century_ITS1.fasta; fi

if [ `grep -c "^>" $TMP/20th_Century_ITS1.fasta` -ne 129 ]; then echo "Record count from NCBI Entrez changed"; false; fi

# Cannot use validation without having some taxonomy entries
set +o pipefail
thapbi_pict ncbi-import -d sqlite:///:memory: $TMP/20th_Century_ITS1.fasta 2>&1 | grep "Taxonomy table empty"
set -o pipefail

export DB=$TMP/20th_Century_ITS1.sqlite
rm -rf $DB
thapbi_pict ncbi-import -x -d $DB $TMP/20th_Century_ITS1.fasta

if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "96" ]; then echo "Wrong its1_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "43" ]; then echo "Wrong taxonomy count"; false; fi
# Other values subject to change

thapbi_pict dump -d $DB -o /dev/null

if [ ! -f "new_taxdump_2018-12-01.zip" ]; then curl -L -O "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/new_taxdump_2018-12-01.zip"; fi
if [ ! -d "new_taxdump_2018-12-01" ]; then unzip new_taxdump_2018-12-01.zip -d new_taxdump_2018-12-01; fi

export DB=$TMP/20th_Century_ITS1_validated.sqlite
rm -rf $DB
thapbi_pict load-tax -d $DB -t new_taxdump_2018-12-01 -a 4783
if [ `sqlite3 $DB "SELECT COUNT(DISTINCT genus) FROM taxonomy;"` -ne "1" ]; then echo "Wrong genus count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(DISTINCT species) FROM taxonomy;"` -ne "251" ]; then echo "Wrong species count"; false; fi
thapbi_pict ncbi-import -d $DB $TMP/20th_Century_ITS1.fasta
if [ `sqlite3 $DB "SELECT COUNT(DISTINCT species) FROM taxonomy;"` -ne "251" ]; then echo "Wrong species count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "96" ]; then echo "Wrong its1_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "251" ]; then echo "Wrong taxonomy count"; false; fi
# Other values subject to change

thapbi_pict dump -d $DB -o /dev/null
thapbi_pict dump -d $DB -o /dev/null -g Phytophthora
thapbi_pict dump -d $DB -o /dev/null -g Phytophthora -s "ilicis, sp. aff. meadii"

echo "$0 passed"
