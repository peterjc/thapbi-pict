#!/bin/bash
IFS=$'\n\t'
set -eux
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "Checking ncbi-import"
thapbi_pict ncbi-import 2>&1 | grep "the following arguments are required"
set -o pipefail

# Cannot use validation without having some taxonomy entries
set +o pipefail
thapbi_pict ncbi-import -d sqlite:///:memory: tests/ncbi-import/20th_Century_ITS1.fasta 2>&1 | grep "Taxonomy table empty"
set -o pipefail

export DB=$TMP/20th_Century_ITS1.sqlite
rm -rf $DB
thapbi_pict ncbi-import -x -d $DB tests/ncbi-import/20th_Century_ITS1.fasta
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_source;"` -ne "107" ]; then echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "90" ]; then echo "Wrong its1_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "42" ]; then echo "Wrong taxonomy count"; false; fi
# Other values subject to change

thapbi_pict dump -d $DB -o /dev/null

if [ ! -f "new_taxdump_2018-12-01.zip" ]; then curl -L -O "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/new_taxdump_2018-12-01.zip"; fi
if [ ! -d "new_taxdump_2018-12-01" ]; then unzip new_taxdump_2018-12-01.zip -d new_taxdump_2018-12-01; fi

export DB=$TMP/20th_Century_ITS1_validated.sqlite
rm -rf $DB
thapbi_pict load-tax -d $DB -t new_taxdump_2018-12-01 -a 4783
if [ `sqlite3 $DB "SELECT COUNT(DISTINCT genus) FROM taxonomy;"` -ne "1" ]; then echo "Wrong genus count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(DISTINCT species) FROM taxonomy;"` -ne "252" ]; then echo "Wrong species count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "252" ]; then echo "Wrong taxonomy count"; false; fi
thapbi_pict ncbi-import -d $DB tests/ncbi-import/20th_Century_ITS1.fasta
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_source;"` -ne "104" ]; then echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "87" ]; then echo "Wrong its1_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "252" ]; then echo "Wrong taxonomy count"; false; fi
# Other values subject to change

thapbi_pict dump -d $DB -o /dev/null
thapbi_pict dump -d $DB -o /dev/null -g Phytophthora
thapbi_pict dump -d $DB -o /dev/null -g Phytophthora -s "ilicis, sp. aff. meadii"

# Now using the -g option,
export DB=$TMP/20th_Century_ITS1_genus_only.sqlite
rm -rf $DB
thapbi_pict load-tax -d $DB -t new_taxdump_2018-12-01
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "604" ]; then echo "Wrong taxonomy count"; false; fi
thapbi_pict ncbi-import -d $DB tests/ncbi-import/20th_Century_ITS1.fasta -g
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_source;"` -ne "106" ]; then echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "89" ]; then echo "Wrong its1_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "604" ]; then echo "Wrong taxonomy count"; false; fi
if [ `thapbi_pict dump -d $DB | grep -c 4783` -ne "106" ]; then echo "Should all be just genus"; false; fi

# Now using the Phytophthora at species level, Peronosporaceae at genus level:
export DB=$TMP/20th_Century_ITS1_mixed.sqlite
rm -rf $DB
thapbi_pict load-tax -d $DB -t new_taxdump_2018-12-01
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "604" ]; then echo "Wrong taxonomy count"; false; fi
thapbi_pict ncbi-import -d $DB tests/ncbi-import/20th_Century_ITS1.fasta
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_source;"` -ne "104" ]; then echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "87" ]; then echo "Wrong its1_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "604" ]; then echo "Wrong taxonomy count"; false; fi
if [ `thapbi_pict dump -d $DB -f fasta | grep "^>" | grep  "species=Phytophthora " -c` -ne 104 ]; then echo "Wrong Phytophthora species count"; false; fi
thapbi_pict ncbi-import -d $DB tests/ncbi-import/20th_Century_ITS1_Peronosporaceae.fasta -g
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "2" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_source;"` -ne "239" ]; then echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "95" ]; then echo "Wrong its1_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "604" ]; then echo "Wrong taxonomy count"; false; fi
# 104 + 106 + 29 = 239
if [ `thapbi_pict dump -d $DB -f fasta | grep "^>" | grep -c "species=Phytophthora "` -ne 104 ]; then echo "Wrong Phytophthora species count";	false; fi
if [ `thapbi_pict dump -d $DB -f fasta | grep "^>" | grep -c "species=Phytophthora\]"` -n 106 ]; then echo "Wrong Phytophthora genus-only count"; false; fi
if [ `thapbi_pict dump -d $DB -f fasta | grep "^>" | grep -c -v "Phytophthora"` -ne 29 ]; then echo "Wrong sister genus count"; false; fi

echo "$0 - test_ncbi-import.sh passed"
