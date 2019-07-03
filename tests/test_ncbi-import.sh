#!/bin/bash

# Copyright 2018-2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eux
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp}

echo "Checking ncbi-import"
thapbi_pict ncbi-import 2>&1 | grep "the following arguments are required"
# Cannot use validation without having some taxonomy entries
thapbi_pict ncbi-import -d sqlite:///:memory: --input tests/ncbi-import/20th_Century_ITS1.fasta 2>&1 | grep "Taxonomy table empty"
set -o pipefail

# examples with multiple HMM matches in the sequence (tandem repeats etc)
export DB=$TMP/multiple_hmm.sqlite
rm -rf $DB
thapbi_pict load-tax -d $DB -t new_taxdump_2018-12-01
# NCBI import at genus level only, as used in bundled ITS1_DB.sqlite
thapbi_pict ncbi-import -d $DB -g -i tests/ncbi-import/multiple_hmm.fasta -n "NCBI examples with multiple HMM matches"
# WARNING: 2 HMM matches in MF370571.1
# WARNING: 2 HMM matches in MH169111.1
# WARNING: Discarding exactly duplicated ITS1 matches in DQ641247.1
# File tests/ncbi-import/multiple_hmm.fasta had 5 sequences. Found 5 ITS1, of which 5 accepted.
# Of 5 potential entries, 0 unparsable, 2 failed sp. validation, 3 OK.
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_source;"` -ne "5" ]; then echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "5" ]; then echo "Wrong its1_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "604" ]; then echo "Wrong taxonomy count"; false; fi
# Debugging output,
# $ sqlite3 $DB "SELECT md5, LENGTH(its1_sequence.sequence), source_accession FROM its1_sequence, its1_source WHERE its1_sequence.id=its1_source.its1_id;"
# 63fa728c0fe76536f13eb593df99bd46|179|MF370571.1
# a42d385f25a9f1b10d0642ad8a72a584|116|MF370571.1
# 4c9e98f437ca0f55d0d8ba3b2928239c|199|MH169111.1
# 158774f117b3e6058b22ba9ef877f346|199|MH169111.1
# 7f27d3a8f7150e0ee7ad64073e6da6b5|170|DQ641247.1
if [ `sqlite3 $DB "SELECT MAX(LENGTH(sequence)) FROM its1_sequence;"` -ne "199" ]; then echo "Wrong max ITS1 sequence length"; false; fi


export DB=$TMP/20th_Century_ITS1.sqlite
rm -rf $DB
thapbi_pict ncbi-import -x -d $DB -i tests/ncbi-import/20th_Century_ITS1.fasta
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
thapbi_pict ncbi-import -d $DB -i tests/ncbi-import/20th_Century_ITS1.fasta
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
thapbi_pict ncbi-import -d $DB -i tests/ncbi-import/20th_Century_ITS1.fasta -g
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
thapbi_pict ncbi-import -d $DB -i tests/ncbi-import/20th_Century_ITS1.fasta
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_source;"` -ne "104" ]; then echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "87" ]; then echo "Wrong its1_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "604" ]; then echo "Wrong taxonomy count"; false; fi
if [ `thapbi_pict dump -d $DB -f fasta | grep "^>" | grep  "species=Phytophthora " -c` -ne 104 ]; then echo "Wrong Phytophthora species count"; false; fi
thapbi_pict ncbi-import -d $DB -i tests/ncbi-import/20th_Century_ITS1_Peronosporaceae.fasta -g
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "2" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_source;"` -ne "239" ]; then echo "Wrong its1_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM its1_sequence;"` -ne "95" ]; then echo "Wrong its1_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "604" ]; then echo "Wrong taxonomy count"; false; fi
# 104 + 106 + 29 = 239
if [ `thapbi_pict dump -d $DB -f fasta | grep "^>" | grep -c "species=Phytophthora "` -ne 104 ]; then echo "Wrong Phytophthora species count";	false; fi
if [ `thapbi_pict dump -d $DB -f fasta | grep "^>" | grep -c "species=Phytophthora\]"` -ne 106 ]; then echo "Wrong Phytophthora genus-only count"; false; fi
if [ `thapbi_pict dump -d $DB -f fasta | grep "^>" | grep -c -v "Phytophthora"` -ne 29 ]; then echo "Wrong sister genus count"; false; fi

echo "$0 - test_ncbi-import.sh passed"
