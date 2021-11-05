#!/bin/bash

# Copyright 2018-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

IFS=$'\n\t'
set -eu
# Note not using "set -o pipefail" until after check error message with grep

export TMP=${TMP:-/tmp/thapbi_pict}/ncbi_import
rm -rf $TMP
mkdir -p $TMP

export LEFT=GAAGGTGAAGTCGTAACAAGG
export RIGHT=GCARRGACTTTCGTCCCYRC
export RIGHT_RC=GYRGGGACGAAAGTCYYTGC

echo "===================="
echo "Checking ncbi-import"
echo "===================="
set -x
thapbi_pict import 2>&1 | grep "the following arguments are required"
# Cannot use validation without having some taxonomy entries
thapbi_pict import -c ncbi -d sqlite:///:memory: --input tests/ncbi-import/20th_Century_ITS1.fasta 2>&1 | grep "Taxonomy table empty"
set -o pipefail

if [ ! -f "new_taxdump_2019-09-01.zip" ]; then curl -L -O "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/new_taxdump_2019-09-01.zip"; fi
if [ ! -d "new_taxdump_2019-09-01" ]; then unzip new_taxdump_2019-09-01.zip nodes.dmp names.dmp -d new_taxdump_2019-09-01; fi


# Check hybrid like "Phytophthora humicola x Phytophthora inundata"
# imports as genus="Phytophthora", species="humicola x inundata"
export DB=$TMP/hybrid.sqlite
rm -rf $DB
thapbi_pict load-tax -d $DB -t new_taxdump_2019-09-01
thapbi_pict import -k ITS1 -l $LEFT -r $RIGHT -c ncbi -d $DB -i tests/ncbi-import/hybrid.fasta
if [ `thapbi_pict dump -d $DB | grep -c "humicola x inundata"` -ne "1" ]; then echo "Expected humicola x inundata"; false; fi


# examples with multiple HMM matches in the sequence (tandem repeats etc)
export DB=$TMP/multiple_hmm.sqlite
rm -rf $DB
thapbi_pict load-tax -d $DB -t new_taxdump_2019-09-01
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "3003" ]; then echo "Wrong taxonomy count"; false; fi
# NCBI import at genus level only, as used in bundled ITS1_DB.sqlite
# 5 sequences all with multiple HMM matches, but after primer trimming most become single HMM entries.
cutadapt --quiet -g $LEFT tests/ncbi-import/multiple_hmm.fasta \
| cutadapt --quiet -a $RIGHT_RC -o $TMP/multiple_hmm.fasta /dev/stdin
thapbi_pict import -k ITS1 -l $LEFT -r $RIGHT -c ncbi -d $DB -g -i $TMP/multiple_hmm.fasta \
            -n "NCBI examples with multiple HMM matches"
# WARNING: 2 HMM matches in MF095142.1
# WARNING: Uncultured, so ignoring 'MF095142.1 Uncultured Peronosporaceae clone MZOo17 small subunit ri...'
# WARNING: 2 HMM matches in KP691407.1
# WARNING: Uncultured, so ignoring 'KP691407.1 Uncultured Phytophthora clone sp3 18S ribosomal RNA gene...'
# File tests/ncbi-import/multiple_hmm.fasta had 5 sequences. Found 5 with ITS1, of which 3 accepted.
# Of 5 potential entries, 0 unparsable, 2 failed sp. validation, 3 OK.
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM sequence_source;"` -ne "3" ]; then echo "Wrong sequence_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_sequence;"` -ne "3" ]; then echo "Wrong marker_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "3003" ]; then echo "Wrong taxonomy count"; false; fi
# Debugging output,
# $ sqlite3 $DB "SELECT md5, LENGTH(marker_sequence.sequence), source_accession FROM marker_sequence, sequence_source WHERE marker_sequence.id=sequence_source.its1_id;"
# 63fa728c0fe76536f13eb593df99bd46|179|MF370571.1
# 4c9e98f437ca0f55d0d8ba3b2928239c|199|MH169111.1
# 7f27d3a8f7150e0ee7ad64073e6da6b5|170|DQ641247.1
if [ `sqlite3 $DB "SELECT MAX(LENGTH(sequence)) FROM marker_sequence;"` -ne "217" ]; then echo "Wrong max ITS1 sequence length"; false; fi

# When importing NCBI files, we no longer assume P. is Phytophthora:
rm -rf $TMP/20th_Century_ITS1.fasta $TMP/20th_Century_ITS1_Peronosporaceae.fasta
cat tests/ncbi-import/20th_Century_ITS1.fasta | sed "s/ P\./ Phytophthora /g" \
| cutadapt --quiet -g $LEFT /dev/stdin \
| cutadapt --quiet -a $RIGHT_RC -o $TMP/20th_Century_ITS1.fasta /dev/stdin
cat tests/ncbi-import/20th_Century_ITS1_Peronosporaceae.fasta | sed "s/ P\./ Phytophthora /g" \
| cutadapt --quiet -g $LEFT /dev/stdin \
| cutadapt --quiet -a $RIGHT_RC -o $TMP/20th_Century_ITS1_Peronosporaceae.fasta /dev/stdin

export DB=$TMP/20th_Century_ITS1.sqlite
rm -rf $DB
thapbi_pict import -k ITS1 -l $LEFT -r $RIGHT -c ncbi -x -d $DB -i $TMP/20th_Century_ITS1.fasta
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM sequence_source;"` -ne "120" ]; then echo "Wrong sequence_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_sequence;"` -ne "103" ]; then echo "Wrong marker_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "42" ]; then echo "Wrong taxonomy count"; false; fi
# Other values subject to change

thapbi_pict dump -d $DB -o /dev/null

export DB=$TMP/20th_Century_ITS1_validated.sqlite
rm -rf $DB
thapbi_pict load-tax -d $DB -t new_taxdump_2019-09-01 -a 4783
if [ `sqlite3 $DB "SELECT COUNT(DISTINCT genus) FROM taxonomy;"` -ne "1" ]; then echo "Wrong genus count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(DISTINCT species) FROM taxonomy;"` -ne "832" ]; then echo "Wrong species count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "832" ]; then echo "Wrong taxonomy count"; false; fi
thapbi_pict import -k ITS1 -l $LEFT -r $RIGHT -c ncbi -d $DB -i $TMP/20th_Century_ITS1.fasta
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM sequence_source;"` -ne "120" ]; then echo "Wrong sequence_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_sequence;"` -ne "103" ]; then echo "Wrong marker_sequence count"; false; fi
# Confirm AF271230.1 Pythium undulatum -> Phytophthora undulatum
thapbi_pict dump -d $DB -k ITS1 | cut -f 1-4 | grep AF271230.1 | grep Phytophthora
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "832" ]; then echo "Wrong taxonomy count"; false; fi
# Other values subject to change

thapbi_pict dump -d $DB -o /dev/null
thapbi_pict dump -d $DB -o /dev/null -g Phytophthora
thapbi_pict dump -d $DB -o /dev/null -g Phytophthora -s "ilicis, sp. aff. meadii"

# Now using the -g option,
export DB=$TMP/20th_Century_ITS1_genus_only.sqlite
rm -rf $DB
thapbi_pict load-tax -d $DB -t new_taxdump_2019-09-01
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "3003" ]; then echo "Wrong taxonomy count"; false; fi
thapbi_pict import -k ITS1 -l $LEFT -r $RIGHT -c ncbi -d $DB -i $TMP/20th_Century_ITS1.fasta -g
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM sequence_source;"` -ne "120" ]; then echo "Wrong sequence_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_sequence;"` -ne "103" ]; then echo "Wrong marker_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "3003" ]; then echo "Wrong taxonomy count"; false; fi
# Confirm AF271230.1 Pythium undulatum -> Phytophthora undulatum -> Phytophthora
thapbi_pict dump -d $DB | cut -f1-3 | grep AF271230.1 | grep Phytophthora
if [ `thapbi_pict dump -d $DB | grep -c 4783` -ne "120" ]; then echo "Should all be just genus"; false; fi

# Now using the Phytophthora at species level, Peronosporaceae at genus level:
export DB=$TMP/20th_Century_ITS1_mixed.sqlite
rm -rf $DB
thapbi_pict load-tax -d $DB -t new_taxdump_2019-09-01
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "3003" ]; then echo "Wrong taxonomy count"; false; fi
thapbi_pict import -k ITS1 -l $LEFT -r $RIGHT -c ncbi -d $DB -i $TMP/20th_Century_ITS1.fasta
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "1" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM sequence_source;"` -ne "120" ]; then echo "Wrong sequence_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_sequence;"` -ne "103" ]; then echo "Wrong marker_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "3003" ]; then echo "Wrong taxonomy count"; false; fi
if [ `thapbi_pict dump -d $DB -f fasta | grep "^>" | grep  " Phytophthora " -c` -ne 103 ]; then echo "Wrong Phytophthora species count"; false; fi
# Adding to existing DB, do not need to give primers again:
thapbi_pict import -c ncbi -k ITS1 -d $DB -i $TMP/20th_Century_ITS1_Peronosporaceae.fasta -g
if [ `sqlite3 $DB "SELECT COUNT(id) FROM data_source;"` -ne "2" ]; then echo "Wrong data_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM sequence_source;"` -ne "269" ]; then echo "Wrong sequence_source count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM marker_sequence;"` -ne "109" ]; then echo "Wrong marker_sequence count"; false; fi
if [ `sqlite3 $DB "SELECT COUNT(id) FROM taxonomy;"` -ne "3003" ]; then echo "Wrong taxonomy count"; false; fi
# 120 + 120 + 29 = 269
if [ `thapbi_pict dump -d $DB | grep -v "^#" | grep -c Phytophthora$'\t'[a-z]` -ne 120 ]; then echo "Wrong Phytophthora species count"; false; fi
if [ `thapbi_pict dump -d $DB | grep -v "^#" | grep -c Phytophthora$'\t\t'` -ne 120 ]; then echo "Wrong Phytophthora genus-only count"; false; fi
if [ `thapbi_pict dump -d $DB | grep -v "^#" | grep -c -v Phytophthora` -ne 29 ]; then echo "Wrong sister genus count"; false; fi

echo "$0 - test_ncbi-import.sh passed"
