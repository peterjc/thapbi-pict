#!/bin/bash
VERSION=`thapbi_pict -v | sed "s/THAPBI PICT //g"`
echo "Using THAPBI PICT $VERSION"
set -euo pipefail
TAX=new_taxdump_2019-01-01
DB=ITS1_DB
rm -rf "$DB.sqlite" "$DB.fasta" "$DB.txt" "$DB.sql"

thapbi_pict load-tax -d "$DB.sqlite" -t "$TAX"
thapbi_pict legacy-import -d "$DB.sqlite" legacy/Phytophthora_ITS_database_v0.005.fasta
thapbi_pict ncbi-import -d "$DB.sqlite" 2019-04-03-ITS_Peronosporales_16394.fasta -g
thapbi_pict seq-import -d "$DB.sqlite" thapbi20180709p1_MetaControls/prepared_reads_${VERSION}/*.fasta thapbi20180709p1_MetaControls/positive_controls/*.known.tsv

# Add the G-BLOCK synthetic controls (in lax mode as not in the taxonomy)
grep -A 1 ">Control_" legacy/Phytophthora_ITS_database_v0.005.fasta > controls.fasta
thapbi_pict legacy-import -x -d "$DB.sqlite" controls.fasta

# Ad-hoc fix for NCBI taxonomy not yet having caught up with community consensus.
# At the 7th Meeting of the International Union of Forest Research Organisations
# Working Party (IUFRO) 7.02.09, Phytophthoras in forests and natural ecosystems,
# a decision was made to adhere to the original and correct version of the species
# name, Phytophthora austrocedri.
sqlite3 ITS1_DB.sqlite "UPDATE taxonomy SET species='austrocedri' WHERE genus='Phytophthora' AND species='austrocedrae'"

thapbi_pict dump -d "$DB.sqlite" -o "$DB.txt"
thapbi_pict dump -f fasta -d "$DB.sqlite" -o "$DB.fasta"

sqlite3 "$DB.sqlite" .dump > "$DB.sql"
