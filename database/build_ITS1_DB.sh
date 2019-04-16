#!/bin/bash
VERSION=`thapbi_pict -v | sed "s/THAPBI PICT //g"`
echo "Using THAPBI PICT $VERSION"
set -euo pipefail
TAX=new_taxdump_2019-01-01
DB=ITS1_DB
rm -rf "$DB.sqlite" "$DB.fasta" "$DB.txt" "$DB.sql"

thapbi_pict load-tax -d "$DB.sqlite" -t "$TAX"
thapbi_pict legacy-import -d "$DB.sqlite" legacy/Phytophthora_ITS_database_v0.005.fasta
thapbi_pict ncbi-import -d "$DB.sqlite" 2019-04-03-ITS_Peronosporaceae_8336.fasta -g
thapbi_pict seq-import -d "$DB.sqlite" thapbi20180709p1_MetaControls/prepared_reads_${VERSION}/*.fasta thapbi20180709p1_MetaControls/positive_controls/*.known.tsv

thapbi_pict dump -d "$DB.sqlite" -o "$DB.txt"
thapbi_pict dump -f fasta -d "$DB.sqlite" -o "$DB.fasta"

sqlite3 "$DB.sqlite" .dump > "$DB.sql"
