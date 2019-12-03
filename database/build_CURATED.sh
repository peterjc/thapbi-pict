#!/bin/bash
VERSION=`thapbi_pict -v | sed "s/THAPBI PICT //g"`
echo "Using THAPBI PICT $VERSION"
set -euo pipefail
CURATED=Phytophthora_ITS1_curated.fasta
TAX=taxdmp_2019-12-01
DB=CURATED
rm -rf "$DB.sqlite" "$DB.fasta" "$DB.txt" "$DB.sql"

thapbi_pict load-tax -d "$DB.sqlite" -t "$TAX"
thapbi_pict curated-import -d "$DB.sqlite" -i "$CURATED" -v

# Ad-hoc fix for NCBI taxonomy not yet having caught up with community consensus.
# At the 7th Meeting of the International Union of Forest Research Organisations
# Working Party (IUFRO) 7.02.09, Phytophthoras in forests and natural ecosystems,
# a decision was made to adhere to the original and correct version of the species
# name, Phytophthora austrocedri.
sqlite3 "$DB.sqlite" "UPDATE taxonomy SET species='austrocedri' WHERE genus='Phytophthora' AND species='austrocedrae'"

# Add the G-BLOCK synthetic controls
sqlite3 "$DB.sqlite" "INSERT INTO taxonomy (ncbi_taxid, genus, species) VALUES (32630, 'synthetic', 'construct C1');"
sqlite3 "$DB.sqlite" "INSERT INTO taxonomy (ncbi_taxid, genus, species) VALUES (32630, 'synthetic', 'construct C2');"
sqlite3 "$DB.sqlite" "INSERT INTO taxonomy (ncbi_taxid, genus, species) VALUES (32630, 'synthetic', 'construct C3');"
sqlite3 "$DB.sqlite" "INSERT INTO taxonomy (ncbi_taxid, genus, species) VALUES (32630, 'synthetic', 'construct C4');"
thapbi_pict curated-import -d "$DB.sqlite" -i controls.fasta

thapbi_pict dump -m -d "$DB.sqlite" -o "$DB.txt"
thapbi_pict dump -m -f fasta -d "$DB.sqlite" -o "$DB.fasta"

sqlite3 "$DB.sqlite" .dump > "$DB.sql"

cp "$DB.sqlite" "$DB-$VERSION.sqlite"
cp "$DB.sql" "$DB-$VERSION.sql"
cp "$DB.txt" "$DB-$VERSION.txt"
cp "$DB.fasta" "$DB-$VERSION.fasta"

thapbi_pict conflicts -d "$DB.sqlite"
