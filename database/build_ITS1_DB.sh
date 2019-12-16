#!/bin/bash
VERSION=`thapbi_pict -v | sed "s/THAPBI PICT //g"`
echo "Using THAPBI PICT $VERSION"
set -euo pipefail
CURATED=Phytophthora_ITS1_curated.fasta
TAX=taxdmp_2019-12-01
DB=ITS1_DB
rm -rf "$DB.sqlite" "$DB.fasta" "$DB.txt" "$DB.sql"

# 4762 = Oomycetes
thapbi_pict load-tax -d "$DB.sqlite" -t "$TAX" -a 4762
thapbi_pict curated-import -d "$DB.sqlite" -i "$CURATED"
thapbi_pict ncbi-import -d "$DB.sqlite" -i 2019-04-03-ITS_Peronosporales_16394.fasta -g --minlen 150 --maxlen 750

# Add hoc fix for two accessions potentially having wrong genus
sqlite3 "$DB.sqlite" "DELETE FROM its1_sequence WHERE id IN (SELECT its1_sequence.id FROM its1_sequence JOIN its1_source ON its1_sequence.id = its1_source.its1_id WHERE source_accession IN ('AY742739.1', 'JX122744.1'));"
sqlite3 "$DB.sqlite" "DELETE FROM its1_source WHERE source_accession IN ('AY742739.1', 'JX122744.1');"

# Ad-hoc fix for NCBI taxonomy not yet having caught up with community consensus.
# At the 7th Meeting of the International Union of Forest Research Organisations
# Working Party (IUFRO) 7.02.09, Phytophthoras in forests and natural ecosystems,
# a decision was made to adhere to the original and correct version of the species
# name, Phytophthora austrocedri.
#
# First, rename 'Phytophthora austrocedri' synonym to 'Phytophthora austrocedrae'
sqlite3 ITS1_DB.sqlite "UPDATE synonym SET name='Phytophthora austrocedrae' WHERE name='Phytophthora austrocedri';"
# Change 'Phytophthora austrocedrae' main entry to 'Phytophthora austrocedri'
sqlite3 "$DB.sqlite" "UPDATE taxonomy SET species='austrocedri' WHERE genus='Phytophthora' AND species='austrocedrae'"
# Should now be able to import data using either name.
# Will use 'Phytophthora austrocedri' in clasffier/assess

# The known value files are now using Phytophthora austrocedri, not P. austrocedrae
thapbi_pict seq-import -d "$DB.sqlite" -i thapbi20180709p1_MetaControls/prepared_reads_${VERSION}/*.fasta thapbi20180709p1_MetaControls/positive_controls/*.known.tsv

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
