#!/bin/bash
VERSION=`thapbi_pict -v | sed "s/THAPBI PICT //g"`
echo "Using THAPBI PICT $VERSION"
set -euo pipefail
TAX=taxdmp_2021-11-01
DB=ITS1_DB
rm -rf "$DB.sqlite" "$DB.fasta" "$DB.txt" "$DB.sql"

# 4762 = Oomycetes
thapbi_pict load-tax -d "$DB.sqlite" -t "$TAX" -a 4762

# Ad-hoc fix for NCBI taxonomy not yet having caught up with community consensus.
# At the 7th Meeting of the International Union of Forest Research Organisations
# Working Party (IUFRO) 7.02.09, Phytophthoras in forests and natural ecosystems,
# a decision was made to adhere to the original and correct version of the species
# name, Phytophthora austrocedri.
#
# First, rename 'Phytophthora austrocedri' synonym to 'Phytophthora austrocedrae'
sqlite3 "$DB.sqlite" "UPDATE synonym SET name='Phytophthora austrocedrae' WHERE name='Phytophthora austrocedri';"
# Change 'Phytophthora austrocedrae' main entry to 'Phytophthora austrocedri'
sqlite3 "$DB.sqlite" "UPDATE taxonomy SET species='austrocedri' WHERE genus='Phytophthora' AND species='austrocedrae'"
# Should now be able to import data using either name.

# ==========================
# G-BLOCK synthetic controls
# ==========================
sqlite3 "$DB.sqlite" "INSERT INTO taxonomy (ncbi_taxid, genus, species) VALUES (32630, 'synthetic', 'construct C1');"
sqlite3 "$DB.sqlite" "INSERT INTO taxonomy (ncbi_taxid, genus, species) VALUES (32630, 'synthetic', 'construct C2');"
sqlite3 "$DB.sqlite" "INSERT INTO taxonomy (ncbi_taxid, genus, species) VALUES (32630, 'synthetic', 'construct C3');"
sqlite3 "$DB.sqlite" "INSERT INTO taxonomy (ncbi_taxid, genus, species) VALUES (32630, 'synthetic', 'construct C4');"
# This also defines the marker name, primers, and length range for preparing reads
thapbi_pict import -d "$DB.sqlite" -i controls.fasta \
            -k ITS1 -l GAAGGTGAAGTCGTAACAAGG -r GCARRGACTTTCGTCCCYRC \
            --minlen 100 --maxlen 1000

# ===============
# Single isolates
# ===============
# FASTA files prepared via thapbi_pict prepare-reads and curated-seq steps:
thapbi_pict import -d "$DB.sqlite" -i single_isolates/*.fasta

# ===================
# NCBI at genus level
# ===================
thapbi_pict import -d "$DB.sqlite" -i 2021-11-22-ITS_Oomycota_w32.fasta -c ncbi -g --minlen 150 --maxlen 450
# Add hoc fix for some accessions apparently with wrong genus
BAD="('MN128447.1', 'MK794853.1', 'MK794848.1', 'MK795051.1', 'HQ237483.1', 'KP183959.1', 'MW426376.1', 'MW426384.1', 'KY785380.1', 'KY785381.1', 'KU715054.1', 'GQ149496.1', 'JF916542.1')"
sqlite3 "$DB.sqlite" "DELETE FROM marker_sequence WHERE id IN (SELECT marker_sequence.id FROM marker_sequence JOIN sequence_source ON marker_sequence.id = sequence_source.marker_seq_id WHERE source_accession IN $BAD);"
sqlite3 "$DB.sqlite" "DELETE FROM sequence_source WHERE source_accession IN $BAD;"

# =================
# Curated sequences
# =================
thapbi_pict import -d "$DB.sqlite" -i Phytophthora_ITS1_curated.fasta Nothophytophthora_ITS1_curated.fasta --maxlen 450 -s ";"

thapbi_pict dump -m -d "$DB.sqlite" -o "$DB.txt"
thapbi_pict dump -m -f fasta -d "$DB.sqlite" -o "$DB.fasta"

sqlite3 "$DB.sqlite" .dump > "$DB.sql"

cp "$DB.sqlite" "$DB-$VERSION.sqlite"
cp "$DB.sql" "$DB-$VERSION.sql"
cp "$DB.txt" "$DB-$VERSION.txt"
cp "$DB.fasta" "$DB-$VERSION.fasta"

thapbi_pict conflicts -d "$DB.sqlite"
