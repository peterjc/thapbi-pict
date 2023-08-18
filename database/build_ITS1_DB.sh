#!/bin/bash
VERSION=`thapbi_pict -v | sed "s/THAPBI PICT //g"`
echo "Using THAPBI PICT $VERSION"
set -euo pipefail
TAX=taxdmp_2023-08-01
DB=ITS1_DB
rm -rf "$DB.sqlite" "$DB.fasta" "$DB.txt" "$DB.sql"

# ========
# Taxonomy
# ========
if [ ! -f "${TAX}.zip" ]; then curl -L -O "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_archive/${TAX}.zip"; fi
if [ ! -d "${TAX}/" ]; then unzip ${TAX}.zip nodes.dmp names.dmp merged.dmp -d ${TAX}/; fi

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
# These are multi-entries per FASTA sequence, using semi-colon:
thapbi_pict import -d "$DB.sqlite" -s ";" -g -s ";" \
            -i 2023-08-18_ITS1_Oomycota_obs.fasta

# Single entry per FASTA sequence, where semi-colon may appear in free text.
# Length limits deliberately more cautious than the pipeline settings above.
thapbi_pict import -d "$DB.sqlite" -c ncbi -g --minlen 150 --maxlen 450 \
            -i 2023-08-18_ITS1_Oomycota_w32.fasta

# Add hoc fixes for some accessions apparently with wrong genus
# MN128447.1 probably not Phytopythium vexans, but Phytophthora
sqlite3 "$DB.sqlite" "DELETE FROM sequence_source WHERE source_accession IN ('MN128447.1');"
# GQ149496/JF916542 probably not Phytophthora, but budding yeast
sqlite3 "$DB.sqlite" "DELETE FROM sequence_source WHERE source_accession IN ('GQ149496.1', 'JF916542.1');"
# KY986679 probably not Peronospora erucastri, but Hyaloperonospora
sqlite3 "$DB.sqlite" "DELETE FROM sequence_source WHERE source_accession IN ('KY986679.X');"
# Assorted entries, mostly Pythium which probably should be Globisporangium or Phytopythium etc
BAD="('MK794853.1', 'MK794848.1', 'MK795051.1', 'KY785380.1', 'KY785381.1', 'KU715054.1', 'MK794725.1', 'MK794726.1', 'ON394668.1', 'ON394667.1', 'ON394669.1', 'KU210557.1', 'ON075205.1', 'MZ799354.1', 'MZ799355.1', 'LC684551.1', 'ON394664.1', 'ON394670.1', 'KY785379.1', 'KU715059.1', 'MK794740.1', 'MK794746.1', 'ON394673.1', 'AY684925.1')"
sqlite3 "$DB.sqlite" "DELETE FROM sequence_source WHERE source_accession IN $BAD;"
# Recent entry OR398863.1 Globisporangium vs five older accessions saying Pythium
sqlite3 "$DB.sqlite" "DELETE FROM sequence_source WHERE source_accession IN ('OR398863.1');"

# =================
# Curated sequences
# =================

# Using lax mode (-x) since we have some entries not
# in the NCBI taxonomy, like putative novel species.
thapbi_pict import -d "$DB.sqlite" -i Phytophthora_ITS1_curated.fasta Nothophytophthora_ITS1_curated.fasta --maxlen 450 -s ";" -x

# ======
# Export
# ======

thapbi_pict dump -d "$DB.sqlite" -o "$DB.tsv"
thapbi_pict dump -m -d "$DB.sqlite" -o "$DB.txt"
thapbi_pict dump -f fasta -d "$DB.sqlite" -o "$DB.fasta"

sqlite3 "$DB.sqlite" .dump > "$DB.sql"

cp "$DB.sqlite" "$DB-$VERSION.sqlite"
cp "$DB.sql" "$DB-$VERSION.sql"
cp "$DB.tsv" "$DB-$VERSION.tsv"
cp "$DB.txt" "$DB-$VERSION.txt"
cp "$DB.fasta" "$DB-$VERSION.fasta"

thapbi_pict conflicts -d "$DB.sqlite"
