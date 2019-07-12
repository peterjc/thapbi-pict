#!/bin/bash
VERSION=`thapbi_pict -v | sed "s/THAPBI PICT //g"`
echo "Using THAPBI PICT $VERSION"
set -euo pipefail
TAX=new_taxdump_2019-01-01
DB=ITS1_DB
rm -rf "$DB.sqlite" "$DB.fasta" "$DB.txt" "$DB.sql"

thapbi_pict load-tax -d "$DB.sqlite" -t "$TAX"
# In strict mode this will ignore the synthetic controls, we add them later:
thapbi_pict legacy-import -d "$DB.sqlite" -i legacy/Phytophthora_ITS_database_v0.005.fasta
thapbi_pict ncbi-import -d "$DB.sqlite" -i 2019-04-03-ITS_Peronosporales_16394.fasta -g

# Ad-hoc fix for NCBI taxonomy not yet having caught up with community consensus.
# At the 7th Meeting of the International Union of Forest Research Organisations
# Working Party (IUFRO) 7.02.09, Phytophthoras in forests and natural ecosystems,
# a decision was made to adhere to the original and correct version of the species
# name, Phytophthora austrocedri.
sqlite3 ITS1_DB.sqlite "UPDATE taxonomy SET species='austrocedri' WHERE genus='Phytophthora' AND species='austrocedrae'"

# The known value files are now using Phytophthora austrocedri, not P. austrocedrae
thapbi_pict seq-import -d "$DB.sqlite" -i thapbi20180709p1_MetaControls/prepared_reads_${VERSION}/*.fasta thapbi20180709p1_MetaControls/positive_controls/*.known.tsv

# Add the G-BLOCK synthetic controls (in lax mode as not in the taxonomy)
grep -A 1 ">Control_" legacy/Phytophthora_ITS_database_v0.005.fasta > controls.fasta
thapbi_pict legacy-import -x -d "$DB.sqlite" -i controls.fasta

# Ad-hoc fix for three unique sequences getting more than one genus in the NCBI,
# should be Hyaloperonospora not Peronospora. Drop the entries saying Peronospora:
sqlite3 ITS1_DB.sqlite "DELETE FROM its1_source WHERE its1_source.id in (SELECT its1_source.id FROM its1_source JOIN its1_sequence JOIN taxonomy WHERE its1_source.current_taxonomy_id=taxonomy.id AND its1_source.its1_id = its1_sequence.id AND taxonomy.genus='Peronospora' AND its1_sequence.md5 in ('71d4e062275062a6ed3863c71f137e77', 'e29a4c3d458c94842a8dc420dcfe946e', '320df1a347406a2eb13fbe329264ceb1'));"

thapbi_pict dump -d "$DB.sqlite" -o "$DB.txt"
thapbi_pict dump -f fasta -d "$DB.sqlite" -o "$DB.fasta"

sqlite3 "$DB.sqlite" .dump > "$DB.sql"
