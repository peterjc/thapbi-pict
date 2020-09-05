#!/bin/bash
set -eup pipeline

echo NOTE: Expected first time run time is about 1 hour,
echo repeat runs about 1 minute just to regenerate reports.
echo
echo =====================
echo Recycled water - ITS1
echo =====================

# Following the documentation would first run the pipeline with
# default settings for primer and database:
#
# mkdir intermediate_defaults/ summary_defaults/
# thapbi_pict pipeline \
#             -i raw_data/ -o summary_defaults/ -s intermediate_defaults/ \
#             -r recycled-water-defaults -t metadata.tsv -x 7 -c 1,2,3,4,5,6

echo "Building ITS1 database"
# Remove any pre-existing DB
rm -rf Redekar_et_al_2019_sup_table_3.sqlite

# Loading NCBI taxonomy for handling synonyms
thapbi_pict load-tax -d Redekar_et_al_2019_sup_table_3.sqlite -t taxdmp_2019-12-01/

# Using -x / --lax (does not insist on taxonomy match)
# Not giving primers, sequences are already trimmed
thapbi_pict curated-import -x \
        -d Redekar_et_al_2019_sup_table_3.sqlite \
        -i Redekar_et_al_2019_sup_table_3.fasta

echo "Drawing edit-graph for database entries alone"
# Using -s / --showdb
thapbi_pict edit-graph --showdb \
        -d Redekar_et_al_2019_sup_table_3.sqlite \
        -o Redekar_et_al_2019_sup_table_3.xgmml

echo "Running analysis"
thapbi_pict pipeline -i raw_data/ -s intermediate/ -o summary/ \
        --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA \
        --right AGCGTTCTTCATCGATGTGC \
        -d Redekar_et_al_2019_sup_table_3.sqlite --showdb -m onebp \
        -r recycled-water-custom -t metadata.tsv -x 7 -c 1,2,3,4,5,6
