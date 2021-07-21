#!/bin/bash
set -euo pipefail

echo "NOTE: Expected first time run time is under an hour,"
echo "repeat runs about 1 minute (mostly on the edit-graph)."
echo

mkdir -p intermediate_defaults/ intermediate_long/ summary/ tmp_merged/

echo =========================
echo Recycled water - Defaults
echo =========================

echo "First with default settings and DB"
thapbi_pict pipeline \
        -i raw_data/ -o summary/recycled-water-defaults \
        -s intermediate_defaults/ \
        --merged-cache tmp_merged/ \
        -t metadata.tsv -x 7 -c 1,2,3,4,5,6

echo ==========================
echo Recycled water - Custom DB
echo ==========================

if [ ! -f Redekar_et_al_2019_sup_table_3.sqlite ]; then
    echo "Building ITS1 database"

    # Loading NCBI taxonomy for handling synonyms
    thapbi_pict load-tax -d Redekar_et_al_2019_sup_table_3.sqlite -t taxdmp_2019-12-01/

    # Using -x / --lax (does not insist on taxonomy match)
    # Adding 32bp conserved TTTCCGTAGGTGAACCTGCGGAAGGATCATTA to left primer
    # Not giving primers, sequences are already trimmed
    thapbi_pict import -x -s ";" \
                -d Redekar_et_al_2019_sup_table_3.sqlite \
                -i Redekar_et_al_2019_sup_table_3.fasta \
                --left GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA \
                --right AGCGTTCTTCATCGATGTGC --marker ITS1-long
fi

echo "Drawing edit-graph for database entries alone"
# Using -k / --marker to show the DB entries too
thapbi_pict edit-graph -k ITS1-long \
        -d Redekar_et_al_2019_sup_table_3.sqlite \
        -o Redekar_et_al_2019_sup_table_3.xgmml

echo "Running analysis"
thapbi_pict pipeline \
        -i raw_data/ -s intermediate_long/ -o summary/recycled-water-custom \
        --merged-cache tmp_merged/ \
        -d Redekar_et_al_2019_sup_table_3.sqlite -m onebp \
        -t metadata.tsv -x 7 -c 1,2,3,4,5,6

echo ====
echo Done
echo ====
