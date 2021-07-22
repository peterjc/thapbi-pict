#!/bin/bash
set -eup pipeline

echo NOTE: Expected first time run time is under 5 minutes,
echo repeat runs under 1 minute just to regenerate reports.
echo
echo ===================
echo British soil - ITS1
echo ===================

# These primers are the same as our defaults, so can use our default DB.
# Negative controls give confidence we can lower the default min abundance.

thapbi_pict pipeline -i raw_data/ expected/ -s intermediate/ \
            --ignore-prefixes Unavailable \
            -o summary/british_soil -a 50 \
            -t metadata.tsv -c 1,2,3,4,5,6,7 -x 8

echo =================
echo Positive Controls
echo =================

mkdir -p controls
thapbi_pict pipeline -i raw_data/SRR13393802_* raw_data/SRR13393813_* \
            raw_data/SRR13393837_* expected/ -o controls/controls-only -a 1

echo ----------------------------------
echo False negatives from P. boehmeriae
echo ----------------------------------
grep -E "(predictions|boehmeriae)" controls/controls-only.ITS1.reads.onebp.tsv | cut -f 1,2,5

echo -----------------------------
echo False negatives from P. idaei
echo -----------------------------
grep -E "(predictions|idaei)" controls/controls-only.ITS1.reads.onebp.tsv | cut -f 1,2,7- | head

echo ====
echo Done
echo ====
