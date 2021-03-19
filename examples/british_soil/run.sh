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

thapbi_pict pipeline -i raw_data/ -s intermediate/ \
            -o summary/ -r british_soil_its1 -a 50 \
            -t metadata.tsv -c 1,2,3,4,5,6,7 -x 8
