#!/bin/bash
set -euo pipefail

echo "NOTE: Expected first time run time is about 5 minutes,"
echo "repeat runs under a minute just to regenerate reports."
echo
echo ==================
echo Fecal sequel - COI
echo ==================

mkdir -p intermediate/ summary/

# Setup two database, FASTA files are both already primer trimmed
if [ ! -f COI_430_bats.sqlite ]; then
    echo "Setting up first database"
    thapbi_pict import -k COI \
        --left GTHACHGCYCAYGCHTTYGTAATAAT --right CTCCWGCRTGDGCWAGRTTTCC \
        -d COI_430_bats.sqlite -i COI_430_bats.fasta -x -s ";"
fi
if [ ! -f COI_ext_bats.sqlite ]; then
    echo "Setting up extended database"
    cp COI_430_bats.sqlite COI_ext_bats.sqlite
    thapbi_pict import -k COI \
        --left GTHACHGCYCAYGCHTTYGTAATAAT --right CTCCWGCRTGDGCWAGRTTTCC \
        -d COI_ext_bats.sqlite -i observed_3_bats.fasta -x -s ";"
fi

echo ---------------------------------------------------------------
echo Fecal sequel - COI - Mock community using just 430 reference DB
echo ---------------------------------------------------------------

# Primer pair SFF_145f (GTHACHGCYCAYGCHTTYGTAATAAT) and SFF_351r (CTCCWGCRTGDGCWAGRTTTCC)
mkdir -p intermediate/COI_430_bats/
thapbi_pict pipeline -d COI_430_bats.sqlite \
    -i raw_data/ expected/ -s intermediate/ -o summary/430_bats \
    -t metadata.tsv -x 1 -c 2,3,4

# Default edit-graph has very few DB nodes, so run another edit-graph
# including all DB marker entries too
thapbi_pict edit-graph -d COI_430_bats.sqlite -k COI \
    -i summary/430_bats.COI.onebp.tsv \
    -o summary/430_bats.COI.edit-graph.xgmml

echo ---------------------------------------------------------------
echo Fecal sequel - COI - Mock community using extended reference DB
echo ---------------------------------------------------------------

# The FASTA intermediate files are the same, no point regenerating.
# Currently no TSV files are kept there so using same path with -s

# Primer pair SFF_145f (GTHACHGCYCAYGCHTTYGTAATAAT) and SFF_351r (CTCCWGCRTGDGCWAGRTTTCC)
thapbi_pict pipeline -d COI_ext_bats.sqlite \
    -i raw_data/ expected/ -s intermediate/ \
    -o summary/ext_bats \
    -t metadata.tsv -x 1 -c 2,3,4

echo ====
echo Done
echo ====
