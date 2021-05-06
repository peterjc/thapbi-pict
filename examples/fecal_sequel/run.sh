#!/bin/bash
set -euo pipefail

echo "NOTE: Expected first time run time is about 5 minues,"
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
# Default edit-graph has very few DB nodes, so using --showdb argument
mkdir -p intermediate/COI_430_bats/
thapbi_pict pipeline -i raw_data/ expected/ -s intermediate/COI_430_bats/ \
            -o summary/ -r mock-community.COI_430_bats \
            -d COI_430_bats.sqlite -t metadata.tsv -x 1 -c 2,3,4 \
            --left GTHACHGCYCAYGCHTTYGTAATAAT --right CTCCWGCRTGDGCWAGRTTTCC

# Default edit-graph has very few DB nodes, so run another edit-graph
# including all DB entries with -s / --showdb argument
thapbi_pict edit-graph -d COI_430_bats.sqlite \
            -i intermediate/COI_430_bats/*.fasta \
            summary/mock-community.COI_430_bats.all_reads.onebp.tsv \
            -s -o mock-community.COI_430_bats.edit-graph.xgmml

echo ---------------------------------------------------------------
echo Fecal sequel - COI - Mock community using extended reference DB
echo ---------------------------------------------------------------

# The FASTA intermediate files are the same, no point regenerating them...
mkdir -p intermediate/COI_ext_bats/
cd intermediate/COI_ext_bats/
for FASTA in ../COI_430_bats/*.fasta; do ln -f -s $FASTA; done
cd ../../

# Primer pair SFF_145f (GTHACHGCYCAYGCHTTYGTAATAAT) and SFF_351r (CTCCWGCRTGDGCWAGRTTTCC)
thapbi_pict pipeline -i raw_data/ expected/ -s intermediate/COI_ext_bats/ \
            -o summary/ -r mock-community.COI_ext_bats \
            -d COI_ext_bats.sqlite -t metadata.tsv -x 1 -c 2,3,4 \
            --left GTHACHGCYCAYGCHTTYGTAATAAT --right CTCCWGCRTGDGCWAGRTTTCC

echo ====
echo Done
echo ====
