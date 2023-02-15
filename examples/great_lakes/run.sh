#!/bin/bash
set -euo pipefail

echo "NOTE: Expected first time run time is about 40 minutes,"
echo "repeat runs are under a minute (regenerating reports etc)."

mkdir -p tmp_merged/ intermediate/ summary/

# Takes arguments via variable names
function import_marker {
    # Pre-trimmed, not validating species names
    thapbi_pict import -d pooled.sqlite \
        -i $NAME.fasta -x -s ";" \
        -k $NAME --left $LEFT --right $RIGHT
}

if [ ! -f pooled.sqlite ]; then
    echo =================
    echo Creating database
    echo =================

    # MOL16S primers, 183–310 bp
    NAME=MOL16S
    LEFT=RRWRGACRAGAAGACCCT
    RIGHT=ARTCCAACATCGAGGT
    import_marker  # calls function defined above

    # SPH16S primers, 299 bp
    NAME=SPH16S
    LEFT=TAGGGGAAGGTATGAATGGTTTG
    RIGHT=ACATCGAGGTCGCAACC
    import_marker  # calls function defined above
fi

echo ================
echo Running pipeline
echo ================

# Note the unusually minimum low abundance threshold
# of 10 is deliberate. This *does* let unwanted noise
# through - see the discussion in the documentation.
#
# We have different expected/$NAME/*.known.tsv files
# so currently cannot use assess via pipeline...
mkdir -p intermediate/ summary/
thapbi_pict pipeline -d pooled.sqlite --synthetic "" \
    -i raw_data/ -a 10 -f 0 \
    --merged-cache tmp_merged/ \
    -s intermediate/ -o summary/ \
    -t metadata.tsv -x 1 -c 4,5,3,2
#   -t PRJNA379165.tsv -x 1 -c 4,8

echo ===============
echo Running reports
echo ===============

for NAME in MOL16S SPH16S; do
    # Run the assess step per marker... can't do this in the pipeline
    # as need to pass in different expected files per marker.
    thapbi_pict assess -d pooled.sqlite --marker $NAME \
        -i expected/$NAME/ summary/$NAME.tally.tsv \
           summary/$NAME.onebp.tsv \
        -o summary/$NAME.assess.onebp.tsv
done

echo ===========
echo Edit graphs
echo ===========

for NAME in MOL16S SPH16S; do
    echo
    echo "Edit graph of $NAME with -a 100"
    # Run an edit graph at the default -a 100 setting,
    # without showing the DB entries
    thapbi_pict edit-graph -d pooled.sqlite -a 100 \
        -i summary/$NAME.tally.tsv -o summary/$NAME.edit-graph.a100.xgmml
done

echo
echo "Edit graph of just the 7 MOL16S mock community samples"
thapbi_pict sample-tally -d pooled.sqlite -k MOL16S -a 10 -f 0 \
    -i intermediate/MOL16S/SRR5534972.* \
       intermediate/MOL16S/SRR5534973.* \
       intermediate/MOL16S/SRR5534974.* \
       intermediate/MOL16S/SRR5534975.* \
       intermediate/MOL16S/SRR5534976.* \
       intermediate/MOL16S/SRR5534977.* \
       intermediate/MOL16S/SRR5534979.* \
    -o summary/MOL16S.mocks.tally.tsv
thapbi_pict edit-graph -d pooled.sqlite -a 100 \
    -i summary/MOL16S.mocks.tally.tsv \
    -o summary/MOL16S.mocks.edit-graph.a100.xgmml

echo
echo "Edit graph of just the 3 SPH16S mock community samples"
thapbi_pict sample-tally -d pooled.sqlite -k SPH16S -a 10 -f 0 \
    -i intermediate/SPH16S/SRR5534978.* \
       intermediate/SPH16S/SRR5534980.* \
       intermediate/SPH16S/SRR5534981.* \
    -o summary/SPH16S.mocks.tally.tsv
thapbi_pict edit-graph -d pooled.sqlite -a 100 \
    -i summary/SPH16S.mocks.tally.tsv \
    -o summary/SPH16S.mocks.edit-graph.a100.xgmml

echo ====
echo Done
echo ====
