#!/bin/bash
set -eup pipeline

echo "NOTE: Expected first time run time is about 40 minutes,"
echo "repeat runs under 5 minutes (most of which is rebuilding"
echo "the edit-graphs)."

mkdir -p intermediate/ summary/

# Takes arguments via variable names
function analyse {
    if [ ! -f ${NAME}.sqlite ]; then
        echo "Building $NAME database"
        # Pre-trimmed, not validating species names
        thapbi_pict curated-import -d ${NAME}.sqlite \
                    -i $NAME.fasta -x
    fi

    echo "Running analysis"
    # Note the unusually minimum low abundance threshold
    # of 10 is deliberate. This *does* let unwanted noise
    # through - see the discussion in the documentation.
    mkdir -p intermediate/$NAME/
    thapbi_pict pipeline -d ${NAME}.sqlite --showdb \
                --left $LEFT --right $RIGHT -a 10 \
                -i raw_data/ expected/$NAME/ \
                -s intermediate/$NAME/ -o summary/ -r $NAME \
                -t metadata.tsv -x 1 -c 4,5,3,2
    #           -t PRJNA379165.tsv -x 1 -c 4,8

    # Run an edit graph at the default -a 100 setting, without
    # the --showdb setting (most of the DB content )
    thapbi_pict edit-graph -d ${NAME}.sqlite \
                -i intermediate/$NAME/ -a 100 \
                -o summary/$NAME.edit-graph.a100.xgmml

    echo "$NAME done"
}

echo ==========================
echo MOL16S primers, 183–310 bp
echo ==========================

NAME=MOL16S
LEFT=RRWRGACRAGAAGACCCT
RIGHT=ARTCCAACATCGAGGT

analyse # call function above

#Edit graph of just the mock community samples:
thapbi_pict edit-graph -d ${NAME}.sqlite -a 100 \
            -i intermediate/MOL16S/SRR5534972.* \
               intermediate/MOL16S/SRR5534973.* \
               intermediate/MOL16S/SRR5534974.* \
               intermediate/MOL16S/SRR5534975.* \
               intermediate/MOL16S/SRR5534976.* \
               intermediate/MOL16S/SRR5534977.* \
               intermediate/MOL16S/SRR5534979.* \
            -o summary/$NAME.edit-graph.a100.mock.xgmml

echo ======================
echo SPH16S primers, 299 bp
echo ======================

NAME=SPH16S
LEFT=TAGGGGAAGGTATGAATGGTTTG
RIGHT=ACATCGAGGTCGCAACC

analyse # call function above

# Edit graph of just the mock community samples:
thapbi_pict edit-graph -d ${NAME}.sqlite -a 100 \
                -i intermediate/SPH16S/SRR5534978.* \
                   intermediate/SPH16S/SRR5534980.* \
                   intermediate/SPH16S/SRR5534981.* \
                -o summary/$NAME.edit-graph.a100.mock.xgmml

echo ==============================
echo Mixed primers for long product
echo ==============================

mkdir -p intermediate/large/

thapbi_pict prepare-reads -a 10 \
            --left TAGGGGAAGGTATGAATGGTTTG \
            --right ARTCCAACATCGAGGT \
            -i raw_data -o intermediate/large/

echo ====
echo Done
echo ====
