#!/bin/bash
set -eup pipeline

echo NOTE: Expected first time run time is about 30 minutes,
echo repeat runs about 1 minute just to regenerate reports.
echo

mkdir -p intermediate/ summary/

# Takes arguments via variable names
function analyse {
    echo "Building $NAME database"
    rm -rf ${NAME}.sqlite

    # Pre-trimmed, not validating species names
    thapbi_pict curated-import -d ${NAME}.sqlite \
                -i $NAME.fasta -x

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
    #           -t PRJNA379165.txt -x 1 -c 4,8
    echo "$NAME done"
}

echo ==========================
echo MOL16S primers, 183–310 bp
echo ==========================

NAME=MOL16S
LEFT=RRWRGACRAGAAGACCCT
RIGHT=ARTCCAACATCGAGGT

analyse # call function above

echo ======================
echo SPH16S primers, 299 bp
echo ======================

NAME=SPH16S
LEFT=TAGGGGAAGGTATGAATGGTTTG
RIGHT=ACATCGAGGTCGCAACC

analyse # call function above

echo ==============================
echo Mixed primers for long product
echo ==============================

mkdir -p intermediate/large/

thapbi_pict prepare-reads -a 10 \
            --left TAGGGGAAGGTATGAATGGTTTG \
            --right ARTCCAACATCGAGGT \
            -i raw_data -o intermediate/large/
