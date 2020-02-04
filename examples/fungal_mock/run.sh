#!/bin/bash
set -eup pipeline

# Takes arguments via variable names
function analyse {
    echo "Building $MARKER database for $NAME"
    rm -rf ${NAME}.sqlite
    thapbi_pict curated-import -i ${MARKER}.fasta -d ${NAME}.sqlite --left $LEFT --right $RIGHT -x

    echo "Running analysis with default abundance threshold"
    mkdir -p $LIBRARY/$NAME
    thapbi_pict pipeline -d ${NAME}.sqlite --left $LEFT --right $RIGHT \
		-i $LIBRARY/raw_data/ $LIBRARY/expected/ \
		-s $LIBRARY/$NAME -o $LIBRARY/ -r $NAME \
		-t $LIBRARY/metadata.tsv -c 5,6,7,3,4,2 -x 1 -g 6
    echo "$NAME done"
}

echo ================================================
echo Amplicon library one - ITS1 - BITS/B58S3 primers
echo ================================================

MARKER=ITS1
LIBRARY=amp_lib_one
NAME=BITS_B58S3
LEFT=ACCTGCGGARGGATC
RIGHT=GAGATCCRTTGYTRAAAGTT

analyse # call function above

echo ================================================
echo Amplicon library two - ITS1 - ITS1f/ITS2 primers
echo ================================================

MARKER=ITS1
LIBRARY=amp_lib_two
NAME=ITS1f_ITS2
LEFT=CTTGGTCATTTAGAGGAAGTAA
RIGHT=GCTGCGTTCTTCATCGATGC

analyse # call function above

echo =========================================================
echo Amplicon library two - ITS2 - ITS3-KYO2/ITS4-KYO3 primers
echo =========================================================

MARKER=ITS2
LIBRARY=amp_lib_two
NAME=ITS3-KYO2_ITS4-KYO3
LEFT=GATGAAGAACGYAGYRAA
RIGHT=CTBTTVCCKCTTCACTCG

analyse # call function above
