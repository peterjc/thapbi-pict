#!/bin/bash
set -euo pipefail

echo "NOTE: Expected first time run time is about 5 minutes,"
echo "repeat runs under a minute just to regenerate reports."
echo

mkdir -p intermediate/ summary/

# Takes arguments via variable names
function analyse {
    if [ ! -f ${NAME}.sqlite ]; then
        echo "Trimming $MARKER sequences for $NAME"
        export RIGHT_RC=`python -c "from Bio.Seq import reverse_complement as rc; print(rc('$RIGHT'))"`
        # Doing the left and right primer trimming separately:
        cutadapt --quiet -g $LEFT $MARKER.fasta \
          | cutadapt --quiet -a $RIGHT_RC -o $NAME.fasta /dev/stdin
        echo "Building $MARKER database for $NAME"
        thapbi_pict import -k $NAME -i $NAME.fasta -d $NAME.sqlite -x
    fi

    echo "Running analysis with minimum abundance threshold ten"
    # No threshold (-a 0 or -a 1) gives 450k total unique entries over samples
    # Using minimum of 2 gives 75k unique, 5 gives 22k, and 10 gives 8.8k unique.
    # Using minimum of 100 (default) gives under 800 unique over samples.
    # [Counts were over both amplicons using the actual primer pairs, 3 runs]
    mkdir -p intermediate/${LIBRARY}_${NAME}/
    for METHOD in identity onebp blast; do
        thapbi_pict pipeline -d ${NAME}.sqlite \
                    --left $LEFT --right $RIGHT \
                    -i raw_data/$LIBRARY/ expected/$LIBRARY/ -m $METHOD \
                    -s intermediate/${LIBRARY}_${NAME}/ \
                    -o summary/${LIBRARY}_${NAME} -a 10 \
                    -t metadata_$LIBRARY.tsv -c 5,6,7,3,4,2 -x 1 -g 6
    done
    # Now run an edit-graph at a higher abundance threshold
    # (works as long as pipeline or prepare-reads was run with
    # the same or lower threshold).
    # Including all DB entries with -s / --showdb argument
    # Do not show the classifier output using -m with "-"
    thapbi_pict edit-graph -d ${NAME}.sqlite --showdb \
                -i intermediate/${LIBRARY}_${NAME}/ -a 75 -m - \
                -o summary/${LIBRARY}_${NAME}.edit-graph.a75.xgmml
    echo "$NAME done"
}

echo ================================================
echo Amplicon library one - ITS1 - BITS/B58S3 primers
echo ================================================

MARKER=ITS1
LIBRARY=AL1
NAME=BITS_B58S3
LEFT=ACCTGCGGARGGATC
RIGHT=GAGATCCRTTGYTRAAAGTT

analyse # call function above

echo ================================================
echo Amplicon library two - ITS1 - BITS/B58S3 primers
echo ================================================
echo Note: This is a blinkered view of this dataset,
echo really used ITS1f/ITS2 primers which amplify a
echo a larger fragment, see below.

MARKER=ITS1
LIBRARY=AL2
NAME=BITS_B58S3
LEFT=ACCTGCGGARGGATC
RIGHT=GAGATCCRTTGYTRAAAGTT

analyse # call function above

echo ================================================
echo Amplicon library two - ITS1 - ITS1f/ITS2 primers
echo ================================================

MARKER=ITS1
LIBRARY=AL2
NAME=ITS1f_ITS2
LEFT=CTTGGTCATTTAGAGGAAGTAA
RIGHT=GCTGCGTTCTTCATCGATGC

analyse # call function above

echo =========================================================
echo Amplicon library two - ITS2 - ITS3-KYO2/ITS4-KYO3 primers
echo =========================================================

MARKER=ITS2
LIBRARY=AL2
NAME=ITS3-KYO2_ITS4-KYO3
LEFT=GATGAAGAACGYAGYRAA
RIGHT=CTBTTVCCKCTTCACTCG

analyse # call function above

echo ====
echo Done
echo ====
