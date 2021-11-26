#!/bin/bash
set -euo pipefail

echo "NOTE: Expected first time run time is under 5 minutes,"
echo "repeat runs under a minute to just to regenerate reports."
echo

mkdir -p intermediate/ summary/

# Takes arguments via variable names
function import_marker {
    echo "Trimming $GENE sequences for $MARKER"
    export RIGHT_RC=`python -c "from Bio.Seq import reverse_complement as rc; print(rc('$RIGHT'))"`
    # Doing the left and right primer trimming separately:
    cutadapt --quiet -g $LEFT $GENE.fasta \
      | cutadapt --quiet -a $RIGHT_RC -o $MARKER.fasta /dev/stdin
    echo "Adding $GENE $MARKER to database"
    thapbi_pict import -d $DB -i $MARKER.fasta  -x \
                -k $MARKER --left $LEFT --right $RIGHT
}

# Takes arguments via variable names
function analyse {
    echo "Running analysis with minimum abundance threshold ten"
    # No threshold (-a 0 or -a 1 with -f 0) gives 450k total unique entries over samples
    # Using minimum of 2 gives 75k unique, 5 gives 22k, and 10 gives 8.8k unique.
    # Using minimum of 100 (default) gives under 800 unique over samples.
    # [Counts were over both amplicons using the actual primer pairs, 3 runs]
    mkdir -p intermediate/${LIBRARY}/
    for METHOD in identity onebp blast; do
        thapbi_pict pipeline -d $DB -m $METHOD \
                    -i raw_data/$LIBRARY/ expected/$LIBRARY/ \
                    -s intermediate/${LIBRARY}/ \
                    -o summary/${LIBRARY} -a 10 -f 0 \
                    -t metadata_$LIBRARY.tsv -c 5,6,7,3,4,2 -x 1 -g 6
    done
    echo "$LIBRARY done"
}

# Takes arguments via variable names
function edit_graph {
    # Now run an edit-graph at a higher abundance threshold
    # (works as long as pipeline or prepare-reads was run with
    # the same or lower threshold).
    # Including relevant DB entries with -k / --marker argument
    # Do not show the classifier output using -m with "-"
    thapbi_pict edit-graph -d $DB -k $MARKER \
                -i intermediate/${LIBRARY}/${MARKER}/ -a 75 -m - \
                -o summary/${LIBRARY}.${MARKER}.edit-graph.a75.xgmml
}


DB=fungi_solo.sqlite
if [ ! -f $DB ]; then
    echo ==================================================
    echo Creating database - ITS1 - BITS/B58S3 primers only
    echo ==================================================

    GENE=ITS1
    MARKER=BITS-B58S3
    LEFT=ACCTGCGGARGGATC
    RIGHT=GAGATCCRTTGYTRAAAGTT
    import_marker  # call function above
fi

echo =====================================================
echo Amplicon library one - ITS1 - BITS/B58S3 primers only
echo =====================================================

LIBRARY=AL1
analyse  # call function above

MARKER=BITS-B58S3
edit_graph  # call function abover

echo ================================================
echo Amplicon library two - ITS1 - BITS/B58S3 primers
echo ================================================
echo Note: This is a blinkered view of this dataset,
echo really used ITS1f/ITS2 primers which amplify a
echo a larger fragment, see below.

LIBRARY=AL2
analyse  # call function above

MARKER=BITS-B58S3
edit_graph  # call function abover

DB=fungi_duo.sqlite
if [ ! -f $DB ]; then
    echo =================================
    echo Creating database - ITS1 and ITS2
    echo =================================

    GENE=ITS1
    MARKER=ITS1f-ITS2
    LEFT=CTTGGTCATTTAGAGGAAGTAA
    RIGHT=GCTGCGTTCTTCATCGATGC
    import_marker  # call function above

    GENE=ITS2
    LIBRARY=AL2
    MARKER=ITS3-KYO2-ITS4-KYO3
    LEFT=GATGAAGAACGYAGYRAA
    RIGHT=CTBTTVCCKCTTCACTCG
    import_marker  # call function above
fi

echo ================================================
echo Amplicon library two, with two primers products:
echo ITS1 - ITS1f/ITS2 primers
echo ITS2 - ITS3-KYO2/ITS4-KYO3 primers
echo ================================================

LIBRARY=AL2
analyse  # call function above

for MARKER in ITS1f-ITS2 ITS3-KYO2-ITS4-KYO3; do
    edit_graph  # call function above
done

echo ====
echo Done
echo ====
