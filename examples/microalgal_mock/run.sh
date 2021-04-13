#!/bin/bash
set -eup pipeline

echo NOTE: Expected first time run time is about 10 minutes,
echo repeat runs about 1 minute just to regenerate reports.
echo

# Takes arguments via variable names
function analyse {
    if [ ! -f ${NAME}.sqlite ]; then
        echo "Building $NAME database for $NAME"
        # FASTA file has full 18S rRNA gene, use primers to trim to targetted region:
        thapbi_pict curated-import -i mock_community.fasta -d ${NAME}.sqlite --left $LEFT --right $RIGHT -x
    fi

    echo "Running analysis"
    mkdir -p intermediate/$NAME/
    # Assume FASTQ already have primers removed!
    thapbi_pict pipeline -d ${NAME}.sqlite --left "" --right "" \
                -i raw_data/$NAME/ expected/$NAME/ \
                -s intermediate/$NAME/ -o summary/ -r $NAME \
                -t metadata.tsv -c 1,2,3,4,5 -x $ID_COL
    echo "$NAME done"
}

echo ============================
echo V4 - Reuk454FWD1/V4r primers
echo ============================

NAME=V4
ID_COL=6
LEFT=CCAGCASCYGCGGTAATTCC
RIGHT=ACTTTCGTTCTTGAT

analyse # call function above

echo =========================
echo V8-V9 - V8f/1510r primers
echo =========================

NAME=V8V9
ID_COL=7
LEFT=ATAACAGGTCTGTGATGCCCT
RIGHT=CCTTCYGCAGGTTCACCTAC

analyse # call function above

echo ====
echo Done
echo ====
