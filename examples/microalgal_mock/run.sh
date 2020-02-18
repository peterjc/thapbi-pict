#!/bin/bash
set -eup pipeline

echo NOTE: Expected first time run time is about 10 minutes,
echo repeat runs about 1 minute just to regenerate reports.
echo

# Takes arguments via variable names
function analyse {
    echo "Building $NAME database for $NAME"
    rm -rf ${NAME}.sqlite
    # FASTA file has full 18S rRNA gene, use primers to trim to targetted region:
    thapbi_pict curated-import -i mock_community.fasta -d ${NAME}.sqlite --left $LEFT --right $RIGHT -x

    echo "Running analysis"
    mkdir -p $NAME/intermediate/
    # thapbi_pict prepare-reads -i $NAME/raw_data/ -o $NAME/intermediate/ --left "" --right ""
    # Assume FASTQ already have primers removed!
    thapbi_pict pipeline -d ${NAME}.sqlite --left "" --right "" \
                -i $NAME/raw_data/ $NAME/expected/ \
                -s $NAME/intermediate/ -o $NAME/ -r $NAME \
                -t $NAME/metadata.tsv -c 3,4 -x 1
    echo "$NAME done"
}

echo ============================
echo V4 - Reuk454FWD1/V4r primers
echo ============================

NAME=V4
LEFT=CCAGCASCYGCGGTAATTCC
RIGHT=ACTTTCGTTCTTGAT

analyse # call function above

echo =========================
echo V8-V9 - V8f/1510r primers
echo =========================

NAME=V8V9
LEFT=ATAACAGGTCTGTGATGCCCT
RIGHT=CCTTCYGCAGGTTCACCTAC

analyse # call function above
