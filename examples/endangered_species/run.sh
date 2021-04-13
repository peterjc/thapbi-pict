#!/bin/bash
set -eup pipeline

echo NOTE: Expected first time run time is about 1.5 hours,
echo repeat runs about 1 minute just to regenerate reports.
echo

# Quoting paper on minimum lengths:
#
#     We implemented a minimum DNA barcode length of 200 nt, except for DNA
#     barcodes with a basic length shorter than 200 nt, in which case the
#     minimum expected DNA barcode length is set to 100 nt for ITS2, 140 nt
#     for mini-rbcL, and 10 nt for the trnL (P6 loop) marker.

rm -rf references/pooled.sqlite

# Takes arguments via variable names
function analyse {
    if [ ! -f references/${NAME}.sqlite ]; then
        echo "Building database for $NAME"
        # Assume pre-trimmed
        thapbi_pict curated-import -i references/${NAME}.fasta -d references/${NAME}.sqlite \
            -x --minlen $MINLEN
    fi
    # This is a big hack, just for the assess command to work on the pool:
    thapbi_pict curated-import -i references/${NAME}.fasta -d references/pooled.sqlite \
        -x --minlen $MINLEN

    echo "Running analysis for $NAME"
    mkdir -p intermediate/$NAME/
    thapbi_pict pipeline -d references/${NAME}.sqlite --showdb \
        -i raw_data/ expected/ -s intermediate/$NAME/ --merged-cache tmp_merged/ \
        -o summary/ -r $NAME --left $LEFT --right $RIGHT --minlen $MINLEN \
        -t metadata.tsv -c 3,4,5 -x 2 -g 4
    # Pipeline now includes the fasta-nr step:
    # thapbi_pict fasta-nr -i intermediate/$NAME/*.fasta -o summary/$NAME.all_reads.fasta
    thapbi_pict classify -i summary/$NAME.all_reads.fasta -o summary/ -d references/$NAME.sqlite
    echo "$NAME done"
}

function pool {
    # This is a hack to make up for limited functionality within
    # THAPBI PICT itself; we may add a 'pool' command to the tool
    # to do this directly.
    rm -rf intermediate_pool/*
    echo "Pooling intermediate FASTA files..."
    # Excluding primer specific header lines with grep,
    # only want a single header
    for S in `cut -f 4 PRJEB18620.tsv | grep -v "sample_alias"`; do
        grep "^#" intermediate/16S/$S.fasta | grep -v -E "(_primer|cutadapt|abundance)" > intermediate_pool/$S.fasta
        cat intermediate/*/$S.fasta | grep -v "^#" >> intermediate_pool/$S.fasta
    done

    echo "Computing species list for combined header..."
    # Quick and dirty pooling by concatenating the intermediate
    # classifier result files per sample is enough for sample-summary,
    echo "Pooling intermediate onebp classifications..."
    for S in `cut -f 4 PRJEB18620.tsv | grep -v "sample_alias"`; do
        echo -e "#sequence-name\ttaxid\tgenus-species\tnote" > intermediate_pool/$S.onebp.tsv
        cat intermediate/*/$S.onebp.tsv | grep -v "^#" >> intermediate_pool/$S.onebp.tsv
    done;

    echo "Generating pooled reports for onebp classifier."
    # Now the reports:
    thapbi_pict summary -m onebp -i intermediate_pool/ \
        -o summary/ -r pooled \
        -t metadata.tsv -c 3,4,5 -x 2 -g 4

    # assessment... as of v0.8.1 need a DB giving possible species
    thapbi_pict assess -i expected/ intermediate_pool/ -m onebp \
         -d references/pooled.sqlite -o summary/pooled.assess.onebp.tsv

    echo "Pooled results done"
}

echo =====================================================
echo Universal animal DNA barcodes and mini-barcodes - 16S
echo =====================================================

NAME=16S
LEFT=CGCCTGTTTATCAAAAACAT
RIGHT=CCGGTCTGAACTCAGATCACGT
MINLEN=200

analyse # call function above

echo ==========================================================
echo Universal animal DNA barcodes and mini-barcodes - Mini-16S
echo ==========================================================

NAME=Mini-16S
LEFT=AYAAGACGAGAAGACCC
RIGHT=GATTGCGCTGTTATTCC
MINLEN=200

analyse # call function above

echo ==========================================================
echo Universal animal DNA barcodes and mini-barcodes - Mini-COI
echo ==========================================================

#TODO - COI, multiple right primers?

NAME=Mini-COI
LEFT=GGWACWGGWTGAACWGTWTAYCCYCC
RIGHT=TAIACYTCIGGRTGICCRAARAAYCA
MINLEN=200

analyse # call function above

echo =======================================================
echo Universal animal DNA barcodes and mini-barcodes - cyt-b
echo =======================================================

NAME=cyt-b
LEFT=CCATCCAACATCTCAGCATGATGAAA
RIGHT=GGCAAATAGGAARTATCATTC
MINLEN=200

echo "Skipping, failed to amplify at default threhold or even 50."
echo "Drop to -a 10 and you get a modest number of sequences."

#analyse

echo ============================================================
echo Universal animal DNA barcodes and mini-barcodes - Mini-cyt-b
echo ============================================================

NAME=Mini-cyt-b
LEFT=CCATCCAACATCTCAGCATGATGAAA
RIGHT=CCCTCAGAATGATATTTGTCCTCA
MINLEN=200

analyse # call function above

echo =====================================================
echo Universal plant DNA barcodes and mini-barcodes - matK
echo =====================================================

NAME=matK
LEFT=ACCCAGTCCATCTGGAAATCTTGGTTC
RIGHT=CGTACAGTACTTTTGTGTTTACGAG
MINLEN=200

echo "Skipping, failed to amplify"
echo "(at least at default threshold)"

#analyse

#Note authors excluded the other matK primers...

echo =====================================================
echo Universal plant DNA barcodes and mini-barcodes - rbcL
echo =====================================================

NAME=rbcL
LEFT=ATGTCACCACAAACAGAGACTAAAGC
RIGHT=GTAAAATCAAGTCCACCRCG
MINLEN=100  # Even shorter than author's 200 default

analyse

echo ==========================================================
echo Universal plant DNA barcodes and mini-barcodes - Mini-rbcL
echo ==========================================================

NAME=Mini-rbcL
LEFT=GTTGGATTCAAAGCTGGTGTTA
RIGHT=CVGTCCAMACAGTWGTCCATGT
MINLEN=140  # Shorter!

analyse # call function above

echo =========================================================
echo Universal plant DNA barcodes and mini-barcodes - trnL-UAA
echo =========================================================

NAME=trnL-UAA
LEFT=CGAAATCGGTAGACGCTACG
RIGHT=GGGGATAGAGGGACTTGAAC
MINLEN=200

analyse # call function above

echo =============================================================
echo Universal plant DNA barcodes and mini-barcodes - trnL-P6-loop
echo =============================================================

NAME=trnL-P6-loop
LEFT=GGGCAATCCTGAGCCAA
RIGHT=CCATTGAGTCTCTGCACCTATC
MINLEN=10  # Much shorter!

analyse

echo =====================================================
echo Universal plant DNA barcodes and mini-barcodes - ITS2
echo =====================================================

NAME=ITS2
LEFT=ATGCGATACTTGGTGTGAAT
RIGHT=GACGCTTCTCCAGACTACAAT
MINLEN=100  # Shorter!

analyse # call function above

echo ===============
echo Pooling markers
echo ===============

pool # call function above

echo ====
echo Done
echo ====

echo "You may wish to delete ``tmp_merged/*.fasta.gz`` now (about 1GB)."
