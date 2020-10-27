#!/bin/bash
set -eup pipeline

echo NOTE: Expected first time run time is about 1.5 hours,
echo repeat runs about 1 minute just to regenerate reports.
echo

# Takes arguments via variable names
function analyse {
    echo "Building database for $NAME"
    rm -rf references/${NAME}.sqlite
    # Assume pre-trimmed
    thapbi_pict curated-import -i references/${NAME}.fasta -d references/${NAME}.sqlite -x

    echo "Running analysis for $NAME"
    mkdir -p intermediate/$NAME/
    thapbi_pict prepare-reads --left $LEFT --right $RIGHT \
        -i raw_data/ --merged-cache tmp_merged/ -o intermediate/$NAME/
    thapbi_pict fasta-nr -i intermediate/$NAME/*.fasta -o summary/$NAME.all.fasta
    thapbi_pict classify -i summary/$NAME.all.fasta -o summary/ -d references/$NAME.sqlite
    thapbi_pict pipeline -d references/${NAME}.sqlite --showdb \
        --merged-cache tmp_merged/ \
        --left $LEFT --right $RIGHT \
                -i raw_data/ expected/ \
                -s intermediate/$NAME/ -o summary/ -r $NAME \
                -t metadata.tsv -c 3,4,5 -x 2 -g 4
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
    for S in `cut -f 4 PRJEB18620.txt | grep -v "sample_alias"`; do
        grep "^#" intermediate/16S/$S.fasta | grep -v -E "(_primer|cutadapt|abundance)" > intermediate_pool/$S.fasta
        cat intermediate/*/$S.fasta | grep -v "^#" >> intermediate_pool/$S.fasta
    done

    echo "Computing species list for combined header..."
    # Quick and dirty pooling by concatenating the intermediate
    # classifier result files per sample is enough for sample-summary,
    # but headers need reworking for the assess command to work...
    #
    # This is largely messing about to get tabs in the right places
    echo "#sequence-name" > header.tmp
    echo "taxid" >> header.tmp
    echo -n "genus-species:" >> header.tmp
    # Get all possible species classifications, one per line, sorted and unique, then semi-colon separate them
    cat intermediate/*/*.onebp.tsv | grep -v '^#' | cut -f 3 | tr ';' '\n' | sort | uniq | tr '\n' ';' >> header.tmp
    echo "" >> header.tmp
    tr '\n' '\t' < header.tmp > intermediate_pool/header.onebp.txt
    echo "note" >> intermediate_pool/header.onebp.txt
    rm header.tmp

    echo "Pooling intermediate onebp classifications..."
    for S in `cut -f 4 PRJEB18620.txt | grep -v "sample_alias"`; do
        cp intermediate_pool/header.onebp.txt intermediate_pool/$S.onebp.tsv
        cat intermediate/*/$S.onebp.tsv | grep -v "^#" >> intermediate_pool/$S.onebp.tsv
    done;

    echo "Generating pooled reports for onebp classifier."
    # Now the reports:
    thapbi_pict summary -m onebp -i intermediate_pool/ \
        -o summary/ -r pooled \
        -t metadata.tsv -c 3,4,5 -x 2 -g 4
    # And the assessment
    thapbi_pict assess -i expected/ intermediate_pool/ -m onebp \
        -o summary/pooled.assess.onebp.tsv

    echo "Pooled results done"
}

echo =====================================================
echo Universal animal DNA barcodes and mini-barcodes - 16S
echo =====================================================

NAME=16S
LEFT=CGCCTGTTTATCAAAAACAT
RIGHT=CCGGTCTGAACTCAGATCACGT

analyse # call function above

echo ==========================================================
echo Universal animal DNA barcodes and mini-barcodes - Mini-16S
echo ==========================================================

NAME=Mini-16S
LEFT=AYAAGACGAGAAGACCC
RIGHT=GATTGCGCTGTTATTCC

analyse # call function above

echo ==========================================================
echo Universal animal DNA barcodes and mini-barcodes - Mini-COI
echo ==========================================================

#TODO - COI, multiple right primers?

NAME=Mini-COI
LEFT=GGWACWGGWTGAACWGTWTAYCCYCC
RIGHT=TAIACYTCIGGRTGICCRAARAAYCA

analyse # call function above

echo =======================================================
echo Universal animal DNA barcodes and mini-barcodes - cyt-b
echo =======================================================

NAME=cyt-b
LEFT=CCATCCAACATCTCAGCATGATGAAA
RIGHT=GGCAAATAGGAARTATCATTC

echo "Skipping, failed to amplify at default threhold or even 50."
echo "Drop to -a 10 and you get a modest number of sequences."

#analyse

echo ============================================================
echo Universal animal DNA barcodes and mini-barcodes - Mini-cyt-b
echo ============================================================

NAME=Mini-cyt-b
LEFT=CCATCCAACATCTCAGCATGATGAAA
RIGHT=CCCTCAGAATGATATTTGTCCTCA

analyse # call function above

echo =====================================================
echo Universal plant DNA barcodes and mini-barcodes - matK
echo =====================================================

NAME=matK
LEFT=ACCCAGTCCATCTGGAAATCTTGGTTC
RIGHT=CGTACAGTACTTTTGTGTTTACGAG

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

analyse

echo ==========================================================
echo Universal plant DNA barcodes and mini-barcodes - Mini-rbcL
echo ==========================================================

NAME=Mini-rbcL
LEFT=GTTGGATTCAAAGCTGGTGTTA
RIGHT=CVGTCCAMACAGTWGTCCATGT

analyse # call function above

echo =========================================================
echo Universal plant DNA barcodes and mini-barcodes - trnL-UAA
echo =========================================================

NAME=trnL-UAA
LEFT=CGAAATCGGTAGACGCTACG
RIGHT=GGGGATAGAGGGACTTGAAC

analyse # call function above

echo =============================================================
echo Universal plant DNA barcodes and mini-barcodes - trnL-P6-loop
echo =============================================================

NAME=trnL-P6-loop
LEFT=GGGCAATCCTGAGCCAA
RIGHT=CCATTGAGTCTCTGCACCTATC

echo "Skipping, failed to amplify - no trimmed sequence found"
echo "more than once, only singletons which are as good as noise."
echo "Need to set --min-len 10 to match the author's analysis."

#analyse

echo =====================================================
echo Universal plant DNA barcodes and mini-barcodes - ITS2
echo =====================================================

NAME=ITS2
LEFT=ATGCGATACTTGGTGTGAAT
RIGHT=GACGCTTCTCCAGACTACAAT

analyse # call function above

echo ===============
echo Pooling markers
echo ===============

pool # call function above

echo ====
echo Done
echo ====

echo "You may wish to delete ``tmp_merged/*.fasta.gz`` now (about 1GB)."
