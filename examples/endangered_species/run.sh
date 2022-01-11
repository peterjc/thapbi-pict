#!/bin/bash
set -euo pipefail

echo "NOTE: Expected first time run time is about 1.5 hours,"
echo "repeat runs under a minute to regenerate reports."
echo

mkdir -p tmp_merged/ intermediate/ intermediate_pool/ summary/

# Quoting paper on minimum lengths:
#
#     We implemented a minimum DNA barcode length of 200 nt, except for DNA
#     barcodes with a basic length shorter than 200 nt, in which case the
#     minimum expected DNA barcode length is set to 100 nt for ITS2, 140 nt
#     for mini-rbcL, and 10 nt for the trnL (P6 loop) marker.

function import_marker {
    # Takes arguments via variable names
    thapbi_pict import -d references/pooled.sqlite \
                -i references/${NAME}.fasta -x -s ";" \
                -k $NAME --left $LEFT --right $RIGHT --minlen $MINLEN
}

if [ ! -f references/pooled.sqlite ]; then
    echo =================
    echo Creating database
    echo =================

    # Universal animal DNA barcodes and mini-barcodes - 16S
    NAME=16S
    LEFT=CGCCTGTTTATCAAAAACAT
    RIGHT=CCGGTCTGAACTCAGATCACGT
    MINLEN=200
    import_marker  # calls function defined above

    # Universal animal DNA barcodes and mini-barcodes - Mini-16S
    NAME=Mini-16S
    LEFT=AYAAGACGAGAAGACCC
    RIGHT=GATTGCGCTGTTATTCC
    MINLEN=200
    import_marker  # calls function defined above

    # Universal animal DNA barcodes and mini-barcodes - Mini-COI
    NAME=Mini-COI
    LEFT=GGWACWGGWTGAACWGTWTAYCCYCC
    RIGHT=TAIACYTCIGGRTGICCRAARAAYCA
    MINLEN=200
    import_marker  # calls function defined above

    #TODO - COI, multiple right primers?

    # Universal animal DNA barcodes and mini-barcodes - cyt-b
    NAME=cyt-b
    LEFT=CCATCCAACATCTCAGCATGATGAAA
    RIGHT=GGCAAATAGGAARTATCATTC
    MINLEN=200
    # Skipping, failed to amplify at default threshold or even 50.
    # Drop to -a 10 and you get a modest number of sequences.

    # Universal animal DNA barcodes and mini-barcodes - Mini-cyt-b
    NAME=Mini-cyt-b
    LEFT=CCATCCAACATCTCAGCATGATGAAA
    RIGHT=CCCTCAGAATGATATTTGTCCTCA
    MINLEN=200
    import_marker  # calls function defined above

    # Universal plant DNA barcodes and mini-barcodes - matK
    NAME=matK
    LEFT=ACCCAGTCCATCTGGAAATCTTGGTTC
    RIGHT=CGTACAGTACTTTTGTGTTTACGAG
    MINLEN=200
    # Skipping, failed to amplify (at least at default threshold)
    # Note authors excluded the other matK primers...

    # Universal plant DNA barcodes and mini-barcodes - rbcL
    NAME=rbcL
    LEFT=ATGTCACCACAAACAGAGACTAAAGC
    RIGHT=GTAAAATCAAGTCCACCRCG
    MINLEN=100  # Even shorter than author's 200 default
    import_marker  # calls function defined above

    # Universal plant DNA barcodes and mini-barcodes - Mini-rbcL
    NAME=Mini-rbcL
    LEFT=GTTGGATTCAAAGCTGGTGTTA
    RIGHT=CVGTCCAMACAGTWGTCCATGT
    MINLEN=140  # Shorter!
    import_marker  # calls function defined above

    # Universal plant DNA barcodes and mini-barcodes - trnL-UAA
    NAME=trnL-UAA
    LEFT=CGAAATCGGTAGACGCTACG
    RIGHT=GGGGATAGAGGGACTTGAAC
    MINLEN=200
    import_marker  # calls function defined above

    # Universal plant DNA barcodes and mini-barcodes - trnL-P6-loop
    NAME=trnL-P6-loop
    LEFT=GGGCAATCCTGAGCCAA
    RIGHT=CCATTGAGTCTCTGCACCTATC
    MINLEN=10  # Much shorter!
    import_marker  # calls function defined above

    # Universal plant DNA barcodes and mini-barcodes - ITS2
    NAME=ITS2
    LEFT=ATGCGATACTTGGTGTGAAT
    RIGHT=GACGCTTCTCCAGACTACAAT
    MINLEN=100  # Shorter!
    import_marker  # calls function defined above
fi

echo ================
echo Running pipeline
echo ================

thapbi_pict pipeline -d references/pooled.sqlite --synthetic "" \
            -i raw_data/ expected/ -s intermediate/ --merged-cache tmp_merged/ \
            -o summary/ -t metadata.tsv -c 3,4,5 -x 2 -g 4

echo ====
echo Done
echo ====

echo "You may wish to delete ``tmp_merged/*.fasta.gz`` now (about 1GB)."
