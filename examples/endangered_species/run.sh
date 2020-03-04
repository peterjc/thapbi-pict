#!/bin/bash
set -eup pipeline

echo NOTE: Expected first time run time is about XXX hours,
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
    thapbi_pict prepare-reads -i raw_data/ -o intermediate/$NAME/ --left $LEFT --right $RIGHT
    thapbi_pict fasta-nr -i intermediate/$NAME/*.fasta -o summary/$NAME.all.fasta
    thapbi_pict classify -i summary/$NAME.all.fasta -o summary/ -d references/$NAME.sqlite
    for METHOD in identity onebp blast; do
        thapbi_pict pipeline -d references/${NAME}.sqlite --left $LEFT --right $RIGHT \
                    -i raw_data/ expected/ -m $METHOD \
                    -s intermediate/$NAME/ -o summary/ -r $NAME \
                    -t metadata.tsv -c 3,4,5 -x 2 -g 4
    done
    # Now run an edit-graph at a higher abundance threshold
    # (works as long as pipeline or prepare-reads was run with
    # the same or lower threshold).
    # Including all DB entries with -s / --showdb argument
    #thapbi_pict edit-graph -d ${NAME}.sqlite -i $NAME/ -s \
    #		-o $NAME.edit-graph.a75.xgmml -a 75
    echo "$NAME done"
}

function pool {
    # Hard coded list of markers, ignores those know to give zero results
    # at default threshold (not being run, will have no output).
    # Maybe "cat intermediate/*/$S.fasta > intermediate_pool/$S.fasta"
    # and "cat intermediate/*/$S.$M.tsv > intermediate_pool/$S.$M.tsv"
    for S in `cut -f 4 PRJEB18620.txt | grep -v "sample_alias"`; do
	cat intermediate/16S/$S.fasta intermediate/Mini-16S/$S.fasta \
	    intermediate/Mini-COI/$S.fasta intermediate/Mini-cyt-b/$S.fasta \
	    intermediate/Mini-rbcL/$S.fasta intermediate/trnL-UAA/$S.fasta \
	    intermediate/ITS2/$S.fasta > intermediate_pool/$S.fasta;
    done
    for M in identity onebp blast; do
	# Quick and dirty pooling by concatenating the intermediate
	# classifier result files per sample (enough for sample-summary,
	# headers need reworking for the assess command to work):
	for S in `cut -f 4 PRJEB18620.txt | grep -v "sample_alias"`; do
	    cat intermediate/16S/$S.$M.tsv intermediate/Mini-16S/$S.$M.tsv \
		intermediate/Mini-COI/$S.$M.tsv intermediate/Mini-cyt-b/$S.$M.tsv\
		intermediate/Mini-rbcL/$S.$M.tsv intermediate/trnL-UAA/$S.$M.tsv \
		intermediate/ITS2/$S.$M.tsv > intermediate_pool/$S.$M.tsv;	    
	done;
	# Now the reports:
	thapbi_pict sample-summary -i intermediate_pool/ \
		    -o summary/pooled.samples.$M.tsv \
		    -r summary/pooled.samples.$M.txt \
		    -e summary/pooled.samples.$M.xlsx \
		    -t metadata.tsv -c 3,4,5 -x 2 -g 4
	thapbi_pict read-summary -i intermediate_pool/ \
                    -o summary/pooled.reads.$M.tsv \
                    -e summary/pooled.reads.$M.xlsx \
                    -t metadata.tsv -c 3,4,5 -x 2 -g 4
	# Can we do assessment too?
    done
}

pool
false

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

echo "Skipping, failed to amplify"
echo "(at least at default threshold)"
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

echo "Skipping, failed to amplify"
echo "(at least at default threshold)"
#analyse

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

echo "Skipping, failed to amplify"
echo "(at least at default threshold)"
#analyse

echo =====================================================
echo Universal plant DNA barcodes and mini-barcodes - ITS2
echo =====================================================

NAME=ITS2
LEFT=ATGCGATACTTGGTGTGAAT
RIGHT=GACGCTTCTCCAGACTACAAT

analyse # call function above
