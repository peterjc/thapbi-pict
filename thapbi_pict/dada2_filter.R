#!/usr/bin/env Rscript
# Copyright 2018, Peter Cock, James Hutton Institute.
#
# R script for use at the command line to use the R package DADA2
# to perform paired Illumina read sequencing error correction,
# based heavily on https://benjjneb.github.io/dada2/tutorial.html
# https://benjjneb.github.io/dada2/bigdata.html and
# https://benjjneb.github.io/dada2/bigdata_paired.html

# Read in package to read in command-line options.
library("optparse")

version <- "1.0" # The wrapper script version

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Directory of input paired FASTQ files.", metavar="path"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Directory for output FASTQ files.", metavar="path"),
  make_option(c("-t", "--threads"), type="integer", default=1,
              help="Number of threads to use (default: 1).", metavar="integer"),
  make_option(c("--verbose"), action="store_true", type="logical", default=FALSE,
              help="Write out status messages (default: FALSE).", metavar="boolean"),
  make_option(c("--version"), action="store_true", type="logical", default=FALSE,
              help="Print out version number and exit.", metavar="boolean")
)

opt_parser <- OptionParser(
  option_list=option_list, 
  usage = "%prog [options] -i INPUT_DIR -o OUTPUT_DIR",
  description = paste(
    "\n",
    "This is a wrapper script for the DADA2 filter step ",
    "based on the authors\' big data tutorial available here: ",
    "https://benjjneb.github.io/dada2/bigdata.html.\n\n",
    "Be sure to cite the DADA2 paper if you use this script:\n",
    "Callahan BJ et al. 2016. DADA2: ",
    "High-resolution sample inference from Illumina amplicon data. ",
    "Nature Methods 13:581-583.\n",
    "https://dx.doi.org/10.1038%2Fnmeth.3869\n\n",
    sep="")
)

opt <- parse_args(opt_parser)

# Load DADA2
# Seems to trigger annoying stderr text: "Loading required package: Rcpp"
library(dada2);

if (opt$verbose || opt$version) {
  cat("Wrapper version:", version, "\n")
  cat("DADA2 version:", as.character(packageVersion("dada2")), "\n")
}

# Printed out version, now quit if --version flag set.
if (opt$version) {
  quit()
}

# Argument checks
if (is.null(opt$input)) {
  stop("Required argument -i or --input missing.")
}

if (is.null(opt$output)) {
  stop("Required argument -o or --output missing.")
}

# If threads is zero or negative, set to false
if (opt$threads > 1) {
  multithread_opt <- opt$threads
} else {
  multithread_opt <- FALSE
}

# File parsing
path <- opt$input
filtpath <- opt$output
fastqFs <- sort(list.files(path, pattern="_R1_001.fastq")) # accepts *_R1_000.fastq.gz too! 
fastqRs <- sort(list.files(path, pattern="_R2_001.fastq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

if (0 == length(fastqFs)) {
  cat("WARNING: Didn't identify any FASTQ pairs in ", path, "\n")
  quit()
}

if (opt$verbose) {
  cat("Processing ", length(fastqFs), " paired FASTQ files\n")
}

cat("Fwd: ", fastqFs, "\n")
cat("Rev: ", fastqRs, "\n")

# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
filterAndTrim(fwd=file.path(path, fastqFs), filt=file.path(filtpath, fastqFs),
              rev=file.path(path, fastqRs), filt.rev=file.path(filtpath, fastqRs),
              truncLen=c(240,200), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=multithread_opt)

cat("Done, processed ", length(fastqFs), " paired FASTQ files\n")
