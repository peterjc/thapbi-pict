#!/usr/bin/env python
"""Python 3 script to make an md5sum check file from ENA table.

Expected usage, where a TSV file has been download from the ENA::

    $ python my5.py ENA_PROJECT_METADATA.tsv > raw_data/MD5SUM.txt
    $ cd raw_data/
    $ md5sum -c MD5SUM.txt
    $ cd ..

This requires at least the "fastq_md" and "fastq_ftp" columns be present.
"""
import os
import sys

if len(sys.argv) != 2 or not os.path.isfile(sys.argv[1]):
    sys.exit("ERROR: Expects a single argument, an ENA project TSV filename.")

checksum_col = None
url_col = None
with open(sys.argv[1]) as handle:
    for line in handle:
        parts = line.strip().split("\t")
        if "fastq_md5" in parts:
            checksum_col = parts.index("fastq_md5")
            url_col = parts.index("fastq_ftp")
            continue  # ignore header
        for md5, url in zip(parts[checksum_col].split(";"), parts[url_col].split(";")):
            filename = os.path.split(url)[1]
            if filename.startswith("EM"):
                filename = filename.replace(".fastq.gz", ".zip")
            print(f"{md5}  {filename}")
