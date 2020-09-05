#!/usr/bin/env python
"""Convert PRJNA417859.txt into metadata.tsv for THAPBI PICT."""
import os

months = {
    "Jan": 1,
    "Feb": 2,
    "Mar": 3,
    "Apr": 4,
    "April": 4,
    "May": 5,
    "Jun": 6,
    "Jul": 7,
    "Aug": 8,
    "Sep": 9,
    "Oct": 10,
    "Nov": 11,
    "Dec": 12,
}

# Columns:
#
# 0 run_accession
# 1 experiment_title
# 2 experiment_alias
# 3 run_alias
# 4 fastq_md5
# 5 fastq_ftp
# 6 sample_title

with open("PRJNA417859.txt") as handle:
    with open("metadata.tsv", "w") as meta:
        # Write our header line
        meta.write("Source\tSite\tProcess\tPeriod\tYear-Month\tSample\tAccession\n")
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            acc = parts[0]
            if acc == "run_accession":
                # Ignore the header line
                continue
            name = parts[2]

            title = parts[1]
            processed = title.split("Processed by ", 1)[1]
            source, site, period, month_year = parts[-1].split(", ")

            assert site.startswith("Site: ")
            site = site[6:]
            assert source == "Irrigation water"
            # e.g. River-A
            # The letter codes are globally unique
            source, site = site.split("-", 1)

            assert period.startswith("Period: ")
            # Add leading zero for sorting if treated as text
            period = "%02i" % int(period[8:].strip())

            assert month_year.startswith("Month-Year: ")
            month_year = month_year[12:]
            month, year = month_year.split("-")
            try:
                yyyy_mm = f"{int(year):04d}-{months[month]:02d}"
            except KeyError:
                print(month_year, month, year)
                raise

            meta.write(
                f"{source}\t{site}\t{processed}\t{period}\t{yyyy_mm}\t{name}\t{acc}\n"
            )

            R1 = f"raw_data/{acc}_1.fastq.gz"
            R2 = f"raw_data/{acc}_2.fastq.gz"

            if not os.path.isfile(R1) or not os.path.isfile(R2):
                print(f"Missing file(s) for {acc}")
print("Wrote metadata TSV file")
