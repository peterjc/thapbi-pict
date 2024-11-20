# Copyright 2019-2024 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Code for sample submission to ENA/SRA.

This implements the ``thapbi_pict ena-submit ...`` command.
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile

from .prepare import find_fastq_pairs
from .utils import load_metadata

# As of Nov 2024 and earlier, there are 12 mandatory fields:
#
# Field Name         Field Label            Permitted Value
# sample             Sample
# study              Study
# instrument_model   Instrument model       Illumina MiSeq
# library_name       Library name
# library_source     Library source         METAGENOMIC
# library_selection  Library selection      PCR
# library_strategy   Library strategy       AMPLICON
# library_layout     Library layout         PAIRED
# forward_file_name  Forward file name
# forward_file_md5   Forward File checksum
# reverse_file_name  Reverse file name
# reverse_file_md5   Reverse file checksum
#
# We will ignore the optional fields (some of which used to be expected
# but accepted blank values):
#
# Field Name                    Field Label
# library_design                Library design
# library_construction_protocol Library construction protocol
# design_description            Design description
# insert_size                   Insert size
# forward_file_unencrypted_md5  Forward file unencrypted checksum
# reverse_file_unencrypted_md5  Reverse file unencrypted checksum

TABLE_HEADER = (
    "sample\tstudy\tinstrument_model\tlibrary_name\tlibrary_source\t"
    "library_selection\tlibrary_strategy\tlibrary_layout\t"
    "forward_file_name\tforward_file_md5\t"
    "reverse_file_name\treverse_file_md5\n"
)
TABLE_TEMPLATE = "%s\t%s\t%s\t%s\tMETAGENOMIC\tPCR\tAMPLICON\tPAIRED\t%s\t%s\t%s\t%s\n"
assert TABLE_HEADER.count("\t") == TABLE_TEMPLATE.count("\t")


def load_md5(file_list: list[str]) -> dict[str, str]:
    """Return a dict mapping given filenames to MD5 digests."""
    assert file_list, "Nothing to do here."
    answer = {}
    base_names = {os.path.split(_)[1] for _ in file_list}
    if len(base_names) < len(file_list):
        # This isn't a problem locally, but will be on upload to ENA
        sys.exit("ERROR: Duplicate FASTQ names once folder dropped")
    checksum_files = {_ + ".md5" for _ in file_list}
    checksum_files.update(
        os.path.join(os.path.split(_)[0], "MD5SUM.txt") for _ in file_list
    )
    for cache in checksum_files:
        if os.path.isfile(cache):
            with open(cache) as handle:
                for line in handle:
                    md5, filename = line.strip().split()
                    # If MD5 file contains relative path, interpret it
                    assert "/" not in filename, filename
                    filename = os.path.join(os.path.split(cache)[0], filename)
                    if filename in file_list:
                        answer[filename] = md5
    if not answer:
        sys.exit("ERROR: Need to provide MD5SUM.txt or example.fastq.gz.md5 files")
    for f in file_list:
        if f not in answer:
            sys.exit(f"ERROR: Need MD5 for {f} and not in {f}.md5 or MD5SUM.txt")
    return answer


def write_table(
    handle,
    pairs: list[tuple[str, str, str]],
    meta: dict[str, str] | None,
    study: str,
    library_name: str,
    instrument_model: str,
    flat: bool = False,
) -> None:
    """Write read file table for ENA upload.

    Leave library name as "-" and the folder name will be used.
    """
    file_list = [_[1] for _ in pairs] + [_[2] for _ in pairs]
    md5_dict = load_md5(file_list)
    lines = []
    for stem, raw_R1, raw_R2 in pairs:
        sample = os.path.split(stem)[1]
        folder = os.path.split(os.path.split(raw_R1)[0])[1]
        lines.append(
            TABLE_TEMPLATE
            % (
                meta[sample] if meta else sample,
                study,
                instrument_model,
                folder if library_name == "-" else library_name,
                os.path.split(raw_R1)[1] if flat else raw_R1.replace("\\", "/"),
                md5_dict[raw_R1],
                os.path.split(raw_R2)[1] if flat else raw_R2.replace("\\", "/"),
                md5_dict[raw_R2],
            )
        )
    handle.write(TABLE_HEADER)
    for line in sorted(lines):
        handle.write(line)


def main(
    fastq: list[str],
    output: str,
    study: str,
    metadata_file: str | None = None,
    metadata_encoding: str | None = None,
    metadata_cols: str | None = None,  # single column with sample name
    metadata_fieldnames: str | None = None,
    metadata_index: str | None = None,
    ignore_prefixes: str | None = None,
    library_name: str = "-",
    instrument_model: str = "Illumina MiSeq",
    ignore_stems: str | None = None,
    tmp_dir: str | None = None,
    flat: bool = False,
    debug: bool = False,
):
    """Implement the ``thapbi_pict ena-submit`` command."""
    fastq_file_pairs = find_fastq_pairs(fastq, debug=debug)
    if not fastq_file_pairs:
        sys.exit("ERROR: No FASTQ pairs found")

    if ignore_stems:
        ignore = set()
        with open(ignore_stems) as handle:
            for line in handle:
                if line and not line.startswith("#"):
                    ignore.add(line.strip())
        before = len(fastq_file_pairs)
        fastq_file_pairs = [
            (stem, left_filename, right_filename)
            for (stem, left_filename, right_filename) in fastq_file_pairs
            if os.path.basename(stem) not in ignore
        ]
        sys.stderr.write(
            f"Ignored {before-len(fastq_file_pairs)} paired FASTQ stems using"
            f" {len(ignore)} potential entries from file {ignore_stems}\n"
        )
        del ignore
        if not fastq_file_pairs:
            sys.exit("ERROR: Ignored all FASTQ pairs found")

    if debug:
        sys.stderr.write("Preparing %i FASTQ pairs\n" % len(fastq_file_pairs))

    if tmp_dir:
        # Up to the user to remove the files
        tmp_obj = None
        shared_tmp = tmp_dir
    else:
        tmp_obj = tempfile.TemporaryDirectory()
        shared_tmp = tmp_obj.name

    if debug:
        sys.stderr.write("DEBUG: Temp folder %s\n" % shared_tmp)

    tmp_output = None
    if output == "-":
        table_handle = sys.stdout
    elif os.path.isdir(output):
        sys.exit(f"ERROR: Output filename {output} is a directory")
    else:
        tmp_output = os.path.join(shared_tmp, "experiment_paired_fastq_spreadsheet.tsv")
        table_handle = open(tmp_output, "w")

    samples = sorted(
        os.path.basename(stem) for stem, _raw_R1, _raw_R2 in fastq_file_pairs
    )

    (
        metadata,
        _,
        meta_names,
        group_col,
    ) = load_metadata(
        metadata_file,
        metadata_encoding,
        metadata_cols,
        None,  # i.e. metadata_groups=None,
        metadata_fieldnames,
        metadata_index,
        ignore_prefixes=ignore_prefixes,
        debug=debug,
    )
    if metadata_file and len(meta_names) != 1:
        sys.exit(
            "Need one and only one column for the -c argument, "
            "giving the ENA sample accession or refname."
        )
    missing_meta = set(samples).difference(metadata)
    if metadata_file and missing_meta:
        if len(missing_meta) == 1:
            msg = "ERROR: Loaded %i samples, %i missing metadata, %s\n" % (
                len(samples),
                len(missing_meta),
                sorted(missing_meta)[0],
            )
        else:
            msg = "ERROR: Loaded %i samples, %i missing metadata, e.g. %s .. %s\n" % (
                len(samples),
                len(missing_meta),
                sorted(missing_meta)[0],
                sorted(missing_meta)[-1],
            )
        sys.exit(msg)
    if metadata_file:
        # Just one column, don't need values as list:
        meta = {k: v[0] for k, v in metadata.items()}
    else:
        meta = None

    write_table(
        table_handle,
        fastq_file_pairs,
        meta,
        study,
        library_name,
        instrument_model,
        flat=flat,
    )

    if output != "-":
        assert table_handle is not None
        table_handle.close()
        assert tmp_output is not None
        shutil.move(tmp_output, output)
