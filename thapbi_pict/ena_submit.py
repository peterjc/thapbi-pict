# Copyright 2019-2020 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Code for sample submission to ENA/SRA.

This implements the ``thapbi_pict ena-submit ...`` command.
"""
import os
import shutil
import sys
import tempfile

from .prepare import find_fastq_pairs
from .utils import load_metadata


TABLE_HEADER = (
    "sample_alias\tinstrument_model\tlibrary_name\tlibrary_source\t"
    "library_selection\tlibrary_strategy\tdesign_description\t"
    "library_construction_protocol\tinsert_size\t"
    "forward_file_name\tforward_file_md5\t"
    "reverse_file_name\treverse_file_md5\n"
)
TABLE_TEMPLATE = "%s\t%s\t%s\tMETAGENOMIC\tPCR\tAMPLICON\t%s\t%s\t%i\t%s\t%s\t%s\t%s\n"
assert TABLE_HEADER.count("\t") == TABLE_TEMPLATE.count("\t")


def load_md5(file_list):
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
    pairs,
    meta,
    library_name,
    instrument_model,
    design_description,
    library_construction_protocol,
    insert_size,
):
    """Write read file table for ENA upload."""
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
                instrument_model,
                folder if library_name == "-" else library_name,
                design_description,
                library_construction_protocol,
                insert_size,
                os.path.split(raw_R1)[1],
                md5_dict[raw_R1],
                os.path.split(raw_R2)[1],
                md5_dict[raw_R2],
            )
        )
    handle.write(TABLE_HEADER)
    for line in sorted(lines):
        handle.write(line)


def main(
    fastq,
    output,
    metadata_file=None,
    metadata_cols=None,
    metadata_fieldnames=None,
    metadata_index=None,
    ignore_prefixes=None,
    library_name="-",
    instrument_model="Illumina MiSeq",
    design_description="",
    library_construction_protocol="",
    insert_size=250,
    tmp_dir=None,
    debug=False,
):
    """Implement the ``thapbi_pict ena-submit`` command."""
    fastq_file_pairs = find_fastq_pairs(fastq, debug=debug)

    if not fastq_file_pairs:
        sys.exit("ERROR: No FASTQ pairs found")
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

    (metadata, _, meta_names, group_col,) = load_metadata(
        metadata_file,
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
        sys.exit(
            "ERROR: Loaded %i samples, %i missing metadata, e.g. %s\n"
            % (len(samples), len(missing_meta), sorted(missing_meta)[0])
        )
    if metadata_file:
        # Just one column, don't need values as list:
        meta = {k: v[0] for k, v in metadata.items()}
    else:
        meta = None

    write_table(
        table_handle,
        fastq_file_pairs,
        meta,
        library_name,
        instrument_model,
        design_description,
        library_construction_protocol,
        insert_size,
    )

    if output != "-":
        table_handle.close()
        shutil.move(tmp_output, output)
