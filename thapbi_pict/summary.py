# Copyright 2020 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Summarise classification results at sample and read level.

This implements the ``thapbi_pict summary ...`` command.
"""
import os
import sys

from .read_summary import main as _read_summary
from .sample_summary import main as _sample_summary


def main(
    inputs,
    out_dir,
    report,
    method,
    min_abundance=1,
    metadata_file=None,
    metadata_cols=None,
    metadata_groups=None,
    metadata_fieldnames=None,
    metadata_index=None,
    require_metadata=False,
    ignore_prefixes=None,
    debug=False,
):
    """Implement the ``thapbi_pict summary`` command.

    The expectation is that the inputs represent all the samples from
    a meaningful group, likely from multiple sequencing runs (plates).
    """
    # TODO - refactor the old separate reporting code
    assert isinstance(inputs, list)

    if report:
        stem = os.path.join(out_dir, report)
    else:
        # Include version number here?
        stem = os.path.join(out_dir, "thapbi-pict")

    return_code = _sample_summary(
        inputs=[_ for _ in inputs if not _.endswith(".fasta")],  # TSV only
        output=f"{stem}.samples.{method}.tsv",
        excel=f"{stem}.samples.{method}.xlsx",
        human_output=f"{stem}.samples.{method}.txt",
        method=method,
        min_abundance=min_abundance,
        metadata_file=metadata_file,
        metadata_cols=metadata_cols,
        metadata_groups=metadata_groups,
        metadata_fieldnames=metadata_fieldnames,
        metadata_index=metadata_index,
        require_metadata=require_metadata,
        ignore_prefixes=ignore_prefixes,
        debug=debug,
    )
    if return_code:
        return return_code
    sys.stderr.write(f"Wrote {stem}.samples.{method}.*\n")

    return_code = _read_summary(
        inputs=inputs,  # needs FASTA and TSV
        output=f"{stem}.reads.{method}.tsv",
        excel=f"{stem}.reads.{method}.xlsx",
        method=method,
        min_abundance=min_abundance,
        metadata_file=metadata_file,
        metadata_cols=metadata_cols,
        metadata_groups=metadata_groups,
        metadata_fieldnames=metadata_fieldnames,
        metadata_index=metadata_index,
        require_metadata=require_metadata,
        ignore_prefixes=ignore_prefixes,
        debug=debug,
    )
    if return_code:
        return return_code
    sys.stderr.write(f"Wrote {stem}.reads.{method}.*\n")

    return return_code
