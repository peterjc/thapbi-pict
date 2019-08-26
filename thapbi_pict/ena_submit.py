# Copyright 2019 by Peter Cock, The James Hutton Institute.
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
from .utils import load_metadata_dict
from .utils import sample_sort


XML_HEADER = """<?xml version="1.0" encoding="UTF-8"?>
<SAMPLE_SET>
"""

XML_SAMPLE_HEADER = """  <SAMPLE alias="%s" center_name="">
"""

XML_SAMPLE_NAME = """    <SAMPLE_NAME>
      <TAXON_ID>1284369</TAXON_ID>
      <SCIENTIFIC_NAME>stomach metagenome</SCIENTIFIC_NAME>
      <COMMON_NAME></COMMON_NAME>
    </SAMPLE_NAME>
"""

XML_ATTRS_HEADER = """    <SAMPLE_ATTRIBUTES>
"""

XML_ATTR = """      <SAMPLE_ATTRIBUTE>
        <TAG>%s</TAG>
        <VALUE>%s</VALUE>
      </SAMPLE_ATTRIBUTE>
"""

XML_ATTRS_FOOTER = """    </SAMPLE_ATTRIBUTES>
"""

XML_SAMPLE_FOOTER = """  </SAMPLE>
"""

XML_FOOTER = """</SAMPLE_SET>
"""


def main(
    fastq,
    samplexml=None,
    mapping=None,
    metadata_file=None,
    metadata_cols=None,
    metadata_fieldnames=None,
    metadata_index=None,
    tmp_dir=None,
    debug=False,
):
    """Implement the ``thapbi_pict ena-submit`` command."""
    fastq_file_pairs = find_fastq_pairs(fastq, debug=debug)

    if debug:
        sys.stderr.write("Preparing %i FASTQ pairs\n" % len(fastq_file_pairs))

    combined_xml = samplexml and (samplexml == "-" or not os.path.isdir(samplexml))

    if tmp_dir:
        # Up to the user to remove the files
        tmp_obj = None
        shared_tmp = tmp_dir
    else:
        tmp_obj = tempfile.TemporaryDirectory()
        shared_tmp = tmp_obj.name

    if debug:
        sys.stderr.write("DEBUG: Temp folder %s\n" % shared_tmp)

    if combined_xml:
        if samplexml == "-":
            sample_xml_handle = sys.stdout
        else:
            sample_xml_handle = open(os.path.join(shared_tmp, "combined.xml"), "w")
        sample_xml_handle.write(XML_HEADER)
    else:
        sample_xml_handle = None

    samples = sample_sort(
        os.path.basename(stem) for stem, _raw_R1, _raw_R2 in fastq_file_pairs
    )
    (metadata, meta_names, meta_default, missing_meta) = load_metadata_dict(
        samples,
        metadata_file,
        metadata_cols,
        metadata_fieldnames,
        metadata_index,
        sequenced_samples=samples,
        metadata_sort=True,
        debug=debug,
    )
    if debug:
        sys.stderr.write(
            "Loaded %i samples, %i missing metadata\n"
            % (len(samples), len(missing_meta))
        )

    for stem, _raw_R1, _raw_R2 in fastq_file_pairs:
        sample = os.path.split(stem)[1]
        if not combined_xml:
            sample_xml_handle = open(os.path.join(shared_tmp, sample + ".xml"), "w")
            sample_xml_handle.write(XML_HEADER)
        if debug:
            sys.stderr.write("DEBUG: %s\n" % stem)
        sample_xml_handle.write(XML_SAMPLE_HEADER % stem)
        if metadata_cols:
            sample_xml_handle.write(XML_ATTRS_HEADER)
            for key, value in zip(meta_names, metadata.get(sample, meta_default)):
                if value:
                    # TODO: Apply mapping
                    sample_xml_handle.write(XML_ATTR % (key, value))
            sample_xml_handle.write(XML_ATTRS_FOOTER)
        sample_xml_handle.write(XML_SAMPLE_FOOTER)
        if not combined_xml:
            sample_xml_handle.write(XML_FOOTER)
            sample_xml_handle.close()
            shutil.move(
                os.path.join(shared_tmp, sample + ".xml"),
                os.path.join(
                    samplexml if samplexml else os.path.split(stem)[0], sample + ".xml"
                ),
            )

    if combined_xml:
        sample_xml_handle.write(XML_FOOTER)
        if samplexml != "-":
            sample_xml_handle.close()
            shutil.move(os.path.join(shared_tmp, "combined.xml"), samplexml)
