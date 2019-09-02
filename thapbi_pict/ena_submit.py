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

# Using https://www.ebi.ac.uk/ena/submit/drop-box/submit/
# Taxon unclassified sequences id 12908 is not submittable.
# Taxon metagenomes id 408169 is not submittable.
# However, e.g. 939928 rhizosphere metagenome works.

XML_HEADER = """<?xml version="1.0" encoding="UTF-8"?>
<SAMPLE_SET>
"""

XML_SAMPLE_HEADER = """  <SAMPLE alias="%s" center_name="">
"""

# XML_SAMPLE_NAME = """    <SAMPLE_NAME>
#       <TAXON_ID>1284369</TAXON_ID>
#       <SCIENTIFIC_NAME>stomach metagenome</SCIENTIFIC_NAME>
#       <COMMON_NAME></COMMON_NAME>
#     </SAMPLE_NAME>
# """

# From the test service, <SAMPLE_NAME> is required, and
# need at least one entry giving the species (e.g. taxid).
XML_SAMPLE_NAME = """    <SAMPLE_NAME>
      <TAXON_ID>%i</TAXON_ID>
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
    output=None,
    shared=None,
    mapping=None,
    metadata_file=None,
    metadata_cols=None,
    metadata_fieldnames=None,
    metadata_index=None,
    metadata_ncbi_taxid=None,
    default_ncbi_taxid=None,
    tmp_dir=None,
    debug=False,
):
    """Implement the ``thapbi_pict ena-submit`` command."""
    if default_ncbi_taxid:
        default_ncbi_taxid = int(default_ncbi_taxid)
    elif not metadata_ncbi_taxid:
        sys.exit("ERROR: Missing taxid configuration.")

    fastq_file_pairs = find_fastq_pairs(fastq, debug=debug)

    if debug:
        sys.stderr.write("Preparing %i FASTQ pairs\n" % len(fastq_file_pairs))

    if shared:
        if debug:
            sys.stderr.write("Loading shared sample metadata from %s\n" % shared)
        with open(shared) as handle:
            shared_attr = []
            for line in handle:
                if line.startswith("#"):
                    continue
                key, value = line.rstrip("\n").split("\t", 1)
                shared_attr.append(XML_ATTR % (key, value))
            shared_attr = "".join(shared_attr)
    else:
        shared_attr = None

    if tmp_dir:
        # Up to the user to remove the files
        tmp_obj = None
        shared_tmp = tmp_dir
    else:
        tmp_obj = tempfile.TemporaryDirectory()
        shared_tmp = tmp_obj.name

    if debug:
        sys.stderr.write("DEBUG: Temp folder %s\n" % shared_tmp)

    if output == "-":
        sample_xml_handle = sys.stdout
    elif os.path.isdir(output):
        sample_xml_handle = open(os.path.join(shared_tmp, "sample.xml"), "w")
    else:
        sys.exit("ERROR: Output directory does not exist: %s\n" % output)

    sample_xml_handle.write(XML_HEADER)

    added_taxid = False
    if metadata_ncbi_taxid:
        metadata_ncbi_taxid = int(metadata_ncbi_taxid)
        cols = [int(_) for _ in metadata_cols.split(",")]
        if metadata_ncbi_taxid not in [int(_) for _ in metadata_cols.split(",")]:
            # Add it to the columns list...
            added_taxid = True
            cols.append(metadata_ncbi_taxid)
            metadata_cols = ",".join(str(_) for _ in cols)
        taxid_offset = cols.index(metadata_ncbi_taxid)
        del cols

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
        meta = list(zip(meta_names, metadata.get(sample, meta_default)))
        if metadata_ncbi_taxid:
            taxid = meta[taxid_offset][1]
            if taxid:
                try:
                    taxid = int(taxid)
                except ValueError:
                    sys.exit(
                        "ERROR: Expected integer for NCBI TAXID for %s, got %s"
                        % (sample, taxid)
                    )
            else:
                taxid = default_ncbi_taxid
        else:
            taxid = default_ncbi_taxid
        if added_taxid:
            meta = meta[:-1]
        if debug:
            sys.stderr.write("DEBUG: %s\n" % stem)
        sample_xml_handle.write(XML_SAMPLE_HEADER % sample)
        sample_xml_handle.write(XML_SAMPLE_NAME % taxid)

        if shared or metadata_cols:
            sample_xml_handle.write(XML_ATTRS_HEADER)
            if shared:
                sample_xml_handle.write(shared_attr)
            for key, value in meta:
                if value:
                    # TODO: Apply mapping
                    sample_xml_handle.write(XML_ATTR % (key, value))
            sample_xml_handle.write(XML_ATTRS_FOOTER)
        sample_xml_handle.write(XML_SAMPLE_FOOTER)

    sample_xml_handle.write(XML_FOOTER)
    if output != "-":
        sample_xml_handle.close()
        shutil.move(
            os.path.join(shared_tmp, "sample.xml"), os.path.join(output, "sample.xml")
        )
