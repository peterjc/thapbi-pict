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
from .utils import load_metadata
from .utils import sample_sort

# Using https://www.ebi.ac.uk/ena/submit/drop-box/submit/
# Taxon unclassified sequences id 12908 is not submittable.
# Taxon metagenomes id 408169 is not submittable.
# However, e.g. 939928 rhizosphere metagenome works.

XML_SAMPLE_SET_HEADER = """<?xml version="1.0" encoding="UTF-8"?>
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

XML_SAMPLE_SET_FOOTER = """</SAMPLE_SET>
"""

######################

XML_EXPR_SET_HEADER = """<EXPERIMENT_SET>
"""

XML_EXPR_TEMPLATE = """   <EXPERIMENT alias="exp_%s">
       <TITLE>%s</TITLE>
       <STUDY_REF accession="%s"/>
       <DESIGN>
           <DESIGN_DESCRIPTION/>
           <SAMPLE_DESCRIPTOR refname="%s"/>
           <LIBRARY_DESCRIPTOR>
               <LIBRARY_NAME/>
               <LIBRARY_STRATEGY>AMPLICON</LIBRARY_STRATEGY>
               <LIBRARY_SOURCE>METAGENOMIC</LIBRARY_SOURCE>
               <LIBRARY_SELECTION>PCR</LIBRARY_SELECTION>
               <LIBRARY_LAYOUT>
                   <PAIRED NOMINAL_LENGTH="250" />
               </LIBRARY_LAYOUT>
               <LIBRARY_CONSTRUCTION_PROTOCOL>Nested PCR with amplicon specific primers.</LIBRARY_CONSTRUCTION_PROTOCOL>
           </LIBRARY_DESCRIPTOR>
       </DESIGN>
       <PLATFORM>
           <ILLUMINA>
               <INSTRUMENT_MODEL>Illumina MiSeq</INSTRUMENT_MODEL>
           </ILLUMINA>
       </PLATFORM>
       <EXPERIMENT_ATTRIBUTES>
           <EXPERIMENT_ATTRIBUTE>
               <TAG>library preparation date</TAG>
               <VALUE>%s</VALUE>
           </EXPERIMENT_ATTRIBUTE>
       </EXPERIMENT_ATTRIBUTES>
   </EXPERIMENT>
"""  # noqa: E501

XML_EXPR_SET_FOOTER = """</EXPERIMENT_SET>
"""

######################

XML_RUN_SET_HEADER = """<RUN_SET>
"""

XML_RUN_TEMPLATE = """    <RUN alias="run_%s" center_name="">
        <EXPERIMENT_REF refname="exp_%s"/>
        <DATA_BLOCK>
            <FILES>
                <FILE filename="%s" filetype="fastq"
                    checksum_method="MD5" checksum="%s"/>
                <FILE filename="%s" filetype="fastq"
                    checksum_method="MD5" checksum="%s"/>
            </FILES>
        </DATA_BLOCK>
    </RUN>
"""

XML_RUN_SET_FOOTER = """</RUN_SET>
"""

######################


def write_expr(handle, pairs):
    """Write experiment.xml to handle."""
    handle.write(XML_EXPR_SET_HEADER)
    for stem, _raw_R1, _raw_R2 in pairs:
        sample = os.path.split(stem)[1]
        title = ""
        study_accession = ""
        seq_date = ""  # get from directory name?
        handle.write(
            XML_EXPR_TEMPLATE % (sample, title, study_accession, sample, seq_date)
        )
    handle.write(XML_EXPR_SET_FOOTER)


def load_md5(file_list):
    """Return a dict mapping given filenames to MD5 digests."""
    answer = {}
    for cache in {os.path.join(os.path.split(_)[0], "MD5SUM.txt") for _ in file_list}:
        if os.path.isfile(cache):
            with open(cache) as handle:
                for line in handle:
                    md5, filename = line.strip().split()
                    filename = os.path.join(os.path.split(cache)[0], filename)
                    if filename in file_list:
                        answer[filename] = md5
    for f in file_list:
        if f not in answer:
            # TODO - do this at run time? Too slow?
            sys.exit(f"ERROR: Need MD5 for {f} and not in MD5SUM.txt")
    return answer


def write_run(handle, pairs):
    """Write run.xml to handle."""
    file_list = [_[1] for _ in pairs] + [_[2] for _ in pairs]
    md5_dict = load_md5(file_list)
    handle.write(XML_RUN_SET_HEADER)
    for stem, raw_R1, raw_R2 in pairs:
        sample = os.path.split(stem)[1]
        handle.write(
            XML_RUN_TEMPLATE
            % (sample, sample, raw_R1, md5_dict[raw_R1], raw_R2, md5_dict[raw_R2])
        )
    handle.write(XML_RUN_SET_FOOTER)


######################


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
    ignore_prefixes=None,
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
        expr_xml_handle = sys.stdout
        run_xml_handle = sys.stdout
    elif os.path.isdir(output):
        sample_xml_handle = open(os.path.join(shared_tmp, "sample.xml"), "w")
        expr_xml_handle = open(os.path.join(shared_tmp, "experiment.xml"), "w")
        run_xml_handle = open(os.path.join(shared_tmp, "run.xml"), "w")
    else:
        sys.exit("ERROR: Output directory does not exist: %s\n" % output)

    sample_xml_handle.write(XML_SAMPLE_SET_HEADER)
    # expr_xml_handle.write(XML_EXPR_SET_HEADER)

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

    # (metadata, meta_names, meta_default, missing_meta) = load_metadata_dict(...)
    (metadata, meta_names, group_col,) = load_metadata(
        metadata_file,
        metadata_cols,
        None,  # i.e. metadata_groups=None,
        metadata_fieldnames,
        metadata_index,
        metadata_sort=True,
        ignore_prefixes=ignore_prefixes,
        debug=debug,
    )
    missing_meta = set(samples).difference(metadata)
    if debug:
        sys.stderr.write(
            "Loaded %i samples, %i missing metadata\n"
            % (len(samples), len(missing_meta))
        )

    for stem, _raw_R1, _raw_R2 in fastq_file_pairs:
        sample = os.path.split(stem)[1]
        meta = metadata[sample]
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
            for key, value in zip(meta_names, meta):
                if value:
                    # TODO: Apply mapping
                    sample_xml_handle.write(XML_ATTR % (key, value))
            sample_xml_handle.write(XML_ATTRS_FOOTER)
        sample_xml_handle.write(XML_SAMPLE_FOOTER)

    sample_xml_handle.write(XML_SAMPLE_SET_FOOTER)

    write_expr(expr_xml_handle, fastq_file_pairs)

    write_run(run_xml_handle, fastq_file_pairs)

    if output != "-":
        sample_xml_handle.close()
        expr_xml_handle.close()
        run_xml_handle.close()
        for name in ("sample.xml", "experiment.xml", "run.xml"):
            shutil.move(os.path.join(shared_tmp, name), os.path.join(output, name))
