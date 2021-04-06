# Copyright 2018-2021 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Defines the ``thapbi_pict`` command line tool.

This works via ``setup.py`` where under ``entry_points`` we define a
``console_scripts`` entry for ``thapbi_pict`` (executable name) pointing to
the ``main()`` function define in this Python file.
"""
import argparse
import os
import sys

from . import __version__
from .classify import method_classify_file as method_classifier
from .utils import find_requested_files
from .utils import primer_clean


# Common command line defaults
# ============================


DEFAULT_METHOD = "onebp"
DEFAULT_MIN_ABUNDANCE = 100
CTRL_A = chr(1)
IGNORE_PREFIXES = ("Undetermined",)

# Argument validation functions
# =============================


def check_input_file(filename):
    """Command line validation of an input filename."""
    if not os.path.isfile(filename):
        sys.exit(f"ERROR: Could not find input file: {filename}")


def check_output_directory(out_dir, must_exist=True):
    """Command line validation of output directory value."""
    if not out_dir:
        sys.exit("ERROR: Output directory name blank\n")
    elif out_dir == "-" or os.path.isdir(out_dir):
        return True
    elif os.path.isfile(out_dir):
        sys.exit(f"ERROR: Output directory name is a file: {out_dir}\n")
    elif must_exist:
        sys.exit(f"ERROR: Output directory does not exist: {out_dir}\n")


def expand_database_argument(text, exist=False, hyphen_default=False):
    """Expand an SQLite3 filename to an SQLalchemy URL."""
    # TODO: Expand this to allow other DB prefixes later
    # Note we are not currently checking file exists,
    # as we might be about to create it.
    if text == "-":
        if hyphen_default:
            # Expand to the default bundled DB
            text = os.path.join(os.path.split(__file__)[0], "ITS1_DB.sqlite")
        else:
            sys.exit("ERROR: Using hyphen as DB default is not supported here.")
    if not text:
        sys.exit("ERROR: The database argument is required.")
    prefix = "sqlite:///"
    if text.startswith(prefix):
        db = text[len(prefix) :]
        assert text == prefix + db
    else:
        db = text
    if exist and db != ":memory:" and not os.path.isfile(db):
        sys.exit(f"ERROR: The database {db} was not found.\n")
    return prefix + db


def expand_hmm_argument(text, hyphen_default=True):
    """Validate/expand an HMMER3 HMM stem."""
    if text == "-":
        if hyphen_default:
            # Expand to the default bundled HMM
            text = os.path.join(os.path.split(__file__)[0], "hmm", "controls.hmm")
        else:
            sys.exit("ERROR: Using hyphen as HMM default is not supported here.")
    if not text:
        sys.stderr.write("WARNING: Applying no HMM filtering!\n")
        return ""
    for ext in (".h3f", ".h3i", ".h3m", ".h3p"):
        if not os.path.isfile(text + ext):
            sys.exit(f"ERROR: The HMM {text} was not found (e.g. {text + ext})")
    return text


# Subcommand dispatch
# ===================


def load_tax(args=None):
    """Subcommand to load an NCBI taxonomy dump into a basebase."""
    from .taxdump import main

    return main(
        tax=args.tax,
        db_url=expand_database_argument(args.database),
        ancestors=args.ancestors,
        debug=args.verbose,
    )


def ncbi_import(args=None):
    """Subcommand to import an NCBI FASTA file into a basebase."""
    from .ncbi import main

    return main(
        fasta_file=args.input,
        db_url=expand_database_argument(args.database),
        min_length=args.minlen,
        max_length=args.maxlen,
        name=args.name,
        validate_species=not args.lax,
        genus_only=args.genus,
        left_primer=primer_clean(args.left),
        right_primer=primer_clean(args.right),
        sep=args.sep,
        tmp_dir=args.temp,
        debug=args.verbose,
    )


def curated_seq(args=None):
    """Subcommand to prepare classified sequences for later curated import."""
    from .seq_import import main

    if args.output:
        check_output_directory(args.output, must_exist=False)
    return main(
        inputs=args.input,
        out_dir=args.output,
        method=args.method,
        min_abundance=args.abundance,
        ignore_prefixes=tuple(args.ignore_prefixes),
        tmp_dir=args.temp,
        debug=args.verbose,
    )


def curated_import(args=None):
    """Subcommand to import a curated marker FASTA file into the basebase."""
    from .curated import main

    return main(
        fasta=args.input,
        db_url=expand_database_argument(args.database),
        min_length=args.minlen,
        max_length=args.maxlen,
        name=args.name,
        validate_species=not args.lax,
        genus_only=args.genus,
        left_primer=primer_clean(args.left),
        right_primer=primer_clean(args.right),
        sep=args.sep,
        ignore_prefixes=tuple(args.ignore_prefixes),
        debug=args.verbose,
    )


def dump(args=None):
    """Subcommand to dump a database to a text file."""
    from .dump import main

    return main(
        db_url=expand_database_argument(args.database, exist=True, hyphen_default=True),
        output_filename=args.output,
        output_format=args.format,
        minimal=args.minimal,
        genus=args.genus,
        species=args.species,
        sep=args.sep,
        debug=args.verbose,
    )


def conflicts(args=None):
    """Subcommand to look for taxonomic conflicts in a database."""
    from .conflicts import main

    return main(
        db_url=expand_database_argument(args.database, exist=True, hyphen_default=True),
        output_filename=args.output,
        debug=args.verbose,
    )


def prepare_reads(args=None):
    """Subcommand to prepare FASTA files from paired FASTQ reads."""
    from .prepare import main

    check_output_directory(args.output, must_exist=False)
    if args.temp:
        check_output_directory(args.temp)
    return_code = main(
        fastq=args.input,
        negative_controls=args.negctrls,
        out_dir=args.output,
        hmm_stem=expand_hmm_argument(args.hmm),
        primer_dir=args.primers,
        left_primer=primer_clean(args.left),
        right_primer=primer_clean(args.right),
        flip=args.flip,
        min_abundance=args.abundance,
        min_length=args.minlen,
        max_length=args.maxlen,
        ignore_prefixes=tuple(args.ignore_prefixes),
        merged_cache=args.merged_cache,
        tmp_dir=args.temp,
        debug=args.verbose,
        cpu=args.cpu,
    )
    if isinstance(return_code, list):
        # Should be a list of FASTA filenames
        return 0
    else:
        return return_code


def fasta_nr(args=None):
    """Subcommand to prepare non-redundant MD5 named FASTA files."""
    from .fasta_nr import main

    if not (args.input or args.revcomp):
        sys.exit("ERROR: Require at least one of -i/--input or -r/--revcomp")

    return main(
        inputs=args.input,
        revcomp=args.revcomp,
        output=args.output,
        min_abundance=args.abundance,
        min_length=args.minlen,
        max_length=args.maxlen,
        debug=args.verbose,
    )


def classify(args=None):
    """Subcommand to classify FASTA sequences using a database."""
    from .classify import main

    if args.output:
        check_output_directory(args.output, must_exist=False)
    if args.temp:
        check_output_directory(args.temp)
    return_code = main(
        fasta=args.input,
        db_url=expand_database_argument(args.database, exist=True, hyphen_default=True),
        method=args.method,
        out_dir=args.output,
        ignore_prefixes=tuple(args.ignore_prefixes),
        min_abundance=args.abundance,
        tmp_dir=args.temp,
        debug=args.verbose,
        cpu=args.cpu,
    )
    if isinstance(return_code, list):
        # Should be a list of *.method.tsv filenames.
        return 0
    else:
        return return_code


def assess_classification(args=None):
    """Subcommand to assess classification using known true taxonomy."""
    from .assess import main

    return main(
        inputs=args.input,
        known=args.known,
        db_url=expand_database_argument(args.database, exist=True, hyphen_default=True),
        method=args.method,
        min_abundance=args.abundance,
        assess_output=args.output,
        map_output=args.table,
        confusion_output=args.confusion,
        ignore_prefixes=tuple(args.ignore_prefixes),
        debug=args.verbose,
    )


def summary(args=None):
    """Subcommand to run sample and read summary reports."""
    from .summary import main

    check_output_directory(args.output)

    if args.metadata:
        check_input_file(args.metadata)
        if not args.metacols:
            sys.exit("ERROR: Must also supply -c / --metacols argument.")

    return main(
        inputs=args.input,
        out_dir=args.output,
        report=args.report,
        method=args.method,
        min_abundance=args.abundance,
        metadata_file=args.metadata,
        metadata_cols=args.metacols,
        metadata_groups=args.metagroups,
        metadata_fieldnames=args.metafields,
        metadata_index=args.metaindex,
        require_metadata=args.requiremeta,
        show_unsequenced=args.unsequenced,
        ignore_prefixes=tuple(args.ignore_prefixes),
        debug=args.verbose,
    )


def edit_graph(args=None):
    """Subcommand to create sequence-level edit-distance graph."""
    from .edit_graph import main

    db = (
        None
        if not args.database
        else expand_database_argument(args.database, exist=True, hyphen_default=True)
    )

    return main(
        graph_output=args.output,
        graph_format=args.format,
        db_url=db,
        inputs=args.input,
        method=args.method,
        min_abundance=args.abundance,
        total_min_abundance=args.total,
        always_show_db=args.showdb,
        max_edit_dist=args.editdist,
        ignore_prefixes=tuple(args.ignore_prefixes),
        debug=args.verbose,
    )


def pipeline(args=None):
    """Subcommand to run the default classification pipeline."""
    from .prepare import main as prepare
    from .fasta_nr import main as fasta_nr
    from .classify import main as classify
    from .summary import main as summary
    from .assess import main as assess
    from .edit_graph import main as edit_graph

    check_output_directory(args.output)
    if args.sampleout:
        check_output_directory(args.sampleout)
        intermediate_dir = args.sampleout
    else:
        intermediate_dir = args.output
    if args.temp:
        check_output_directory(args.temp)
    if args.metadata:
        check_input_file(args.metadata)
        if not args.metacols:
            sys.exit("ERROR: Must also supply -c / --metacols argument.")
    db = expand_database_argument(args.database, exist=True, hyphen_default=True)

    hmm = expand_hmm_argument(args.hmm)

    if args.report:
        stem = os.path.join(args.output, args.report)
    else:
        # Include version number here?
        stem = os.path.join(args.output, "thapbi-pict")

    # TODO - apply require_metadata=True to the prepare and classify steps?

    fasta_files = prepare(
        fastq=args.input,
        negative_controls=args.negctrls,
        out_dir=intermediate_dir,
        hmm_stem=hmm,
        primer_dir=None,
        left_primer=primer_clean(args.left),
        right_primer=primer_clean(args.right),
        flip=args.flip,
        min_abundance=args.abundance,
        min_length=args.minlen,
        max_length=args.maxlen,
        ignore_prefixes=tuple(args.ignore_prefixes),
        merged_cache=args.merged_cache,
        tmp_dir=args.temp,
        debug=args.verbose,
        cpu=args.cpu,
    )
    if isinstance(fasta_files, int):
        return_code = fasta_files
        if return_code:
            sys.stderr.write("ERROR: Pipeline aborted during prepare-reads\n")
            sys.exit(return_code)

    all_fasta = stem + ".all_reads.fasta"
    fasta_nr(
        inputs=fasta_files,
        revcomp=None,
        output=all_fasta,
        min_abundance=args.abundance,
        min_length=args.minlen,
        max_length=args.maxlen,
        debug=args.verbose,
    )

    classified_files = classify(
        fasta=fasta_files,
        db_url=db,
        method=args.method,
        out_dir=intermediate_dir,
        ignore_prefixes=tuple(args.ignore_prefixes),  # not really needed
        min_abundance=args.abundance,
        tmp_dir=args.temp,
        debug=args.verbose,
        cpu=args.cpu,
    )
    if isinstance(classified_files, int):
        return_code = classified_files
        if return_code:
            sys.stderr.write("ERROR: Pipeline aborted during classify\n")
            sys.exit(return_code)
    if len(fasta_files) != len(classified_files):
        sys.exit(
            f"ERROR: {len(fasta_files)} FASTA files "
            f"but {len(classified_files)} classified"
        )

    method = args.method

    return_code = summary(
        inputs=fasta_files + classified_files,
        out_dir=args.output,
        report=args.report,
        method=args.method,
        min_abundance=args.abundance,
        metadata_file=args.metadata,
        metadata_cols=args.metacols,
        metadata_groups=args.metagroups,
        metadata_fieldnames=args.metafields,
        metadata_index=args.metaindex,
        require_metadata=args.requiremeta,
        show_unsequenced=args.unsequenced,
        ignore_prefixes=tuple(args.ignore_prefixes),
        debug=args.verbose,
    )
    if return_code:
        sys.stderr.write("ERROR: Pipeline aborted during summary\n")
        sys.exit(return_code)

    # TODO - Support known setting...
    known_files = find_requested_files(
        args.input,
        ".known.tsv",
        ignore_prefixes=tuple(args.ignore_prefixes),
        debug=args.verbose,
    )
    if known_files:
        return_code = assess(
            inputs=known_files + classified_files,
            known="known",  # =args.known,
            db_url=db,
            method=args.method,
            min_abundance=args.abundance,
            assess_output=f"{stem}.assess.{method}.tsv",
            map_output=f"{stem}.assess.tally.{method}.tsv",
            confusion_output=f"{stem}.assess.confusion.{method}.tsv",
            ignore_prefixes=tuple(args.ignore_prefixes),
            debug=args.verbose,
        )
        if return_code:
            sys.stderr.write("ERROR: Pipeline aborted during assess\n")
            sys.exit(return_code)
        sys.stderr.write(f"Wrote {stem}.assess*.{method}.*\n")

    edit_graph_filename = f"{stem}.edit-graph.{method}.xgmml"
    if os.path.isfile(edit_graph_filename):
        # This is slow to compute on complex sample sets
        sys.stderr.write(f"WARNING: Skipping {edit_graph_filename} as already exists\n")
    else:
        # The XGMML output has minimal dependencies compared to PDF output
        return_code = edit_graph(
            graph_output=edit_graph_filename,
            graph_format="xgmml",
            db_url=db,
            inputs=fasta_files + classified_files,
            method=args.method,
            min_abundance=args.abundance,
            # total_min_abundance=args.total,
            always_show_db=args.showdb,
            # max_edit_dist=args.editdist,
            ignore_prefixes=tuple(args.ignore_prefixes),
            debug=args.verbose,
        )
        if return_code == 2:
            # Special value indicating graph skipped as too large
            pass
        elif return_code:
            sys.stderr.write("ERROR: Pipeline aborted during edit-graph\n")
            sys.exit(return_code)
        sys.stderr.write(f"Wrote {edit_graph_filename}\n")

    sys.stderr.write("All done!\n")


def ena_submit(args=None):
    """Subcommand to run multiple-output-folder summary at sample level."""
    from .ena_submit import main

    if args.metadata:
        check_input_file(args.metadata)
        if not args.metacols:
            sys.exit("ERROR: Must also supply -c / --metacols argument.")
    if "\t" in args.library:
        sys.exit("ERROR: Can't use tabs in --library argument")
    if "\t" in args.instrument:
        sys.exit("ERROR: Can't use tabs in --instrument argument")
    if "\t" in args.design:
        sys.exit("ERROR: Can't use tabs in --design argument")
    if "\t" in args.protocol:
        sys.exit("ERROR: Can't use tabs in --protocol argument")
    return main(
        fastq=args.input,
        output=args.output,
        metadata_file=args.metadata,
        metadata_cols=args.metacols,
        metadata_fieldnames=args.metafields,
        metadata_index=args.metaindex,
        ignore_prefixes=args.ignore_prefixes,
        library_name=args.library,
        instrument_model=args.instrument,
        design_description=args.design,
        library_construction_protocol=args.protocol,
        insert_size=args.insert,
        tmp_dir=args.temp,
        debug=args.verbose,
    )


# Common arguments
# ================

# "-d", "--database",
ARG_DB_INPUT = dict(  # noqa: C408
    type=str, default="-", help="Marker database to use, default '-' for bundled DB."
)

# "-i", "--input",
ARG_INPUT_FASTA = dict(  # noqa: C408
    type=str,
    required=True,
    nargs="+",
    help="One or more prepared FASTA filenames or folder names "
    "(containing files named *.fasta).",
)

# "--ignore-prefixes",
ARG_IGNORE_PREFIXES = dict(  # noqa: C408
    type=str,
    nargs="+",
    metavar="STEM",
    default=IGNORE_PREFIXES,
    help="One or more filename prefixes to ignore, default %s"
    % ", ".join(repr(_) for _ in IGNORE_PREFIXES),
)

# "-m", "--method",
ARG_METHOD_OUTPUT = dict(  # noqa: C408
    type=str,
    default=DEFAULT_METHOD,
    choices=sorted(method_classifier),
    help=f"Classify method to run, default is '{DEFAULT_METHOD}'.",
)

# "-m", "--method",
ARG_METHOD_INPUT = dict(  # noqa: C408
    type=str,
    default=DEFAULT_METHOD,
    choices=sorted(method_classifier),
    help=f"Classify method (to infer filenames), default '{DEFAULT_METHOD}'.",
)

# "--merged-cache",
ARG_MERGED_CACHE = dict(  # noqa: C408
    type=str,
    required=False,
    metavar="DIRNAME",
    help="Advanced option. Cache directory for temporary sample FASTA files after "
    "trimming and merging, but before primer trimming and abundance threshold. "
    "Intended for multi-primer analyses from a pooled sequencing run.",
)

# "-t", "--temp",
ARG_TEMPDIR = dict(  # noqa: C408
    type=str,
    required=False,
    metavar="DIRNAME",
    help="Debug option. Specify an (ideally empty) directory to "
    "use for temporary files, which will not be deleted.",
)

# "-v", "--verbose",
ARG_VERBOSE = dict(action="store_true", help="Verbose logging.")  # noqa: C408

# "--cpu",
ARG_CPU = dict(  # noqa: C408
    type=int, default=0, help="Number of parallel threads to use in called tools."
)

# Common import arguments
# =======================

# "-s", "--sep",
ARG_FASTA_SEP = dict(  # noqa: C408
    type=str,
    default=";",
    metavar="CHAR",
    help="FASTA description entry separator, default semi-colon.",
)

# "-d", "--database",
ARG_DB_WRITE = dict(  # noqa: C408
    type=str, required=True, help="Which database to write to (or create)."
)

# "-n", "--name",
ARG_NAME = dict(  # noqa: C408
    type=str, default="", help="Name to record for this data source (string)."
)

# "-x", "--lax",
ARG_LAX = dict(  # noqa: C408
    default=False,
    action="store_true",
    help="Accept species names without pre-loaded taxonomy.",
)

# "-g", "--genus"
ARG_GENUS_ONLY = dict(  # noqa: C408
    default=False,
    action="store_true",
    help="Record at genus level only (and only validate at genus level, "
    "unless using -x / --lax in which case anything is accepted as a genus).",
)

# Prepare reads arguments
# =======================

# "--minlen",
ARG_MIN_LENGTH = dict(  # noqa: C408
    type=int,
    default=100,
    metavar="LENGTH",
    help="Minimum length sequence to accept (default 100).",
)

# "--maxlen",
ARG_MAX_LENGTH = dict(  # noqa: C408
    type=int,
    default=1000,
    metavar="LENGTH",
    help="Maximum length sequence to accept (default 1000).",
)

# "--flip",
ARG_FLIP = dict(  # noqa: C408
    default=False,
    action="store_true",
    help="Also check reverse complement strand for primers.",
)

# "-l", "--left",
ARG_PRIMER_LEFT = dict(  # noqa: C408
    type=str,
    default="GAAGGTGAAGTCGTAACAAGG",
    metavar="PRIMER",
    help="Left primer sequence, find and remove from start of "
    "marker sequence. Can use IUPAC ambiguity codes. "
    "Default 21bp ITS6 'GAAGGTGAAGTCGTAACAAGG' from Cooke "
    "et al. 2000 https://doi.org/10.1006/fgbi.2000.1202",
)
ARG_PRIMER_LEFT_BLANK = dict(  # noqa: C408
    type=str, default="", metavar="PRIMER", help="Left primer sequence, default blank."
)

# "-r", "--right",
ARG_PRIMER_RIGHT = dict(  # noqa: C408
    type=str,
    default="GCARRGACTTTCGTCCCYRC",
    metavar="PRIMER",
    help="Right primer sequence, find and remove reverse "
    "complement from end of marker sequence. Can use "
    "IUPAC ambiguity codes. Default 20bp 5.8S-1R primer "
    "'GCARRGACTTTCGTCCCYRC' from Scibetta et al. 2012 "
    "https://doi.org/10.1016/j.mimet.2011.12.012 - meaning "
    "look for 'GYRGGGACGAAAGTCYYTGC' after marker.",
)
ARG_PRIMER_RIGHT_BLANK = dict(  # noqa: C408
    type=str, default="", metavar="PRIMER", help="Right primer sequence, default blank."
)

# "--hmm",
ARG_HMM = dict(  # noqa: C408
    type=str,
    default="-",
    metavar="PATH",
    help="Location of negative control HMMER3 Hidden Markov Model "
    "file, filename stem without the '.h3i', '.h3f', etc extension. "
    "Use '' for none, or '-' for default four synthetic controls (default).",
)

# Common pipeline arguments
# =========================

# "-i", "--input",
ARG_INPUT_FASTQ = dict(  # noqa: C408
    type=str,
    required=True,
    nargs="+",
    help="One or more paired FASTQ filenames or folder names "
    "(containing files named *.fastq or *.fastq.gz).",
)

# "-n", "--negctrls",
ARG_CONTROLS = dict(  # noqa: C408
    type=str,
    nargs="+",
    # Does accept folder names, but kind of pointless...
    help="One or more negative control paired FASTQ filenames "
    "(may also appear in the inputs). High marker levels "
    "will increase the minimum abundance threshold of other "
    "FASTQ files in the same folder. Can use '-' for none.",
)

# "-a", "--abundance",
ARG_FASTQ_MIN_ABUNDANCE = dict(  # noqa: C408
    type=int,
    default=str(DEFAULT_MIN_ABUNDANCE),
    help="Mininum abundance applied to unique marker sequences in each sample"
    f" (i.e. each FASTQ pair), default {DEFAULT_MIN_ABUNDANCE}."
    " May be increased based on negative controls.",
)

# Common metadata arguments
# ========================

# "-t", "--metadata",
ARG_METADATA = dict(  # noqa: C408
    type=str,
    default="",
    metavar="METADATAFILE",
    help="Optional tab separated table containing metadata indexed by "
    "sample name. Must also specify the columns with -c / --metacols. ",
)

# "-c", "--metacols",
ARG_METACOLS = dict(  # noqa: C408
    type=str,
    default="",
    metavar="COLUMNS",
    help="Comma separated list (e.g, '1,3,5') of columns from the metadata "
    "table specified with the -m / --metadata argument to be included in the "
    "report header, and use to sort the samples. Use in conjunction with "
    "-m / --metadata argument.",
)

# "-x", "--metaindex",
ARG_METAINDEX = dict(  # noqa: C408
    type=int,
    default="1",
    metavar="COL",
    help="If using metadata, which column contains the sequenced sample "
    "names. Default 1. Field can contain multiple semi-colon separated names "
    "catering to the fact that a field sample could be sequenced multiple "
    "times with technical replicates.",
)

# "-g", "--metagroups",
ARG_METAGROUPS = dict(  # noqa: C408
    type=int,
    default="0",
    metavar="COL",
    help="If using metadata, which column values should be used for applying "
    "background color bands.  All samples with the same metadata value must "
    "be grouped together after sorting, as the colors are reused. "
    "Zero (default) is interpretted as the first column requested as metadata "
    "output (the primary sorting key, ensuring all members of the same group "
    "will be together).",
)

# "-f", "--metafields",
ARG_METAFIELDS = dict(  # noqa: C408
    type=int,
    default="1",
    metavar="ROW",
    help="If using metadata, which row should be used as the field names? "
    "Default 1, use 0 for no labels. "
    "Use in conjunction with -m / --metadata argument.",
)

# "-q", "--requiremeta",
ARG_REQUIREMETA = dict(  # noqa: C408
    action="store_true",
    help="Ignore any input files without metadata for report.",
)

# "-u", "--unsequenced",
ARG_UNSEQUENCED = dict(  # noqa: C408
    action="store_true",
    help="Show unsequenced metadata entries in sample report.",
)
# Command line definition
# =======================


def main(args=None):
    """Execute the command line script thapbi_pict.

    Note this needs to cope with no command line information via the Python
    function arguments because that is how it is invoked by ``thapbi_pict``
    via setuptools.

    To facilitate automated testing etc, it can also be called with an
    argument list.
    """
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(
        prog="thapbi_pict",
        description=(
            f"THAPBI Phytophthora ITS1 Classifier Tool (PICT), v{__version__}."
        ),
        epilog="e.g. run 'thapbi_pict pipeline -h' for the pipeline subcommand help. "
        "Please see https://thapbi-pict.readthedocs.io/ for documentation.",
    )
    parser.add_argument(
        "-v", "--version", action="version", version=f"THAPBI PICT v{__version__}"
    )
    subparsers = parser.add_subparsers(
        title="subcommands", help="Each subcommand has its own additional help."
    )

    # pipeline (listing first as likely to be the most used subcommand)
    subcommand_parser = subparsers.add_parser(
        "pipeline",
        description="Run default classification pipeline on paired FASTQ files.",
        epilog="This is equivalent to running the individual stages (prepare-reads, "
        "classify, summary, edit-graph) with their defaults, with only a minority of "
        "settings available here.",
    )
    subcommand_parser.add_argument("-i", "--input", **ARG_INPUT_FASTQ)
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument("-n", "--negctrls", **ARG_CONTROLS)
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        metavar="DIRNAME",
        help="Output directory. Required.",
    )
    subcommand_parser.add_argument(
        "-r",
        "--report",
        type=str,
        metavar="STEM",
        help="Stem for generating report filenames.",
    )
    subcommand_parser.add_argument(
        "-s",
        "--sampleout",
        type=str,
        default="",
        metavar="DIRNAME",
        help="Output directory for intermediate files for each sample "
        "(FASTQ pair). Defaults to -o / --output.",
    )
    subcommand_parser.add_argument("-a", "--abundance", **ARG_FASTQ_MIN_ABUNDANCE)
    subcommand_parser.add_argument("-d", "--database", **ARG_DB_INPUT)
    subcommand_parser.add_argument("--hmm", **ARG_HMM)
    subcommand_parser.add_argument("--minlen", **ARG_MIN_LENGTH)
    subcommand_parser.add_argument("--maxlen", **ARG_MAX_LENGTH)
    # Not using -l and -r for primers as used -r for report:
    subcommand_parser.add_argument("--left", **ARG_PRIMER_LEFT)
    subcommand_parser.add_argument("--right", **ARG_PRIMER_RIGHT)
    subcommand_parser.add_argument("--flip", **ARG_FLIP)
    subcommand_parser.add_argument("-m", "--method", **ARG_METHOD_OUTPUT)
    subcommand_parser.add_argument("-t", "--metadata", **ARG_METADATA)
    subcommand_parser.add_argument("-c", "--metacols", **ARG_METACOLS)
    subcommand_parser.add_argument("-x", "--metaindex", **ARG_METAINDEX)
    subcommand_parser.add_argument("-g", "--metagroups", **ARG_METAGROUPS)
    subcommand_parser.add_argument("-f", "--metafields", **ARG_METAFIELDS)
    # Can't use -s for --showdb as already used for sample intermediates
    subcommand_parser.add_argument(
        "--showdb",
        action="store_true",
        help="Show DB entries in edit-graph, regardless of their abundance "
        "in samples. Very slow with large database.",
    )
    subcommand_parser.add_argument("--merged-cache", **ARG_MERGED_CACHE)
    subcommand_parser.add_argument("-q", "--requiremeta", **ARG_REQUIREMETA)
    subcommand_parser.add_argument("-u", "--unsequenced", **ARG_UNSEQUENCED)
    # Can't use -t for --temp as already using for --metadata:
    subcommand_parser.add_argument("--temp", **ARG_TEMPDIR)
    subcommand_parser.add_argument("--cpu", **ARG_CPU)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=pipeline)
    del subcommand_parser

    # load-tax
    subcommand_parser = subparsers.add_parser(
        "load-tax", description="Load an NCBI taxonomy dump into an ITS1 database."
    )
    subcommand_parser.add_argument(
        "-t",
        "--tax",
        type=str,
        required=True,
        metavar="DIRNAME",
        help="Folder containing NCBI taxonomy files 'names.dmp' etc.",
    )
    subcommand_parser.add_argument("-d", "--database", **ARG_DB_WRITE)
    subcommand_parser.add_argument(
        "-a",
        "--ancestors",
        type=str,
        default="4762",
        help="Comma separated list of NCBI taxids at genus level or higher. "
        "Default is 4762 for Oomycetes, use 4776 for Peronosporales, or "
        "4783 for Phytophthora only.",
    ),
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=load_tax)
    del subcommand_parser  # To prevent acidentally adding more

    # ncbi-import
    subcommand_parser = subparsers.add_parser(
        "ncbi-import",
        description="Load an NCBI format ITS1 FASTA file into a database. "
        "By default verifies species names against a pre-loaded taxonomy, "
        "non-matching entries are rejected. In lax mode uses heuristics to "
        "split the species and any following free text - see also the "
        "curated-import command which avoids this ambiguity.",
    )
    subcommand_parser.add_argument(
        "-i", "--input", type=str, required=True, help="One ITS1 fasta filename."
    )
    subcommand_parser.add_argument("-d", "--database", **ARG_DB_WRITE)
    subcommand_parser.add_argument("--minlen", **ARG_MIN_LENGTH)
    subcommand_parser.add_argument("--maxlen", **ARG_MAX_LENGTH)
    subcommand_parser.add_argument("-n", "--name", **ARG_NAME)
    subcommand_parser.add_argument("-x", "--lax", **ARG_LAX)
    subcommand_parser.add_argument("-g", "--genus", **ARG_GENUS_ONLY)
    subcommand_parser.add_argument("-l", "--left", **ARG_PRIMER_LEFT)
    subcommand_parser.add_argument("-r", "--right", **ARG_PRIMER_RIGHT)
    subcommand_parser.add_argument(
        "-s",
        "--sep",
        type=str,
        default=CTRL_A,
        metavar="CHAR",
        help="FASTA description entry separator, default Ctrl+A.",
    )
    subcommand_parser.add_argument("-t", "--temp", **ARG_TEMPDIR)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=ncbi_import)
    del subcommand_parser  # To prevent acidentally adding more

    # curated-seq
    subcommand_parser = subparsers.add_parser(
        "curated-seq",
        description="Prepare FASTA files with curated species name (for "
        "later import to a database with 'thapbi_pict curated-import) "
        "from prepared single isolate controls (FASTA and TSV files).",
    )
    subcommand_parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        nargs="+",
        help="One or more ITS1 FASTA and classifier filenames or folders "
        "(names containing files named *.fasta and *.method.tsv, where "
        "method is set via the -m / --method argument).",
    )
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="DIRNAME",
        help="Directory for FASTA files with species names in records. "
        "Default '-' for stdout.",
    )
    subcommand_parser.add_argument(
        "-m",
        "--method",
        type=str,
        default="known",
        help="Method name used to determines TSV filenames from which the "
        "species classification will be read. Default is 'known' matching "
        "the convention used in the classifier assessment command for "
        "known trusted species assignments (i.e. positive controls).",
    )
    subcommand_parser.add_argument(
        "-a",
        "--abundance",
        type=int,
        default=str(DEFAULT_MIN_ABUNDANCE * 10),
        help=(
            "Mininum abundance to require before importing a sequence, "
            "over-and-above whatever was used to prepare the FASTA file. "
            f"Default here is {DEFAULT_MIN_ABUNDANCE * 10}, ten times the default "
            f" of {DEFAULT_MIN_ABUNDANCE} used for the classification pipeline - "
            "be cautious what goes in your ITS1 database)."
        ),
    )
    subcommand_parser.add_argument("-t", "--temp", **ARG_TEMPDIR)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=curated_seq)
    del subcommand_parser  # To prevent acidentally adding more

    # curated-import
    subcommand_parser = subparsers.add_parser(
        "curated-import",
        description="Load curated marker FASTA file(s) into a database. "
        "Treats the FASTA title lines as identitier, space, species. "
        "Also supports multiple entries in a single record using a FASTA "
        "title line of identifier1, space, species1, separator, identifier2, "
        "space, species2, etc. "
        "By default verifies species names against a pre-loaded taxonomy, "
        "non-matching entries are rejected.",
    )
    subcommand_parser.add_argument("-i", "--input", **ARG_INPUT_FASTA)
    subcommand_parser.add_argument("-d", "--database", **ARG_DB_WRITE)
    subcommand_parser.add_argument("--minlen", **ARG_MIN_LENGTH)
    subcommand_parser.add_argument("--maxlen", **ARG_MAX_LENGTH)
    subcommand_parser.add_argument("-n", "--name", **ARG_NAME)
    subcommand_parser.add_argument("-x", "--lax", **ARG_LAX)
    subcommand_parser.add_argument("-g", "--genus", **ARG_GENUS_ONLY)
    subcommand_parser.add_argument("-l", "--left", **ARG_PRIMER_LEFT_BLANK)
    subcommand_parser.add_argument("-r", "--right", **ARG_PRIMER_RIGHT_BLANK)
    subcommand_parser.add_argument("-s", "--sep", **ARG_FASTA_SEP)
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=curated_import)
    del subcommand_parser  # To prevent acidentally adding more

    # dump
    subcommand_parser = subparsers.add_parser(
        "dump",
        description="Export an ITS1 database to a text file.",
        epilog="e.g. 'thapbi_pict dump -d ... -g Peronospora -o Peronospora.txt'",
    )
    subcommand_parser.add_argument("-d", "--database", **ARG_DB_INPUT)
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        help="File to write to (default '-' meaning stdout)",
    )
    subcommand_parser.add_argument(
        "-m",
        "--minimal",
        action="store_true",
        help="Minimal output, one record per trimmed sequence with MD5 and species "
        "list - rather than all entries with original ID and full sequence.",
    )
    subcommand_parser.add_argument(
        "-f",
        "--format",
        type=str,
        default="txt",
        choices=["txt", "fasta"],
        help="Format to write out (default 'txt' for debugging).",
    )
    subcommand_parser.add_argument(
        "-g",
        "--genus",
        type=str,
        default="",
        help="Which genus (or genera) to export (comma separated list). "
        "Default is not to filter by genus.",
    )
    subcommand_parser.add_argument(
        "-s",
        "--species",
        type=str,
        default="",
        help="Which species to export (comma separated list, ignores any "
        "spaces after each comma). Requires a single genus argument be given. "
        "Default is not to filter by species.",
    )
    subcommand_parser.add_argument("--sep", **ARG_FASTA_SEP)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=dump)
    del subcommand_parser  # To prevent acidentally adding more

    # conflicts
    subcommand_parser = subparsers.add_parser(
        "conflicts",
        description="Count genus or species conflicts in a marker database.",
        epilog="Number of genus level conflicts is used as the return code. "
        "e.g. 'thapbi_pict conflicts -d ... -o conflicts.txt ; echo $?'",
    )
    subcommand_parser.add_argument("-d", "--database", **ARG_DB_INPUT)
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        help="File to write to (default '-' meaning stdout)",
    )
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=conflicts)
    del subcommand_parser

    # prepare reads
    subcommand_parser = subparsers.add_parser(
        "prepare-reads",
        description="Trim and merge paired FASTQ files of marker amplicons.",
        epilog="Each pair of input files should follow the naming style"
        "XXX_1.fastq[.gz] and XXX_2.fastq[.gz], or "
        "XXX_R1.fastq[.gz] and XXX_R2.fastq[.gz], or "
        "XXX_R1_001.fastq[.gz] and XXX_R2_001.fastq[.gz], and will "
        "result an output file XXX.fasta in the specified ouput"
        "directory (which defaults to the FASTQ directory).\n\n"
        "The output FASTA files are non-redundant, records named by "
        "checksum and their abundance, and sorted by decreasing "
        "abundance then alphabetically by sequence.",
    )
    subcommand_parser.add_argument("-i", "--input", **ARG_INPUT_FASTQ)
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument("-n", "--negctrls", **ARG_CONTROLS)
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="DIRNAME",
        help="Directory to write output FASTA files to, "
        "default is next to each input file.",
    )
    subcommand_parser.add_argument("-a", "--abundance", **ARG_FASTQ_MIN_ABUNDANCE)
    subcommand_parser.add_argument("--hmm", **ARG_HMM)
    subcommand_parser.add_argument("--minlen", **ARG_MIN_LENGTH)
    subcommand_parser.add_argument("--maxlen", **ARG_MAX_LENGTH)
    subcommand_parser.add_argument(
        "-p",
        "--primers",
        type=str,
        default="",
        metavar="DIRNAME",
        help="Where to write optional failed primer FASTA files.",
    )
    subcommand_parser.add_argument("-l", "--left", **ARG_PRIMER_LEFT)
    subcommand_parser.add_argument("-r", "--right", **ARG_PRIMER_RIGHT)
    subcommand_parser.add_argument("--flip", **ARG_FLIP)
    subcommand_parser.add_argument("--merged-cache", **ARG_MERGED_CACHE)
    subcommand_parser.add_argument("-t", "--temp", **ARG_TEMPDIR)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.add_argument("--cpu", **ARG_CPU)
    subcommand_parser.set_defaults(func=prepare_reads)
    del subcommand_parser  # To prevent acidentally adding more

    subcommand_parser = subparsers.add_parser(
        "fasta-nr",
        description="Prepare non-redundant FASTA file using MD5 naming.",
        epilog="Each unique sequence will be named <MD5>_<count> using "
        "the MD5 checksum of the upper case sequence and its total "
        "abundance. Input files may use <prefix>_<count> naming, this "
        "count will be used for the associated sequence.",
    )
    subcommand_parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=False,
        nargs="+",
        metavar="INPUT",
        help="One or more FASTA files.",
    )
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument(
        "-r",
        "--revcomp",
        type=str,
        required=False,
        metavar="INPUT",
        help="One or more FASTA files as input "
        "to be reverse-complemented on loading.",
    )
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="PATH",
        help="Single output filename, '-' for stdout (default). "
        "Can be a directory if a single -i/-r input file given.",
    )
    subcommand_parser.add_argument(
        "-a",
        "--abundance",
        type=int,
        default=0,
        help="Mininum abundance to require before outputing a sequence. "
        "Default no minimum.",
    )
    subcommand_parser.add_argument("--minlen", **ARG_MIN_LENGTH)
    subcommand_parser.add_argument("--maxlen", **ARG_MAX_LENGTH)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=fasta_nr)

    # classify
    subcommand_parser = subparsers.add_parser(
        "classify",
        description="Classify FASTA file of marker sequences by species.",
        epilog="Each input file XXX.fasta will result in an output file "
        "named XXX.method.tsv in the specified output directory (default "
        "input dir).",
    )
    subcommand_parser.add_argument("-i", "--input", **ARG_INPUT_FASTA)
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument(
        "-a",
        "--abundance",
        type=int,
        default=0,
        help="Mininum abundance applied to unique marker sequences in each "
        "FASTA sample file, default 0 (classify all).",
    )
    subcommand_parser.add_argument("-d", "--database", **ARG_DB_INPUT)
    subcommand_parser.add_argument("-m", "--method", **ARG_METHOD_OUTPUT)
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="",
        metavar="DIRNAME",
        help="Directory for output reports. Default '' means next to "
        "each input file. Use '-' for stdout.",
    )
    subcommand_parser.add_argument("-t", "--temp", **ARG_TEMPDIR)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.add_argument("--cpu", **ARG_CPU)
    subcommand_parser.set_defaults(func=classify)
    del subcommand_parser  # To prevent acidentally adding more

    # assess-classification
    subcommand_parser = subparsers.add_parser(
        "assess",
        description="Assess accuracy of marker sequence classification.",
        epilog="Takes as input predictions named XXX.method.tsv "
        "and matching expected classifications in XXX.known.tsv "
        "(which can be in different directories) to produce a "
        "multi-species confusion matrix (output on request) and "
        "classifier performance metrics (to stdout by default). "
        "You can deliberately compare two prediction methods to "
        "each other using this, but a known set of positive controls "
        "is the expected benchmark. The assessment is at sample level "
        "(taking the union of species predicted by all sequences from "
        "each sample).",
    )
    subcommand_parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        nargs="+",
        help="One or more prediction file or folder names. Expects to "
        "find matching files *.method.tsv to be assessed against "
        "*.known.tsv, where these filenames can be set via "
        "-m / --method and -k / --known arguments. ",
    )
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument(
        "-k",
        "--known",
        type=str,
        default="known",
        help="Replaces the string used in filenames for the truth "
        "against which the method in -m / --method is assessed. "
        "This could be any defined method, default is 'known'.",
    )
    arg = subcommand_parser.add_argument("-d", "--database", **ARG_DB_INPUT)
    arg.help += " Used for species list to infer true negatives."
    del arg
    subcommand_parser.add_argument("-m", "--method", **ARG_METHOD_INPUT)
    subcommand_parser.add_argument(
        "-a",
        "--abundance",
        type=int,
        default="1",
        help="Mininum abundance to require before considering a classification. "
        "Default is one meaning look at everything, but rather than re-running "
        "the classifier with a stricter minimum abundance you can apply it here.",
    )
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="FILENAME",
        help="File to write species level classification assessment table to. "
        "Default is '-' meaning to stdout.",
    )
    subcommand_parser.add_argument(
        "-t",
        "--table",
        type=str,
        metavar="FILENAME",
        help="File to write expected-to-predicted mapping tally table to. "
        "Can use '-' meaning to stdout. Default is not to write this file.",
    )
    subcommand_parser.add_argument(
        "-c",
        "--confusion",
        type=str,
        metavar="FILENAME",
        help="File to write species level confusion matrix to. "
        "Can use '-' meaning to stdout. Default is not to write this file.",
    )
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=assess_classification)
    del subcommand_parser  # To prevent acidentally adding more

    # summary
    subcommand_parser = subparsers.add_parser(
        "summary",
        description="Sample and sequence summary reports on classifier output.",
        epilog="Assumes you've run prepare-reads and classify, and have "
        "folders with XXX.fasta and XXX.method.tsv files from your plate(s). "
        "The output is two sets of tables. The read tables have one row per "
        "unique sequence (as trimmed by the prepare-reads step, can be "
        "thousands of rows) and one column per sample (typically 96 samples "
        "per plate). The sample tables have one row per sample, and one "
        "column per genus and species. Plus a plain text sample report.",
    )
    subcommand_parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        nargs="+",
        help="One or more prepared read files (*.fasta), prediction "
        "files (*.method.tsv) or folder names. If passing folder names, "
        "it expects to find paired files using these extensions. "
        "The classifier method extension can be set via -m / --method.",
    )
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        metavar="DIRNAME",
        help="Output directory. Required.",
    )
    subcommand_parser.add_argument(
        "-r",
        "--report",
        type=str,
        metavar="STEM",
        help="Stem for generating report filenames.",
    )
    subcommand_parser.add_argument(
        "-m",
        "--method",
        type=str,
        default=DEFAULT_METHOD,
        help="Classifier method(s) to report, comma separated list "
        f"(used to infer filenames), default is '{DEFAULT_METHOD}' (only).",
    )
    subcommand_parser.add_argument(
        "-a",
        "--abundance",
        type=int,
        default=str(DEFAULT_MIN_ABUNDANCE),
        help="Mininum sample level abundance to require for the report. "
        f"Default {DEFAULT_MIN_ABUNDANCE} reflects default in prepare-reads. "
        "Rather than re-running the prepare or classifier steps with a stricter "
        "minimum abundance you can apply it here. Use zero or one to look at "
        "everything (but beware that negative control samples will include low "
        "abundance entries).",
    )
    subcommand_parser.add_argument("-t", "--metadata", **ARG_METADATA)
    subcommand_parser.add_argument("-c", "--metacols", **ARG_METACOLS)
    subcommand_parser.add_argument("-x", "--metaindex", **ARG_METAINDEX)
    subcommand_parser.add_argument("-g", "--metagroups", **ARG_METAGROUPS)
    subcommand_parser.add_argument("-f", "--metafields", **ARG_METAFIELDS)
    subcommand_parser.add_argument("-q", "--requiremeta", **ARG_REQUIREMETA)
    subcommand_parser.add_argument("-u", "--unsequenced", **ARG_UNSEQUENCED)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=summary)
    del subcommand_parser  # To prevent acidentally adding more

    # edit-graph
    subcommand_parser = subparsers.add_parser(
        "edit-graph",
        description="Draw network graph of marker sequences using edit distance.",
        epilog="Takes an ITS1 database and/or prepared FASTA files as input. "
        "The output is a network graph (in a choice of format) with unique "
        "sequences as nodes (in the PDF labelled by the database taxonomy, "
        "colored by genus, size set by total abundance in the FASTA files), "
        "and short edit distances as edges between nodes. For Cytoscape "
        "we recommend generating XGMML output here, the start Cytoscape, "
        "menu 'File', 'Import', 'Import from file', and then run a layout. "
        "Both 'Perfuse Force Directed' and 'Edge-weighted Spring Embedded' "
        "work well.",
    )
    arg = subcommand_parser.add_argument("-d", "--database", **ARG_DB_INPUT)
    arg.help += " Used for labels and colors. Use '' to mean no DB."
    del arg
    # Currently ARG_INPUT_FASTA uses required=True, but we need to change that:
    arg = subcommand_parser.add_argument("-i", "--input", **ARG_INPUT_FASTA)
    arg.required = False
    del arg
    subcommand_parser.add_argument(
        "-m",
        "--method",
        type=str,
        default=DEFAULT_METHOD,
        choices=sorted(method_classifier) + ["-"],
        help="Optional classifier method to annotate sequences sequences with. "
        f"Default is {DEFAULT_METHOD}, use '-' for none.",
    )
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument(
        "-a",
        "--abundance",
        type=int,
        default=str(DEFAULT_MIN_ABUNDANCE),
        help="Mininum sample level abundance for FASTA sequences. "
        f"Default {DEFAULT_MIN_ABUNDANCE} reflects default in prepare-reads.",
    )
    subcommand_parser.add_argument(
        "-t",
        "--total",
        type=int,
        default="0",
        help="Mininum total abundance for FASTA sequences. "
        "Applied after per-sample level minimum (-a / --abundance). "
        "Offered as a way to simplify the final graph.",
    )
    subcommand_parser.add_argument(
        "-s",
        "--showdb",
        action="store_true",
        help="Show DB entries, regardless of their abundance in the "
        "FASTA inputs. Required if only drawing DB entries.",
    )
    subcommand_parser.add_argument(
        "-e",
        "--editdist",
        type=int,
        default="3",
        choices=[1, 2, 3],
        help="Maximum edit distance to draw. Default and maximum 3.",
    )
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="FILENAME",
        help="Write graph here. Default is '-' meaning stdout.",
    )
    subcommand_parser.add_argument(
        "-f",
        "--format",
        type=str,
        default="xgmml",
        choices=["graphml", "gexf", "gml", "xgmml", "pdf"],
        help="Format to write out (default 'xgmml' for Cytoscape).",
    )
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=edit_graph)
    del subcommand_parser  # To prevent accidentally adding more

    # submit samples to ENA
    subcommand_parser = subparsers.add_parser(
        "ena-submit",
        description="Write paired FASTQ table for upload to ENA (and thus SRA).",
        epilog="Facilitate interactive upload of paired FASTQ files to the ENA/SRA. "
        "First upload the raw FASTQ files, and interactively define your samples. "
        "Then use the output of this command to match the FASTQ files to the samples. "
        "If your ENA sample names match your FASTQ prefixes, no metadata file is "
        "needed here. However, you can provide a metadata table to map the FASTQ "
        "stem to your sample aliases or ENA sample accessions.",
    )
    subcommand_parser.add_argument("-i", "--input", **ARG_INPUT_FASTQ)
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        help="File to write to (default '-' meaning stdout)",
    )
    subcommand_parser.add_argument("-t", "--metadata", **ARG_METADATA)
    subcommand_parser.add_argument("-c", "--metacols", **ARG_METACOLS)
    subcommand_parser.add_argument("-x", "--metaindex", **ARG_METAINDEX)
    subcommand_parser.add_argument("-f", "--metafields", **ARG_METAFIELDS)
    subcommand_parser.add_argument(
        "--library",
        type=str,
        default="-",
        help="Value for library_name field, default is '-' meaning take the "
        "FASTQ file's parent folder name, which is intended to support "
        "configurations like using multiple 96-well plates.",
    )
    subcommand_parser.add_argument(
        "--instrument",
        type=str,
        default="Illumina MiSeq",
        help="Value for instrument_model field, default 'Illumina MiSeq'.",
    )
    subcommand_parser.add_argument(
        "--protocol",
        type=str,
        default="",
        help="Value for optional library_construction_protocol field, default blank.",
    )
    subcommand_parser.add_argument(
        "--design",
        type=str,
        default="",
        help="Value for optional design_description field, default blank.",
    )
    subcommand_parser.add_argument(
        "--insert",
        type=int,
        default=250,
        help="Value for nominal/average insert_size, default 250.",
    )
    # Can't use -t for --temp as already using for --metadata:
    subcommand_parser.add_argument("--temp", **ARG_TEMPDIR)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=ena_submit)
    del subcommand_parser

    # What have we been asked to do?
    options = parser.parse_args(args)
    if hasattr(options, "func"):
        # Invoke the subcommand
        sys.exit(options.func(options))
    else:
        # Called without a subcommand
        parser.print_help()


if __name__ == "__main__":
    main()
