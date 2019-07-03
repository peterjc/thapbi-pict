# Copyright 2018-2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

"""This defines the thapbi_pict command line tool.

This works via ``setup.py`` where under ``entry_points`` we define a
``console_scripts`` entry for ``thapbi_pict`` (executable name) pointing to
the ``main()`` function define in this Python file.
"""

import argparse
import os
import sys

from . import __version__
from .classify import method_classify_file as method_classifier


# Common command line defaults
# ============================


DEFAULT_METHOD = "onebp"
DEFAULT_MIN_ABUNDANCE = 100


# Argument validation functions
# =============================


def check_input_file(filename):
    """Command line validation of an input filename."""
    if not os.path.isfile(filename):
        sys.exit("ERROR: Could not find input file: %s" % filename)


def check_output_directory(out_dir):
    """Command line validation of output directory value."""
    if out_dir == "-" or os.path.isdir(out_dir):
        return True
    elif os.path.isfile(out_dir):
        sys.exit("ERROR: Output directory name is a file: %s\n" % out_dir)
    else:
        sys.exit("ERROR: Output directory does not exist")


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
            sys.exit(
                "ERROR: Using hyphen as a database default is not supported here.\n"
            )
    if not text:
        sys.exit("ERROR: The database argument is required.\n")
    prefix = "sqlite:///"
    if text.startswith(prefix):
        db = text[len(prefix) :]
        assert text == prefix + db
    else:
        db = text
    if exist and db != ":memory:" and not os.path.isfile(db):
        sys.exit("ERROR: The database %s was not found.\n" % db)
    return prefix + db


# Subcommand dispatch
# ===================


def load_tax(args=None):
    """Subcommand to load an NCBI taxonomy dump into an ITS1 database."""
    from .taxdump import main

    return main(
        tax=args.tax,
        db_url=expand_database_argument(args.database),
        ancestors=args.ancestors,
        debug=args.verbose,
    )


def ncbi_import(args=None):
    """Subcommand to import an NCBI ITS1 FASTA file into a database."""
    from .ncbi import main

    return main(
        fasta_file=args.input,
        db_url=expand_database_argument(args.database),
        name=args.name,
        validate_species=not args.lax,
        genus_only=args.genus,
        debug=args.verbose,
    )


def seq_import(args=None):
    """Subcommand to import classified ITS1 sequences into a database."""
    from .seq_import import main

    return main(
        inputs=args.input,
        method=args.method,
        db_url=expand_database_argument(args.database),
        min_abundance=args.abundance,
        name=args.name,
        validate_species=not args.lax,
        genus_only=args.genus,
        debug=args.verbose,
    )


def legacy_import(args=None):
    """Subcommand to import a legacy ITS1 FASTA file into a database."""
    from .legacy import main

    return main(
        fasta_file=args.input,
        db_url=expand_database_argument(args.database),
        name=args.name,
        validate_species=not args.lax,
        genus_only=args.genus,
        debug=args.verbose,
    )


def dump(args=None):
    """Subcommand to dump a database to a text file."""
    from .dump import main

    return main(
        db_url=expand_database_argument(args.database, exist=True, hyphen_default=True),
        output_filename=args.output,
        output_format=args.format,
        clade=args.clade,
        genus=args.genus,
        species=args.species,
        debug=args.verbose,
    )


def prepare_reads(args=None):
    """Subcommand to prepare FASTA paired reads."""
    from .prepare import main

    check_output_directory(args.output)
    if args.temp:
        check_output_directory(args.temp)
    return_code = main(
        fastq=args.input,
        negative_controls=args.negctrls,
        out_dir=args.output,
        primer_dir=args.primers,
        left_primer=args.left,
        right_primer=args.right,
        min_abundance=args.abundance,
        tmp_dir=args.temp,
        debug=args.verbose,
        cpu=args.cpu,
    )
    if isinstance(return_code, list):
        # Should be a list of FASTA filenames
        return 0
    else:
        return return_code


def classify(args=None):
    """Subcommand to classify ITS1 sequences using a database."""
    from .classify import main

    if args.output:
        check_output_directory(args.output)
    if args.temp:
        check_output_directory(args.temp)
    return_code = main(
        fasta=args.input,
        db_url=expand_database_argument(args.database, exist=True, hyphen_default=True),
        method=args.method,
        out_dir=args.output,
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
        level=args.level,
        known=args.known,
        method=args.method,
        min_abundance=args.abundance,
        assess_output=args.output,
        map_output=args.table,
        confusion_output=args.confusion,
        debug=args.verbose,
    )


def read_summary(args=None):
    """Subcommand to run per-output-folder summary at sequence level."""
    from .read_summary import main

    if args.metadata:
        check_input_file(args.metadata)
        if not args.metacols:
            sys.exit("ERROR: Must also supply -c / --metacols argument.")
    return main(
        inputs=args.input,
        output=args.output,
        excel=args.excel,
        method=args.method,
        min_abundance=args.abundance,
        metadata_file=args.metadata,
        metadata_cols=args.metacols,
        metadata_groups=args.metagroups,
        metadata_fieldnames=args.metafields,
        metadata_index=args.metaindex,
        debug=args.verbose,
    )


def sample_summary(args=None):
    """Subcommand to run multiple-output-folder summary at sample level."""
    from .sample_summary import main

    if args.metadata:
        check_input_file(args.metadata)
        if not args.metacols:
            sys.exit("ERROR: Must also supply -c / --metacols argument.")
    return main(
        inputs=args.input,
        output=args.output,
        human_output=args.human,
        method=args.method,
        min_abundance=args.abundance,
        metadata_file=args.metadata,
        metadata_cols=args.metacols,
        # metadata_groups=args.metagroups,
        metadata_fieldnames=args.metafields,
        metadata_index=args.metaindex,
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
        min_abundance=args.abundance,
        total_min_abundance=args.total,
        always_show_db=args.showdb,
        max_edit_dist=args.editdist,
        debug=args.verbose,
    )


def pipeline(args=None):
    """Subcommand to run the default classification pipeline."""
    from .prepare import main as prepare
    from .classify import main as classify
    from .sample_summary import main as sample_summary
    from .read_summary import main as read_summary
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

    fasta_files = prepare(
        fastq=args.input,
        negative_controls=args.negctrls,
        out_dir=intermediate_dir,
        primer_dir=None,
        left_primer=ARG_PRIMER_LEFT["default"],
        right_primer=ARG_PRIMER_RIGHT["default"],
        min_abundance=args.abundance,
        tmp_dir=args.temp,
        debug=args.verbose,
        cpu=args.cpu,
    )
    if isinstance(fasta_files, int):
        return_code = fasta_files
        if return_code:
            sys.stderr.write("ERROR: Pipeline aborted during prepare-reads\n")
            sys.exit(return_code)

    classified_files = classify(
        fasta=fasta_files,
        db_url=db,
        method=args.method,
        out_dir=intermediate_dir,
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
            "ERROR: %i FASTA files but %i classified"
            % (len(fasta_files), len(classified_files))
        )

    if args.report:
        stem = os.path.join(args.output, args.report)
    else:
        # Include version number here?
        stem = os.path.join(args.output, "thapbi-pict")

    return_code = sample_summary(
        inputs=classified_files,
        output=stem + ".samples.tsv",
        human_output=stem + ".samples.txt",
        method=args.method,
        min_abundance=args.abundance,
        metadata_file=args.metadata,
        metadata_cols=args.metacols,
        # metadata_groups=args.metagroups
        metadata_fieldnames=args.metafields,
        metadata_index=args.metaindex,
        debug=args.verbose,
    )
    if return_code:
        sys.stderr.write("ERROR: Pipeline aborted during sample-summary\n")
        sys.exit(return_code)
    sys.stderr.write("Wrote %s.samples.*\n" % stem)

    return_code = read_summary(
        inputs=fasta_files + classified_files,
        output=stem + ".reads.tsv",
        excel=stem + ".reads.xlsx",
        method=args.method,
        min_abundance=args.abundance,
        metadata_file=args.metadata,
        metadata_cols=args.metacols,
        metadata_groups=args.metagroups,
        metadata_fieldnames=args.metafields,
        metadata_index=args.metaindex,
        debug=args.verbose,
    )
    if return_code:
        sys.stderr.write("ERROR: Pipeline aborted during read-summary\n")
        sys.exit(return_code)
    sys.stderr.write("Wrote %s.reads.*\n" % stem)

    edit_graph_filename = stem + ".edit-graph.xgmml"
    if os.path.isfile(edit_graph_filename):
        # This is slow to compute on complex sample sets
        sys.stderr.write(
            "WARNING: Skipping %s as already exists\n" % edit_graph_filename
        )
    else:
        # The XGMML output has minimal dependencies compared to PDF output
        return_code = edit_graph(
            graph_output=edit_graph_filename,
            graph_format="xgmml",
            db_url=db,
            inputs=fasta_files,
            min_abundance=args.abundance,
            # total_min_abundance=args.total,
            # always_show_db=args.showdb,
            # max_edit_dist=args.editdist,
            debug=args.verbose,
        )
        if return_code:
            sys.stderr.write("ERROR: Pipeline aborted during edit-graph\n")
            sys.exit(return_code)
        sys.stderr.write("Wrote %s\n" % edit_graph_filename)

    sys.stderr.write("All done!")


# Common arguments
# ================

# "-d", "--database",
ARG_DB_INPUT = dict(  # noqa: C408
    type=str, default="-", help="ITS1 database to use, default '-' is bundled database."
)

# "-i", "--input",
ARG_INPUT_FASTA = dict(  # noqa: C408
    type=str,
    required=True,
    nargs="+",
    help="One or more ITS1 FASTA filenames or folder names "
    "(containing files named *.fasta).",
)

# "-m", "--method",
ARG_METHOD_OUTPUT = dict(  # noqa: C408
    type=str,
    default=DEFAULT_METHOD,
    choices=list(method_classifier),
    help="Classify method to run, default is '%s'." % DEFAULT_METHOD,
)

# "-m", "--method",
ARG_METHOD_INPUT = dict(  # noqa: C408
    type=str,
    default=DEFAULT_METHOD,
    choices=list(method_classifier),
    help="Classify method (to infer filenames), default '%s'." % DEFAULT_METHOD,
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

# "-l", "--left",
ARG_PRIMER_LEFT = dict(  # noqa: C408
    type=str,
    default="GAAGGTGAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTA",
    metavar="PRIMER",
    help="Left primer sequence, finds and removes from start of "
    "merged read pairs. Can use IUPAC ambiguity codes. "
    "Default 21bp ITS6 'GAAGGTGAAGTCGTAACAAGG' from Cooke "
    "et al. 2000 https://doi.org/10.1006/fgbi.2000.1202 and "
    "conserved 32bp 'TTTCCGTAGGTGAACCTGCGGAAGGATCATTA'.",
)

# "-r", "--right",
ARG_PRIMER_RIGHT = dict(  # noqa: C408
    type=str,
    default="GCARRGACTTTCGTCCCYRC",
    metavar="PRIMER",
    help="Right primer sequence, finds and removes reverse "
    "complement from end of merged read pairs. Can use "
    "IUPAC ambiguity codes. Default 20bp 5.8S-1R primer "
    "'GCARRGACTTTCGTCCCYRC' from Scibetta et al. 2012 "
    "https://doi.org/10.1016/j.mimet.2011.12.012 - meaning "
    "looks for 'GYRGGGACGAAAGTCYYTGC' in merged reads.",
)

# Common pipeline arguments
# =========================

# "-i", "--input",
ARG_INPUT_FASTQ = dict(  # noqa: C408
    type=str,
    required=True,
    nargs="+",
    help="One or more ITS1 FASTQ filenames or folder names "
    "(containing files named *.fastq or *.fastq.gz).",
)

# "-n", "--negctrls",
ARG_CONTROLS = dict(  # noqa: C408
    type=str,
    nargs="+",
    help="One or more negative control FASTQ filenames or folder "
    "names (which can be duplicated in the FASTQ argument). "
    "ITS1 levels in these paired reads are used to increase "
    "the minimum abundance threshold automatically.",
)

# "-a", "--abundance",
ARG_FASTQ_MIN_ABUNDANCE = dict(  # noqa: C408
    type=int,
    default=str(DEFAULT_MIN_ABUNDANCE),
    help="Mininum abundance applied to the unique ITS1 sequences "
    "in each sample (i.e. each FASTQ pair), default %i. "
    "This may be increased based on any FASTQ controls." % DEFAULT_MIN_ABUNDANCE,
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
    default="0",
    metavar="COL",
    help="If using metadata, which column contains the sequenced sample "
    "names. Default is the first column requested as metadata output "
    "with the -c / --metacols argument. This column can contain multiple "
    "semi-colon separated names catering to the fact that a field sample "
    "could be sequenced multiple times with technical replicates.",
)

# "-g", "--metagroups",
ARG_METAGROUPS = dict(  # noqa: C408
    type=int,
    default="0",
    metavar="COL",
    help="If using metadata, which column values should be used for applying "
    "background color bands.  All samples with the same metadata value must "
    "be grouped together after sorting, as the colors are reused. "
    "Default is the first column requested as metadata output (which would "
    "also be the primary sorting key, and thus ensures all members of the "
    "same group will be together).",
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
            "THAPBI Phytophthora ITS1 Classifier Tool (PICT), v%s." % __version__
        ),
        epilog="e.g. run 'thapbi_pict pipeline -h' for the pipeline subcommand help. "
        "Please see https://thapbi-pict.readthedocs.io/ for documentation.",
    )
    parser.add_argument(
        "-v", "--version", action="version", version="THAPBI PICT v%s" % __version__
    )
    subparsers = parser.add_subparsers(
        title="subcommands", help="Each subcommand has its own additional help."
    )

    # pipeline (listing first as likely to be the most used subcommand)
    parser_pipeline = subparsers.add_parser(
        "pipeline",
        description="Run default classification pipeline on FASTQ files.",
        epilog="This is equivalent to running the individual stages (prepare-reads, "
        "classify, sample-summary, read-summary, edit-graph) with their defaults, "
        "with only a minority of settings available here.",
    )
    parser_pipeline.add_argument("-i", "--input", **ARG_INPUT_FASTQ)
    parser_pipeline.add_argument("-n", "--negctrls", **ARG_CONTROLS)
    parser_pipeline.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        metavar="DIRNAME",
        help="Output directory. Required.",
    )
    parser_pipeline.add_argument(
        "-r",
        "--report",
        type=str,
        metavar="STEM",
        help="Stem for generating report filenames.",
    )
    parser_pipeline.add_argument(
        "-s",
        "--sampleout",
        type=str,
        default="",
        metavar="DIRNAME",
        help="Output directory for intermediate files for each sample "
        "(FASTQ pair). Defaults to -o / --output.",
    )
    parser_pipeline.add_argument("-a", "--abundance", **ARG_FASTQ_MIN_ABUNDANCE)
    parser_pipeline.add_argument("-d", "--database", **ARG_DB_INPUT)
    parser_pipeline.add_argument("-m", "--method", **ARG_METHOD_OUTPUT)
    parser_pipeline.add_argument("-t", "--metadata", **ARG_METADATA)
    parser_pipeline.add_argument("-c", "--metacols", **ARG_METACOLS)
    parser_pipeline.add_argument("-x", "--metaindex", **ARG_METAINDEX)
    parser_pipeline.add_argument("-g", "--metagroups", **ARG_METAGROUPS)
    parser_pipeline.add_argument("-f", "--metafields", **ARG_METAFIELDS)
    # Can't use -t for --temp as already using for --metadata:
    parser_pipeline.add_argument("--temp", **ARG_TEMPDIR)
    parser_pipeline.add_argument("--cpu", **ARG_CPU)
    parser_pipeline.add_argument("-v", "--verbose", **ARG_VERBOSE)
    parser_pipeline.set_defaults(func=pipeline)
    del parser_pipeline

    # load-tax
    parser_load_tax = subparsers.add_parser(
        "load-tax", description="Load an NCBI taxonomy dump into an ITS1 database."
    )
    parser_load_tax.add_argument(
        "-t",
        "--tax",
        type=str,
        required=True,
        metavar="DIRNAME",
        help="Folder containing NCBI taxonomy files 'names.dmp' etc.",
    )
    parser_load_tax.add_argument("-d", "--database", **ARG_DB_WRITE)
    parser_load_tax.add_argument(
        "-a",
        "--ancestors",
        type=str,
        default="4776",
        help="Comma separated list of NCBI taxids at genus level or higher. "
        "Default is 4776 for Peronosporales, use 4783 for Phytophthora only.",
    ),
    parser_load_tax.add_argument("-v", "--verbose", **ARG_VERBOSE)
    parser_load_tax.set_defaults(func=load_tax)
    del parser_load_tax  # To prevent acidentally adding more

    # ncbi-import
    parser_ncbi_import = subparsers.add_parser(
        "ncbi-import",
        description="Load an NCBI format ITS1 FASTA file into a database. "
        "By default verifies species names against a pre-loaded taxonomy, "
        "non-matching entries are rejected.",
    )
    parser_ncbi_import.add_argument(
        "-i", "--input", type=str, required=True, help="One ITS1 fasta filename."
    )
    parser_ncbi_import.add_argument("-d", "--database", **ARG_DB_WRITE)
    parser_ncbi_import.add_argument("-n", "--name", **ARG_NAME)
    parser_ncbi_import.add_argument("-x", "--lax", **ARG_LAX)
    parser_ncbi_import.add_argument("-g", "--genus", **ARG_GENUS_ONLY)
    parser_ncbi_import.add_argument("-v", "--verbose", **ARG_VERBOSE)
    parser_ncbi_import.set_defaults(func=ncbi_import)
    del parser_ncbi_import  # To prevent acidentally adding more

    # seq-import
    parser_seq_import = subparsers.add_parser(
        "seq-import",
        description="Load classified sequences from one or more processed "
        "FASTA files into an ITS1 database. e.g. Using 'known' classifier "
        "results from single species culture samples, files which can also "
        "be used for running classifier assessment. "
        "By default verifies species names against a pre-loaded taxonomy, "
        "non-matching entries are rejected.",
    )
    parser_seq_import.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        nargs="+",
        help="One or more ITS1 FASTA and classifier filenames or folders "
        "(names containing files named *.fasta and *.method.tsv, where "
        "method is set via the -m / --method argument).",
    )
    parser_seq_import.add_argument(
        "-m",
        "--method",
        type=str,
        default="known",
        help="Method name used to determines TSV filenames from which the "
        "species classification will be read. Default is 'known' matching "
        "the convention used in the classifier assessment command for "
        "known trusted species assignments (i.e. positive controls).",
    )
    parser_seq_import.add_argument(
        "-a",
        "--abundance",
        type=int,
        default=str(DEFAULT_MIN_ABUNDANCE * 10),
        help="Mininum abundance to require before importing a sequence, "
        "over-and-above whatever was used to prepare the FASTA file. "
        "Default here is %i, ten times the default of %i used for the "
        "classification pipeline - be cautious what goes in your "
        "ITS1 database)." % (DEFAULT_MIN_ABUNDANCE * 10, DEFAULT_MIN_ABUNDANCE),
    )
    parser_seq_import.add_argument("-d", "--database", **ARG_DB_WRITE)
    parser_seq_import.add_argument("-n", "--name", **ARG_NAME)
    parser_seq_import.add_argument("-x", "--lax", **ARG_LAX)
    parser_seq_import.add_argument("-g", "--genus", **ARG_GENUS_ONLY)
    parser_seq_import.add_argument("-v", "--verbose", **ARG_VERBOSE)
    parser_seq_import.set_defaults(func=seq_import)
    del parser_seq_import  # To prevent acidentally adding more

    # legacy-import
    parser_legacy_import = subparsers.add_parser(
        "legacy-import",
        description="Load one of our legacy ITS1 FASTA files into a database. "
        "By default verifies species names against a pre-loaded taxonomy, "
        "non-matching entries are rejected.",
    )
    parser_legacy_import.add_argument(
        "-i", "--input", type=str, required=True, help="One ITS1 fasta filename."
    )
    parser_legacy_import.add_argument("-d", "--database", **ARG_DB_WRITE)
    parser_legacy_import.add_argument("-n", "--name", **ARG_NAME)
    parser_legacy_import.add_argument("-x", "--lax", **ARG_LAX)
    parser_legacy_import.add_argument("-g", "--genus", **ARG_GENUS_ONLY)
    parser_legacy_import.add_argument("-v", "--verbose", **ARG_VERBOSE)
    parser_legacy_import.set_defaults(func=legacy_import)
    del parser_legacy_import  # To prevent acidentally adding more

    # dump
    parser_dump = subparsers.add_parser(
        "dump",
        description="Export an ITS1 database to a text file.",
        epilog="e.g. 'thapbi_pict dump -d ... -c 8a,8b -o clade_8a_8b.txt'",
    )
    parser_dump.add_argument("-d", "--database", **ARG_DB_INPUT)
    parser_dump.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        help="File to write to (default '-' meaning stdout)",
    )
    parser_dump.add_argument(
        "-f",
        "--format",
        type=str,
        default="txt",
        choices=["txt", "fasta"],
        help="Format to write out (default 'txt' for debugging).",
    )
    parser_dump.add_argument(
        "-c",
        "--clade",
        type=str,
        default="",
        help="Which clade(s) to export (comma separated list, "
        "with '-' meaning no clade defined). "
        "Default is not to filter by clade.",
    )
    parser_dump.add_argument(
        "-g",
        "--genus",
        type=str,
        default="",
        help="Which genus (or genera) to export (comma separated list). "
        "Default is not to filter by genus.",
    )
    parser_dump.add_argument(
        "-s",
        "--species",
        type=str,
        default="",
        help="Which species to export (comma separated list, ignores any "
        "spaces after each comma). Requires a single genus argument be given. "
        "Default is not to filter by species.",
    )
    parser_dump.add_argument("-v", "--verbose", **ARG_VERBOSE)
    parser_dump.set_defaults(func=dump)
    del parser_dump  # To prevent acidentally adding more

    # prepare reads
    parser_prepare_reads = subparsers.add_parser(
        "prepare-reads",
        description="Trim and merge paired FASTQ files of ITS1 amplicons.",
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
    parser_prepare_reads.add_argument("-i", "--input", **ARG_INPUT_FASTQ)
    parser_prepare_reads.add_argument("-n", "--negctrls", **ARG_CONTROLS)
    parser_prepare_reads.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="DIRNAME",
        help="Directory to write output FASTA files to, "
        "default is next to each input file.",
    )
    parser_prepare_reads.add_argument("-a", "--abundance", **ARG_FASTQ_MIN_ABUNDANCE)
    parser_prepare_reads.add_argument(
        "-p",
        "--primers",
        type=str,
        default="",
        metavar="DIRNAME",
        help="Where to write optional failed primer FASTA files.",
    )
    parser_prepare_reads.add_argument("-l", "--left", **ARG_PRIMER_LEFT)
    parser_prepare_reads.add_argument("-r", "--right", **ARG_PRIMER_RIGHT)
    parser_prepare_reads.add_argument("-t", "--temp", **ARG_TEMPDIR)
    parser_prepare_reads.add_argument("-v", "--verbose", **ARG_VERBOSE)
    parser_prepare_reads.add_argument("--cpu", **ARG_CPU)
    parser_prepare_reads.set_defaults(func=prepare_reads)
    del parser_prepare_reads  # To prevent acidentally adding more

    # classify
    parser_classify = subparsers.add_parser(
        "classify",
        description="Classify FASTA file of ITS1 sequences by species.",
        epilog="Each input file XXX.fasta will result in an output file "
        "named XXX.method.tsv in the specified output directory (default "
        "input dir).",
    )
    parser_classify.add_argument("-i", "--input", **ARG_INPUT_FASTA)
    parser_classify.add_argument("-d", "--database", **ARG_DB_INPUT)
    parser_classify.add_argument("-m", "--method", **ARG_METHOD_OUTPUT)
    parser_classify.add_argument(
        "-o",
        "--output",
        type=str,
        default="",
        metavar="DIRNAME",
        help="Directory to write output reports to, default (empty "
        "string) is next to each input file. Use '-' for stdout.",
    )
    parser_classify.add_argument("-t", "--temp", **ARG_TEMPDIR)
    parser_classify.add_argument("-v", "--verbose", **ARG_VERBOSE)
    parser_classify.add_argument("--cpu", **ARG_CPU)
    parser_classify.set_defaults(func=classify)
    del parser_classify  # To prevent acidentally adding more

    # assess-classification
    parser_assess = subparsers.add_parser(
        "assess",
        description="Assess accuracy of ITS1 read classification.",
        epilog="Takes as input predictions named XXX.method.tsv "
        "and matching expected classifications in XXX.known.tsv "
        "(which can be in different directories) to produce a "
        "multi-species confusion matrix (output on request) and "
        "classifier performance metrics (to stdout by default). "
        "You can deliberately compare two prediction methods to "
        "each other using this, but a known set of positive controls "
        "is the expected benchmark.",
    )
    parser_assess.add_argument(
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
    parser_assess.add_argument(
        "-l",
        "--level",
        type=str,
        default="sample",
        choices=["sample", "sseq", "useq"],
        help="Assess at sample level (taking union of species predicted "
        "by sequences from each sample), sequence level within samples, "
        "or at unique sequence level (over all samples).",
    )
    parser_assess.add_argument(
        "-k",
        "--known",
        type=str,
        default="known",
        help="Replaces the string used in filenames for the truth "
        "against which the method in -m / --method is assessed. "
        "This could be any defined method, default is 'known'.",
    )
    parser_assess.add_argument("-m", "--method", **ARG_METHOD_INPUT)
    parser_assess.add_argument(
        "-a",
        "--abundance",
        type=int,
        default="1",
        help="Mininum abundance to require before considering a classification. "
        "Default is one meaning look at everything, but rather than re-running "
        "the classifier with a stricter minimum abundance you can apply it here.",
    )
    parser_assess.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="FILENAME",
        help="File to write species level classification assessment table to. "
        "Default is '-' meaning to stdout.",
    )
    parser_assess.add_argument(
        "-t",
        "--table",
        type=str,
        metavar="FILENAME",
        help="File to write expected-to-predicted mapping tally table to. "
        "Can use '-' meaning to stdout. Default is not to write this file.",
    )
    parser_assess.add_argument(
        "-c",
        "--confusion",
        type=str,
        metavar="FILENAME",
        help="File to write species level confusion matrix to. "
        "Can use '-' meaning to stdout. Default is not to write this file.",
    )
    parser_assess.add_argument("-v", "--verbose", **ARG_VERBOSE)
    parser_assess.set_defaults(func=assess_classification)
    del parser_assess  # To prevent acidentally adding more

    # read-summary
    parser_read_summary = subparsers.add_parser(
        "read-summary",
        description="Sequence-level summary report on classifier output.",
        epilog="Assumes you've run prepare-reads and classify, and have "
        "folders with XXX.fasta and XXX.method.tsv files from your plate(s). "
        "The output is a table with one row per unique sequence (as trimmed "
        "by the prepare-reads step, can be thousands of rows) and one column "
        "per sample (typically 96 samples).",
    )
    parser_read_summary.add_argument(
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
    parser_read_summary.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="FILENAME",
        help="File to write summary sequence vs samples TSV table to. "
        "Default is '-' meaning to stdout.",
    )
    parser_read_summary.add_argument(
        "-e",
        "--excel",
        type=str,
        default="",
        metavar="FILENAME",
        help="File to write summary sequence vs samples Excel table to. "
        "Default is '' meaning no output.",
    )
    parser_read_summary.add_argument(
        "-m",
        "--method",
        type=str,
        default=DEFAULT_METHOD,
        help="Classifier method(s) to report, comma separaed list (used to infer "
        "filenames), default is '%s' (only)." % DEFAULT_METHOD,
    )
    parser_read_summary.add_argument(
        "-a",
        "--abundance",
        type=int,
        default=str(DEFAULT_MIN_ABUNDANCE),
        help="Mininum sample level abundance to require for the report. "
        "Default %i reflects default in prepare-reads. Rather than re-running "
        "the prepare or classifier steps with a stricter minimum abundance you "
        "can apply it here. Use zero or one to look at everything (but beware "
        "that negative control samples will include low abundance entries)."
        % DEFAULT_MIN_ABUNDANCE,
    )
    parser_read_summary.add_argument("-t", "--metadata", **ARG_METADATA)
    parser_read_summary.add_argument("-c", "--metacols", **ARG_METACOLS)
    parser_read_summary.add_argument("-x", "--metaindex", **ARG_METAINDEX)
    parser_read_summary.add_argument("-g", "--metagroups", **ARG_METAGROUPS)
    parser_read_summary.add_argument("-f", "--metafields", **ARG_METAFIELDS)
    parser_read_summary.add_argument("-v", "--verbose", **ARG_VERBOSE)
    parser_read_summary.set_defaults(func=read_summary)
    del parser_read_summary  # To prevent acidentally adding more

    # sample-summary
    parser_sample_summary = subparsers.add_parser(
        "sample-summary",
        description="Sample-level summary report on classifier output.",
        epilog="Assumes you've run prepare-reads and classify, and have "
        "folders with XXX.method.tsv files from your samples. The output "
        "is a table with rows for each sample (XXX), describing the species "
        "predicted to be present, and the associated sequence count. "
        "Intended to be used for samples over multiple sequencing plates.",
    )
    parser_sample_summary.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        nargs="+",
        help="One or more prediction files (*.method.tsv) or folder names. "
        "The files should follow this naming convention, where the classifer "
        "method appearing in the extension can be set via -m / --method.",
    )
    parser_sample_summary.add_argument("-m", "--method", **ARG_METHOD_INPUT)
    parser_sample_summary.add_argument(
        "-a",
        "--abundance",
        type=int,
        default=str(DEFAULT_MIN_ABUNDANCE),
        help="Mininum sample level abundance to require for the report. "
        "Default %i reflects default in prepare-reads. Rather than re-running "
        "the prepare or classifier steps with a stricter minimum abundance you "
        "can apply it here. Use zero or one to look at everything (but beware "
        "that negative control samples will include low abundance entries)."
        % DEFAULT_MIN_ABUNDANCE,
    )
    parser_sample_summary.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="FILENAME",
        help="File to write sample species classification summary table to. "
        "Default is '-' meaning to stdout.",
    )
    parser_sample_summary.add_argument(
        "-r",
        "--human",
        type=str,
        metavar="FILENAME",
        help="File to write human readable smaple level species predictions to. "
        "Can use '-' meaning to stdout. Default is not to write this file.",
    )
    parser_sample_summary.add_argument("-t", "--metadata", **ARG_METADATA)
    parser_sample_summary.add_argument("-c", "--metacols", **ARG_METACOLS)
    parser_sample_summary.add_argument("-x", "--metaindex", **ARG_METAINDEX)
    parser_sample_summary.add_argument("-f", "--metafields", **ARG_METAFIELDS)
    parser_sample_summary.add_argument("-v", "--verbose", **ARG_VERBOSE)
    parser_sample_summary.set_defaults(func=sample_summary)
    del parser_sample_summary  # To prevent acidentally adding more

    # edit-graph
    parser_edit_graph = subparsers.add_parser(
        "edit-graph",
        description="Draw network graph of sequences using edit distance.",
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
    arg = parser_edit_graph.add_argument("-d", "--database", **ARG_DB_INPUT)
    arg.help += " Use '' to mean no database."
    del arg
    # Currently ARG_INPUT_FASTA uses required=True, but we need to change thant:
    arg = parser_edit_graph.add_argument("-i", "--input", **ARG_INPUT_FASTA)
    arg.required = False
    del arg
    parser_edit_graph.add_argument(
        "-a",
        "--abundance",
        type=int,
        default=str(DEFAULT_MIN_ABUNDANCE),
        help="Mininum sample level abundance for FASTA sequences. "
        "Default %i reflects default in prepare-reads." % DEFAULT_MIN_ABUNDANCE,
    )
    parser_edit_graph.add_argument(
        "-t",
        "--total",
        type=int,
        default="0",
        help="Mininum total abundance for FASTA sequences. "
        "Applied after per-sample level minimum (-a / --abundance). "
        "Offered as a way to simplify the final graph.",
    )
    parser_edit_graph.add_argument(
        "-s",
        "--showdb",
        action="store_true",
        help="Show DB entries, regardless of their abundance in the "
        "FASTA inputs. Required if only drawing DB entries.",
    )
    parser_edit_graph.add_argument(
        "-e",
        "--editdist",
        type=int,
        default="3",
        choices=[1, 2, 3],
        help="Maximum edit distance to draw. Default and maximum 3.",
    )
    parser_edit_graph.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="FILENAME",
        help="Write graph here. Default is '-' meaning stdout.",
    )
    parser_edit_graph.add_argument(
        "-f",
        "--format",
        type=str,
        default="xgmml",
        choices=["graphml", "gexf", "gml", "xgmml", "pdf"],
        help="Format to write out (default 'xgmml' for use with Cytoscape).",
    )
    parser_edit_graph.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose logging"
    )
    parser_edit_graph.set_defaults(func=edit_graph)
    del parser_edit_graph  # To prevent accidentally adding more

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
