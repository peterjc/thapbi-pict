# Copyright 2018-2023 by Peter Cock, The James Hutton Institute.
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
from .db_import import DEF_MAX_LENGTH
from .db_import import DEF_MIN_LENGTH
from .db_import import fasta_parsing_function as fasta_conventions
from .db_orm import connect_to_db
from .db_orm import MarkerDef
from .prepare import find_fastq_pairs
from .utils import file_to_sample_name
from .utils import find_requested_files

# Common command line defaults
# ============================


DEFAULT_METHOD = "onebp"
DEFAULT_MIN_ABUNDANCE = 100
CTRL_A = chr(1)
IGNORE_PREFIXES = ("Undetermined_", "unknown_")

# Argument validation functions
# =============================


def check_cpu(cpu):
    """If zero use all available CPUs, otherwise cap the given value."""
    try:
        # Don't need to check $SLURM_CPUS_PER_TASK,
        # probably don't need to check $NSLOTS on SGE either
        available = len(os.sched_getaffinity(0))
    except AttributeError:
        # Unavailable on macOS or Windows
        available = os.cpu_count()
    if not cpu:
        return available
    elif cpu > 0:
        return min(available, cpu)  # cap requested number
    else:
        sys.exit(
            f"ERROR: CPU argument {cpu} should be positive, or zero for all available."
        )


def check_input_file(filename):
    """Command line validation of an input filename."""
    if not os.path.isfile(filename):
        sys.exit(f"ERROR: Could not find input file: {filename}")


def check_output_stem(out_stem, dir_only_ok=False, dir_must_exist=True):
    """Command line validation of output stem value.

    Returns the output directory, or aborts.
    """
    if not out_stem:
        sys.exit("ERROR: Output stem is blank")
    elif out_stem.endswith(os.path.sep) or os.path.isdir(out_stem):
        if dir_only_ok:
            return out_stem
        sys.exit(f"ERROR: Output stem {out_stem!r} is a directory")
    out_dir, basename = os.path.split(out_stem)
    if not basename:
        sys.exit(f"ERROR: Output stem needs a partial filename: {out_stem!r}")
    elif not out_dir or os.path.isdir(out_dir):
        return out_dir
    elif os.path.isfile(out_dir):
        sys.exit(f"ERROR: Output stem directory name is a file: {out_dir!r}")
    elif dir_must_exist:
        sys.exit(f"ERROR: Output stem directory does not exist: {out_dir!r}")
    return None


def check_output_directory(out_dir, must_exist=True):
    """Command line validation of output directory value."""
    if not out_dir:
        sys.exit("ERROR: Output directory name blank")
    elif out_dir == "-" or os.path.isdir(out_dir):
        return True
    elif os.path.isfile(out_dir):
        sys.exit(f"ERROR: Output directory name is a file: {out_dir!r}")
    elif must_exist:
        sys.exit(f"ERROR: Output directory does not exist: {out_dir!r}")
    return None


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


def db_import(args=None):
    """Subcommand to import a FASTA file into a basebase."""
    from .db_import import main

    return main(
        fasta=args.input,
        db_url=expand_database_argument(args.database),
        marker=args.marker,
        left_primer=args.left,
        right_primer=args.right,
        min_length=args.minlen,
        max_length=args.maxlen,
        name=args.name,
        convention=args.convention,
        sep=args.sep,
        validate_species=not args.lax,
        genus_only=args.genus,
        ignore_prefixes=tuple(args.ignore_prefixes),
        tmp_dir=args.temp,
        debug=args.verbose,
    )


def dump(args=None):
    """Subcommand to dump a database to a text file."""
    from .dump import main

    return main(
        db_url=expand_database_argument(args.database, exist=True, hyphen_default=True),
        output_filename=args.output,
        output_format=args.format,
        marker=args.marker,
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

    # Connect to the DB,
    db = expand_database_argument(args.database, exist=True, hyphen_default=True)
    Session = connect_to_db(db)
    session = Session()

    return_code = main(
        fastq=args.input,
        out_dir=args.output,
        session=session,
        flip=args.flip,
        min_abundance=args.abundance,
        min_abundance_fraction=args.abundance_fraction,
        ignore_prefixes=tuple(args.ignore_prefixes),
        merged_cache=args.merged_cache,
        tmp_dir=args.temp,
        debug=args.verbose,
        cpu=check_cpu(args.cpu),
    )

    session.close()
    if isinstance(return_code, int):
        return return_code
    else:
        # Should be lists of FASTA filenames
        return 0


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


def denoise(args=None):
    """Subcommand to denoise FASTA file(s) using UNOISE read-correction."""
    from .denoise import main

    return main(
        inputs=args.input,
        output=args.output,
        denoise_algorithm=args.denoise,
        total_min_abundance=args.total,
        min_length=args.minlen,
        max_length=args.maxlen,
        unoise_alpha=args.unoise_alpha,
        unoise_gamma=args.unoise_gamma,
        tmp_dir=args.temp,
        debug=args.verbose,
        cpu=check_cpu(args.cpu),
    )


def sample_tally(args=None):
    """Subcommand to tally per-sample FASTA files using MD5 naming."""
    from .sample_tally import main

    # Connect to the DB,
    db = expand_database_argument(args.database, exist=True, hyphen_default=True)
    Session = connect_to_db(db)
    session = Session()

    return main(
        inputs=args.input,
        synthetic_controls=args.synctrls,
        negative_controls=args.negctrls,
        output=args.output,
        session=session,
        marker=args.marker,
        spike_genus=args.synthetic,
        fasta=args.fasta,
        min_abundance=args.abundance,
        min_abundance_fraction=args.abundance_fraction,
        total_min_abundance=args.total,
        min_length=args.minlen,
        max_length=args.maxlen,
        unoise_alpha=args.unoise_alpha,
        unoise_gamma=args.unoise_gamma,
        denoise_algorithm=args.denoise,
        tmp_dir=args.temp,
        debug=args.verbose,
        cpu=check_cpu(args.cpu),
    )


def classify(args=None):
    """Subcommand to classify FASTA sequences using a database."""
    from .classify import main

    if args.output:
        check_output_directory(args.output, must_exist=False)
    if args.temp:
        check_output_directory(args.temp)

    # Connect to the DB,
    Session = connect_to_db(
        expand_database_argument(args.database, exist=True, hyphen_default=True)
    )
    session = Session()

    return_code = main(
        inputs=args.input,
        session=session,
        marker_name=args.marker,
        method=args.method,
        out_dir=args.output,
        ignore_prefixes=tuple(args.ignore_prefixes),
        min_abundance=args.abundance,
        tmp_dir=args.temp,
        biom=args.biom,
        debug=args.verbose,
        cpu=check_cpu(args.cpu),
    )

    session.close()
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
        marker=args.marker,
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

    check_output_stem(args.output, dir_only_ok=False)

    if args.metadata:
        check_input_file(args.metadata)
        if not args.metacols:
            sys.exit("ERROR: Must also supply -c / --metacols argument.")

    return main(
        inputs=args.input,
        report_stem=args.output,
        method=args.method,
        min_abundance=args.abundance,
        metadata_file=args.metadata,
        metadata_encoding=args.metaencoding,
        metadata_cols=args.metacols,
        metadata_groups=args.metagroups,
        metadata_fieldnames=args.metafields,
        metadata_index=args.metaindex,
        require_metadata=args.requiremeta,
        show_unsequenced=args.unsequenced,
        ignore_prefixes=tuple(args.ignore_prefixes),
        biom=args.biom,
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
        input_file=args.input,
        min_abundance=args.abundance,
        total_min_abundance=args.total,
        min_samples=args.min_samples,
        show_db_marker=args.marker,
        max_edit_dist=args.editdist,
        ignore_prefixes=tuple(args.ignore_prefixes),
        debug=args.verbose,
    )


def pipeline(args=None):
    """Subcommand to run the default classification pipeline."""
    from .prepare import main as prepare
    from .sample_tally import main as sample_tally
    from .classify import main as classify
    from .summary import main as summary
    from .assess import main as assess

    check_output_stem(args.output, dir_only_ok=True)
    if args.temp:
        check_output_directory(args.temp)
    if args.sampleout:
        check_output_directory(args.sampleout)
        intermediate_dir = args.sampleout
    else:
        intermediate_dir = ""
    if args.metadata:
        check_input_file(args.metadata)
        if not args.metacols:
            sys.exit("ERROR: Must also supply -c / --metacols argument.")
    method = args.method

    # Connect to the DB,
    db = expand_database_argument(args.database, exist=True, hyphen_default=True)
    Session = connect_to_db(db)
    session = Session()
    markers = sorted(_.name for _ in session.query(MarkerDef))

    # TODO - apply require_metadata=True to the prepare and classify steps?

    # We will silently ignore controls not also in the inputs argument!
    neg_stems = []
    if args.negctrls:
        neg_stems = [
            os.path.split(stem)[1]
            for (stem, R1, R2) in find_fastq_pairs(
                args.negctrls, ignore_prefixes=tuple(args.ignore_prefixes)
            )
        ]
    syn_stems = []
    if args.synctrls:
        syn_stems = [
            os.path.split(stem)[1]
            for (stem, R1, R2) in find_fastq_pairs(
                args.synctrls, ignore_prefixes=tuple(args.ignore_prefixes)
            )
        ]

    # This will do all the markers itself
    return_code = prepare(
        fastq=args.input,
        out_dir=intermediate_dir,
        session=session,
        flip=args.flip,
        # Unless pipeline explicitly using min abundance 0 or 1
        # (which are equivalent ways to retain singletons),
        # at this point only want to exclude singletons:
        min_abundance=min(2, args.abundance),
        min_abundance_fraction=0.0,
        ignore_prefixes=tuple(args.ignore_prefixes),
        merged_cache=args.merged_cache,
        tmp_dir=args.temp,
        debug=args.verbose,
        cpu=check_cpu(args.cpu),
    )
    if not isinstance(return_code, list):
        session.close()
        sys.stderr.write("ERROR: Pipeline aborted during prepare-reads\n")
        sys.exit(return_code)
    # If not an integer, should be a list of filenames:
    all_fasta_files = return_code
    # TODO - Support known setting...
    # TODO - Can we specify different expected results from diff markers?
    known_files = find_requested_files(
        args.input,
        ".known.tsv",
        ignore_prefixes=tuple(args.ignore_prefixes),
        debug=args.verbose,
    )

    all_tally_files = []
    all_classified_files = []
    for marker in markers:
        if args.output.endswith(os.path.sep) or os.path.isdir(args.output):
            # Just a directory
            stem = os.path.join(args.output, marker)
        else:
            # Have a filename stem (possibly with a directory)
            stem = f"{args.output}.{marker}"
        fasta_files = [
            _
            for _ in all_fasta_files
            if _.startswith(os.path.join(intermediate_dir, marker) + os.path.sep)
        ]
        if len(markers) > 1:
            sys.stderr.write("\n")
            sys.stderr.write(f"Processesing {marker}\n")
            sys.stderr.write("~" * (13 + len(marker)) + "\n")
            sys.stderr.write("\n")
        tally_seqs_file = f"{stem}.tally.tsv"
        sample_tally(
            inputs=fasta_files,
            synthetic_controls=[
                _
                for _ in fasta_files
                if file_to_sample_name(os.path.split(_)[1]) in syn_stems
            ],
            negative_controls=[
                _
                for _ in fasta_files
                if file_to_sample_name(os.path.split(_)[1]) in neg_stems
            ],
            output=tally_seqs_file,
            session=session,
            marker=marker,
            spike_genus=args.synthetic,
            min_abundance=args.abundance,
            min_abundance_fraction=args.abundance_fraction,
            # Historical behaviour, discards rare control-only ASVs:
            total_min_abundance=args.abundance,
            # min_length=args.minlen,
            # max_length=args.maxlen,
            denoise_algorithm=args.denoise,
            unoise_alpha=args.unoise_alpha,
            unoise_gamma=args.unoise_gamma,
            tmp_dir=args.temp,
            debug=args.verbose,
            cpu=check_cpu(args.cpu),
        )
        all_tally_files.append(tally_seqs_file)
        classified_files = classify(
            inputs=[tally_seqs_file],
            session=session,
            marker_name=marker,
            method=args.method,
            out_dir=os.path.split(stem)[0],  # i.e. next to FASTA all_reads input,
            ignore_prefixes=tuple(args.ignore_prefixes),  # not really needed
            min_abundance=args.abundance,
            biom=args.biom,
            tmp_dir=args.temp,
            debug=args.verbose,
            cpu=check_cpu(args.cpu),
        )
        if isinstance(classified_files, int):
            return_code = classified_files
            if return_code:
                sys.stderr.write(f"ERROR: Pipeline aborted during {marker} classify\n")
                sys.exit(return_code)
        if 1 != len(classified_files):
            sys.exit(
                f"ERROR: {len(classified_files)} classifier output files, expected one"
            )
        all_classified_files.extend(classified_files)

        return_code = summary(
            inputs=[tally_seqs_file, *classified_files],
            report_stem=stem,
            method=args.method,
            min_abundance=args.abundance,
            metadata_file=args.metadata,
            metadata_encoding=args.metaencoding,
            metadata_cols=args.metacols,
            metadata_groups=args.metagroups,
            metadata_fieldnames=args.metafields,
            metadata_index=args.metaindex,
            require_metadata=args.requiremeta,
            show_unsequenced=args.unsequenced,
            ignore_prefixes=tuple(args.ignore_prefixes),
            biom=args.biom,
            debug=args.verbose,
        )
        if return_code:
            sys.stderr.write(f"ERROR: Pipeline aborted during {marker} summary\n")
            session.close()
            sys.exit(return_code)

        if known_files:
            sys.stderr.write(f"Assessing {marker} classification...\n")
            return_code = assess(
                inputs=known_files + classified_files,
                known="known",  # =args.known,
                db_url=db,
                marker=marker,
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
                session.close()
                sys.exit(return_code)
            sys.stderr.write(f"Wrote {stem}.assess.*.{method}.*\n")

    if len(markers) > 1:
        # Pooled marker report
        if len(markers) > 1:
            sys.stderr.write("\n")
            sys.stderr.write("Processesing pooled markers\n")
            sys.stderr.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
            sys.stderr.write("\n")
        if args.output.endswith(os.path.sep) or os.path.isdir(args.output):
            # Just a directory
            stem = os.path.join(args.output, "pooled")
        else:
            # Have a filename stem (possibly with a directory)
            stem = f"{args.output}.pooled"
        return_code = summary(
            inputs=all_classified_files,
            report_stem=stem,
            method=args.method,
            min_abundance=args.abundance,
            metadata_file=args.metadata,
            metadata_encoding=args.metaencoding,
            metadata_cols=args.metacols,
            metadata_groups=args.metagroups,
            metadata_fieldnames=args.metafields,
            metadata_index=args.metaindex,
            require_metadata=args.requiremeta,
            show_unsequenced=args.unsequenced,
            ignore_prefixes=tuple(args.ignore_prefixes),
            biom=args.biom,
            debug=args.verbose,
        )
        if return_code:
            sys.stderr.write("ERROR: Pipeline aborted during pooled marker summary\n")
            session.close()
            sys.exit(return_code)

        if known_files:
            sys.stderr.write("Assessing pooled classification...\n")
            return_code = assess(
                inputs=known_files + all_classified_files,
                known="known",  # =args.known,
                db_url=db,
                marker=None,  # all of them!
                method=args.method,
                min_abundance=args.abundance,
                assess_output=f"{stem}.assess.{method}.tsv",
                map_output=f"{stem}.assess.tally.{method}.tsv",
                confusion_output=f"{stem}.assess.confusion.{method}.tsv",
                ignore_prefixes=tuple(args.ignore_prefixes),
                debug=args.verbose,
            )
            if return_code:
                sys.stderr.write("ERROR: Pipeline aborted during pooled assess\n")
                session.close()
                sys.exit(return_code)
            sys.stderr.write(f"Wrote {stem}.assess.*.{method}.*\n")

    session.close()

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
        metadata_encoding=args.metaencoding,
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
ARG_INPUT_FASTA_OR_TSV = dict(  # noqa: C408
    type=str,
    required=True,
    nargs="+",
    help="One or more prepared FASTA or sample tally TSV filenames, or folder "
    "names (containing files named *.fasta or *.tsv).",
)

# "--ignore-prefixes",
ARG_IGNORE_PREFIXES = dict(  # noqa: C408
    type=str,
    nargs="+",
    metavar="STEM",
    default=IGNORE_PREFIXES,
    help="One or more filename prefixes to ignore, default "
    + " ".join(IGNORE_PREFIXES),
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

# "-b", "--biom",
ARG_BIOM = dict(action="store_true", help="Enable optional BIOM output.")  # noqa: C408

# "--merged-cache",
ARG_MERGED_CACHE = dict(  # noqa: C408
    type=str,
    required=False,
    default="",
    metavar="DIRNAME",
    help="Advanced option. Cache directory for temporary per-sample FASTA "
    "file after overlap merging (usually the slowest step), but before "
    "primer trimming and abundance threshold.",
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
    type=int,
    default=1,
    help="Number of parallel threads to use in called tools. "
    "Use zero for all available. Default 1.",
)

# Common import arguments
# =======================

# "-k", "--marker",
ARG_MARKER = dict(  # noqa: C408
    type=str,
    help="Unique marker name for a barcode, linked to a PCR primer pair for "
    "a diagnostic amplicon. e.g. ITS1, COI, 12S. If not already in the DB, "
    "it will be added but this requires both left and right primers, and min "
    "and max product length are recommended. Required unless DB has exactly "
    "one marker already defined in it.",
)
ARG_MARKER_PICK = dict(  # noqa: C408
    type=str,
    default="",
    # Comma separated?
    help="If DB has multiple amplicon markers, which one should be used?",
)

# "-d", "--database",
ARG_DB_WRITE = dict(  # noqa: C408
    type=str, required=True, help="Which database to write to (or create)."
)

# Prepare reads arguments
# =======================

# "--minlen",
ARG_MIN_LENGTH = dict(  # noqa: C408
    type=int,
    default=DEF_MIN_LENGTH,
    metavar="LENGTH",
    help=f"Minimum length sequence to accept (default {DEF_MIN_LENGTH}).",
)

# "--maxlen",
ARG_MAX_LENGTH = dict(  # noqa: C408
    type=int,
    default=DEF_MAX_LENGTH,
    metavar="LENGTH",
    help=f"Maximum length sequence to accept (default {DEF_MAX_LENGTH}).",
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

# Read-correction / denoise arguments
# ===================================

# "--denoise",
ARG_DENOISE = dict(  # noqa: C408
    choices=["-", "unoise-l", "usearch", "vsearch"],
    default="-",
    help="Optional read-correction algorithm, default '-' for none. "
    "Use 'unoise-l' for built-in reimplementation of the Levenshtein distance "
    "UNOISE algorithm as described in Edgar (2016) Use 'usearch' to "
    "call external tool 'usearch -unoise3 ...' and the original "
    "author's UNOISE3 implementation. Use 'vsearch' to call "
    "external tool 'vsearch --cluster_unoise ...' and their UNOISE3 "
    "reimplementation using pairwise alignment based distance.",
)

# "-α", "--unoise_alpha",
ARG_UNOISE_ALPHA = dict(  # noqa: C408
    type=float,
    default=None,
    metavar="FLOAT",
    help="UNOISE read-correction alpha parameter (α), used in "
    "skew threshold function beta (β). Default 2.0 for UNOISE-L, "
    "tool defaults for USEARCH and VSEARCH.",
)
# "-γ", "--unoise_gamma"
ARG_UNOISE_GAMMA = dict(  # noqa: C408
    type=int,
    default=None,
    metavar="INT",
    help="UNOISE read-correction gamma parameter (γ). Variants "
    "below this total abundance are discarded before denoising. Default 4 "
    "for UNOISE-L, tool defaults for USEARCH and VSEARCH.",
)

# Common pipeline arguments
# =========================

# "-i", "--input",
ARG_INPUT_FASTQ = dict(  # noqa: C408
    type=str,
    required=True,
    nargs="+",
    metavar="FASTQ",
    help="One or more paired FASTQ filenames or folder names "
    "(containing files named *.fastq or *.fastq.gz).",
)


# "--synthetic",
ARG_SYNTHETIC_SPIKE = dict(  # noqa: C408
    type=str,
    default="synthetic",
    metavar="GENUS",
    help="Comma separated genus list of any spike-in sequences in the "
    "negative controls, used with the database synthetic controls settings. "
    "Default 'synthetic'.",
)

# "-y", "--synctrls",
ARG_SYN_CONTROLS = dict(  # noqa: C408
    type=str,
    nargs="+",
    metavar="FASTQ",
    # Does accept folder names, but kind of pointless
    # (as would be applied only to that folder)
    help="One or more synthetic control paired FASTQ filenames (which must "
    "also appear in the inputs or they will be ignored). High non-synthetic "
    "marker levels will increase the fractional abundance threshold of other "
    "FASTQ files in the folder. Can use '-' for none.",
)
ARG_SYN_CONTROLS_FASTA = dict(  # noqa: C408
    type=str,
    nargs="+",
    metavar="FASTA",
    help="One or more synthetic control FASTA filenames (which must also "
    "appear in inputs or they will be ignored). High non-synthetic marker "
    "levels will increase the fractional abundance threshold of other FASTA "
    "files from the same threshold pool (set in FASTA header metadata). "
    "Can use '-' for none.",
)

# "-n", "--negctrls",
ARG_NEG_CONTROLS = dict(  # noqa: C408
    type=str,
    nargs="+",
    metavar="FASTQ",
    # Does accept folder names, but kind of pointless
    # (as would be applied only to that folder)
    help="One or more negative control paired FASTQ filenames (which must "
    "also appear in the inputs or they will be ignored). High non-synthetic "
    "levels will increase the absolute minimum abundance threshold of other "
    "FASTQ files in the same folder. Can use '-' for none.",
)
ARG_NEG_CONTROLS_FASTA = dict(  # noqa: C408
    type=str,
    nargs="+",
    metavar="FASTA",
    help="One or more synthetic control FASTA filenames (which must also "
    "appear in the inputs or they will be ignored). High non-synthetic marker "
    "levels will increase the absolute abundance threshold of other FASTA "
    "files from the same threshold pool (set in FASTA header metadata). "
    "Can use '-' for none.",
)

# "-a", "--abundance",
ARG_FASTQ_MIN_ABUNDANCE = dict(  # noqa: C408
    type=int,
    metavar="INTEGER",
    default=str(DEFAULT_MIN_ABUNDANCE),
    help="Minimum abundance applied to unique marker sequences in each sample"
    f" (i.e. each FASTQ pair). Default {DEFAULT_MIN_ABUNDANCE}."
    " May be increased based on negative controls."
    " Half this value is applied to synthetic controls.",
)
ARG_FASTQ_MIN_ABUNDANCE_TWO = dict(  # noqa: C408
    type=int,
    metavar="INTEGER",
    default=str(DEFAULT_MIN_ABUNDANCE),
    help="Minimum abundance applied to unique marker sequences in each sample"
    " (i.e. each FASTQ pair). Default 2 meaning only exclude singletons.",
)
# "-f", "--abundance-fraction",
ARG_FASTQ_NOISE_PERC = dict(  # noqa: C408
    type=float,
    metavar="FLOAT",
    default="0.001",
    help="Minimum abundance fraction, low frequency noise threshold applied"
    " to unique marker sequences in each sample. Default 0.001 meaning 0.1%%."
    " May be increased based on synthetic controls."
    " Half this value is applied to negative controls.",
)
ARG_FASTQ_NOISE_PERC_ZERO = dict(  # noqa: C408
    type=float,
    metavar="FLOAT",
    default="0",
    help="Minimum abundance fraction, low frequency noise threshold applied"
    " to unique marker sequences in each sample. e.g. 0.001 meaning 0.1%%."
    " Default zero, meaning not used.",
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

# "-e", "--metaencoding",
ARG_METAENCODING = dict(  # noqa: C408
    type=str,
    default="",
    metavar="ENCODING",
    help="Optional encoding of the metadata table, e.g. 'UTF-8', 'latin1', "
    "or if exporting from Microsoft Excel on macOS use 'macintosh'.",
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
    "Zero (default) is interpreted as the first column requested as metadata "
    "output (the primary sorting key, ensuring all members of the same group "
    "will be together).",
)

# "--metafields",
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
        "sample-tally, classify, summary, assess) with their defaults, with only a "
        "minority of settings available here.",
    )
    subcommand_parser.add_argument("-i", "--input", **ARG_INPUT_FASTQ)
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument("-y", "--synctrls", **ARG_SYN_CONTROLS)
    subcommand_parser.add_argument("-n", "--negctrls", **ARG_NEG_CONTROLS)
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        metavar="STEM",
        help="Output directory and/or filename stem for generating reports.",
    )
    subcommand_parser.add_argument(
        "-s",
        "--sampleout",
        type=str,
        default="",
        metavar="DIRNAME",
        help="Advanced option. Cache directory for temporary per-sample FASTA "
        "files after primer trimming.",
    )
    subcommand_parser.add_argument("-a", "--abundance", **ARG_FASTQ_MIN_ABUNDANCE)
    subcommand_parser.add_argument("-f", "--abundance-fraction", **ARG_FASTQ_NOISE_PERC)
    subcommand_parser.add_argument("-d", "--database", **ARG_DB_INPUT)
    subcommand_parser.add_argument("--synthetic", **ARG_SYNTHETIC_SPIKE)
    subcommand_parser.add_argument("--flip", **ARG_FLIP)
    subcommand_parser.add_argument("--denoise", **ARG_DENOISE)
    subcommand_parser.add_argument("-α", "--unoise_alpha", **ARG_UNOISE_ALPHA)
    subcommand_parser.add_argument("-γ", "--unoise_gamma", **ARG_UNOISE_GAMMA)
    subcommand_parser.add_argument("-m", "--method", **ARG_METHOD_OUTPUT)
    subcommand_parser.add_argument("-t", "--metadata", **ARG_METADATA)
    subcommand_parser.add_argument("-e", "--metaencoding", **ARG_METAENCODING)
    subcommand_parser.add_argument("-c", "--metacols", **ARG_METACOLS)
    subcommand_parser.add_argument("-x", "--metaindex", **ARG_METAINDEX)
    subcommand_parser.add_argument("-g", "--metagroups", **ARG_METAGROUPS)
    subcommand_parser.add_argument("--metafields", **ARG_METAFIELDS)
    subcommand_parser.add_argument("--merged-cache", **ARG_MERGED_CACHE)
    subcommand_parser.add_argument("-q", "--requiremeta", **ARG_REQUIREMETA)
    subcommand_parser.add_argument("-u", "--unsequenced", **ARG_UNSEQUENCED)
    subcommand_parser.add_argument("-b", "--biom", **ARG_BIOM)
    # Can't use -t for --temp as already using for --metadata:
    subcommand_parser.add_argument("--temp", **ARG_TEMPDIR)
    subcommand_parser.add_argument("--cpu", **ARG_CPU)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=pipeline)
    del subcommand_parser

    # load-tax
    subcommand_parser = subparsers.add_parser(
        "load-tax", description="Load an NCBI taxonomy dump into a marker database."
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
    )
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=load_tax)
    del subcommand_parser  # To prevent accidentally adding more

    # import
    subcommand_parser = subparsers.add_parser(
        "import",
        description="Load a FASTA file into a database. "
        "Supports multiple entries in a single record using a FASTA title "
        "split by the specified separator character.\n\n"
        "Will use the NCBI taxid alone if given and in the database. If not,"
        "by default verifies genus/species names against the pre-loaded "
        "taxonomy, failing that a matching genus, otherwise rejected.\n\n"
        "By default assumes FASTA title line entries are an accesion followed"
        "by an exact species name only, optionally with taxid=... to finish."
        "With NCBI heuristics enabled, tries to split the species and any "
        "following free text using the taxonomy names in the DB. For SINTAX "
        "format only the genus and species are used. For ObiTools format, the "
        "NCBI taxid alone is used if in the DB, falling back on the genus and "
        "species given. For NCBI taxonomy mode, looks for text like "
        "taxid=123456 and takes genus/species from the DB, ignoring all the "
        "other text.\n\n"
        "See also the 'load-tax' command.",
    )
    subcommand_parser.add_argument("-i", "--input", **ARG_INPUT_FASTA)
    subcommand_parser.add_argument("-d", "--database", **ARG_DB_WRITE)
    subcommand_parser.add_argument("-k", "--marker", **ARG_MARKER)
    subcommand_parser.add_argument(
        "-l",
        "--left",
        type=str,
        metavar="PRIMER",
        help="Left primer sequence to record in DB for this marker. "
        "Can use IUPAC ambiguity codes. Default DB uses 21bp ITS6 "
        "'GAAGGTGAAGTCGTAACAAGG' from Cooke et al. 2000 "
        "https://doi.org/10.1006/fgbi.2000.1202",
    )
    subcommand_parser.add_argument(
        "-r",
        "--right",
        type=str,
        metavar="PRIMER",
        help="Right primer sequence to record in the DB for this marker. "
        "Will find and remove the reverse complement from end of marker "
        "sequences. Can use IUPAC ambiguity codes. Default DB uses 20bp "
        "5.8S-1R primer 'GCARRGACTTTCGTCCCYRC' from Scibetta et al. 2012 "
        "https://doi.org/10.1016/j.mimet.2011.12.012 - meaning "
        "look for 'GYRGGGACGAAAGTCYYTGC' after the marker.",
    )
    subcommand_parser.add_argument(
        "--minlen",
        type=int,
        default=None,
        metavar="LENGTH",
        help="Minimum acceptable amplicon length for this import. "
        f"Default {DEF_MIN_LENGTH} if marker not already defined. "
        "Will be recorded for a new marker, if marker already defined "
        "can use a stricter value for this import.",
    )
    subcommand_parser.add_argument(
        "--maxlen",
        type=int,
        default=None,
        metavar="LENGTH",
        help="Maximum acceptable amplicon length for this import. "
        f"Default {DEF_MAX_LENGTH} if marker not already defined. "
        "Will be recorded for a new marker, if marker already defined "
        "can use a stricter value for this import.",
    )
    subcommand_parser.add_argument(
        "-n",
        "--name",
        type=str,
        default="",
        help="Name to record for this data source (string).",
    )
    subcommand_parser.add_argument(
        "-x",
        "--lax",
        default=False,
        action="store_true",
        help="Accept species names without pre-loaded taxonomy.",
    )
    subcommand_parser.add_argument(
        "-g",
        "--genus",
        default=False,
        action="store_true",
        help="Record at genus level only (and only validate at genus level, "
        "unless using -x / --lax in which case any word is accepted as a genus).",
    )
    subcommand_parser.add_argument(
        "-c",
        "--convention",
        type=str,
        default="simple",
        choices=list(fasta_conventions),
        help="Which naming convention does the FASTA file follow.",
    )
    subcommand_parser.add_argument(
        "-s",
        "--sep",
        type=str,
        default=None,
        metavar="CHAR",
        help="FASTA description multi-entry separator character. Default none "
        "meaning assume single entries.",
    )
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument("-t", "--temp", **ARG_TEMPDIR)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=db_import)
    del subcommand_parser  # To prevent accidentally adding more

    # dump
    subcommand_parser = subparsers.add_parser(
        "dump",
        description="Export a marker database to a text file.",
        epilog="e.g. 'thapbi_pict dump -d ... -g Peronospora -o Peronospora.txt'",
    )
    subcommand_parser.add_argument("-d", "--database", **ARG_DB_INPUT)
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="FILENAME",
        help="File to write to (default '-' meaning stdout)",
    )
    subcommand_parser.add_argument(
        "-k",
        "--marker",
        type=str,
        default=None,
        help="Which amplicon marker from the DB wanted? Default all. "
        "Recommend for minimal and FASTA output.",
    )
    subcommand_parser.add_argument(
        "-m",
        "--minimal",
        action="store_true",
        help="Minimal output, one record per trimmed sequence with MD5 and species "
        "list - rather than all entries with original ID, species, and taxid.",
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
    subcommand_parser.add_argument(
        "--sep",
        type=str,
        default=";",
        metavar="CHAR",
        help="FASTA description entry separator, default semi-colon.",
    )
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=dump)
    del subcommand_parser  # To prevent accidentally adding more

    # conflicts
    subcommand_parser = subparsers.add_parser(
        "conflicts",
        description="Count genus or species conflicts in a marker database.",
        epilog="Number of marker or genus level conflicts is used as the return "
        "code. e.g. 'thapbi_pict conflicts -d ... -o conflicts.txt ; echo $?'",
    )
    subcommand_parser.add_argument("-d", "--database", **ARG_DB_INPUT)
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="FILENAME",
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
        "XXX_R1_001.fastq[.gz] and XXX_R2_001.fastq[.gz], and will given an "
        "output file XXX.fasta. These are non-redundant, entries named by "
        "checksum and their abundance, and sorted by decreasing abundance "
        "then alphabetically by sequence.",
    )
    subcommand_parser.add_argument("-i", "--input", **ARG_INPUT_FASTQ)
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        metavar="DIRNAME",
        help="Output directory. Required.",
    )
    subcommand_parser.add_argument("-a", "--abundance", **ARG_FASTQ_MIN_ABUNDANCE_TWO)
    subcommand_parser.add_argument(
        "-f", "--abundance-fraction", **ARG_FASTQ_NOISE_PERC_ZERO
    )
    subcommand_parser.add_argument("-d", "--database", **ARG_DB_INPUT)
    subcommand_parser.add_argument("--flip", **ARG_FLIP)
    subcommand_parser.add_argument("--merged-cache", **ARG_MERGED_CACHE)
    subcommand_parser.add_argument("-t", "--temp", **ARG_TEMPDIR)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.add_argument("--cpu", **ARG_CPU)
    subcommand_parser.set_defaults(func=prepare_reads)
    del subcommand_parser  # To prevent accidentally adding more

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
        metavar="FILENAME",
        help="Single output filename, '-' for stdout (default). "
        "Can be a directory if a single -i/-r input file given.",
    )
    subcommand_parser.add_argument(
        "-a",
        "--abundance",
        type=int,
        default=0,
        help="Minimum abundance to require before outputting a sequence. "
        "Default no minimum.",
    )
    subcommand_parser.add_argument("--minlen", **ARG_MIN_LENGTH)
    subcommand_parser.add_argument("--maxlen", **ARG_MAX_LENGTH)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=fasta_nr)

    subcommand_parser = subparsers.add_parser(
        "denoise",
        description="Apply UNOISE based read-correction to denoise FASTA file(s).",
        epilog="This is a simplified version of the sample-tally command. "
        "Input FASTA files should use <prefix>_<count> naming. In the output "
        "FASTA file each unique sequence will be named <MD5>_<count> using "
        "the upper case sequence MD5 checksum, and its total abundance "
        "including contributions from any reads corrected to that sequences.",
    )
    subcommand_parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        nargs="+",
        metavar="FASTA",
        help="One or more input FASTA files.",
    )
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="FASTA",
        help="Single output FASTA filename, '-' for stdout (default).",
    )
    subcommand_parser.add_argument(
        "--denoise",  # Named to match sample-tally and pipeline setting
        choices=["unoise-l", "usearch", "vsearch"],
        default="unoise-l",
        help="Choice of read-correction algorithm, default 'unoise-l' for "
        "built-in reimplementation of the Edgar (2016) UNOISE algorithm "
        "using Levenshtein distance. Use 'usearch' or 'vsearch' to call "
        "external tools 'usearch -unoise3 ...' or 'vsearch --cluster_unoise "
        "...' respectively.",
    )
    subcommand_parser.add_argument("-α", "--unoise_alpha", **ARG_UNOISE_ALPHA)
    subcommand_parser.add_argument("-γ", "--unoise_gamma", **ARG_UNOISE_GAMMA)
    subcommand_parser.add_argument(
        "-t",
        "--total",
        type=int,
        default="0",
        help="Minimum total abundance for each output sequence. Default 0 (not used).",
    )
    subcommand_parser.add_argument("--minlen", **ARG_MIN_LENGTH)
    subcommand_parser.add_argument("--maxlen", **ARG_MAX_LENGTH)
    subcommand_parser.add_argument("--temp", **ARG_TEMPDIR)
    subcommand_parser.add_argument("--cpu", **ARG_CPU)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=denoise)

    subcommand_parser = subparsers.add_parser(
        "sample-tally",
        description="Prepare ASV x sample TSV sequence file using MD5 naming.",
        epilog="Each unique sequence will be named <Marker>/<MD5>_<count> "
        "using the upper case sequence MD5 checksum, and its total abundance. "
        "Per-sample input FASTA files should use <prefix>_<count> naming. "
        "Output is a plain text tab-separated table with one line per unique "
        "sequence (ASV), named in the first column, then one column per "
        "sample, and the uppercase sequence as the final column. The ASV vs "
        "sample values are the counts from the input FASTA files (pooling any "
        "read-corrected sequences if using denoising).",
    )
    subcommand_parser.add_argument(
        "-k",
        "--marker",
        type=str,
        default=None,
        help="Which amplicon marker to process. Required if the DB has more "
        "than one defined.",
    )
    subcommand_parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        nargs="+",
        metavar="INPUT",
        help="One or more per-sample FASTA files.",
    )
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument("-y", "--synctrls", **ARG_SYN_CONTROLS_FASTA)
    subcommand_parser.add_argument("-n", "--negctrls", **ARG_NEG_CONTROLS_FASTA)
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="FILENAME",
        help="Single output filename, '-' for stdout (default). The pipeline "
        "will use the .tally.tsv extension which is assumed for the summary "
        "and edit-graph commands.",
    )
    subcommand_parser.add_argument(
        "--fasta",
        type=str,
        metavar="FILENAME",
        help="Optional output FASTA filename, '-' for stdout. ",
    )
    subcommand_parser.add_argument("-d", "--database", **ARG_DB_INPUT)
    subcommand_parser.add_argument("--synthetic", **ARG_SYNTHETIC_SPIKE)
    subcommand_parser.add_argument("-a", "--abundance", **ARG_FASTQ_MIN_ABUNDANCE)
    subcommand_parser.add_argument("-f", "--abundance-fraction", **ARG_FASTQ_NOISE_PERC)
    subcommand_parser.add_argument(
        "-t",
        "--total",
        type=int,
        default="0",
        help="Minimum total abundance for sequences. "
        "Applied after per-sample thresholds. Default 0 (not used).",
    )
    subcommand_parser.add_argument("--minlen", **ARG_MIN_LENGTH)
    subcommand_parser.add_argument("--maxlen", **ARG_MAX_LENGTH)
    subcommand_parser.add_argument("--denoise", **ARG_DENOISE)
    subcommand_parser.add_argument("-α", "--unoise_alpha", **ARG_UNOISE_ALPHA)
    subcommand_parser.add_argument("-γ", "--unoise_gamma", **ARG_UNOISE_GAMMA)
    subcommand_parser.add_argument("--temp", **ARG_TEMPDIR)
    subcommand_parser.add_argument("--cpu", **ARG_CPU)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=sample_tally)

    # classify
    subcommand_parser = subparsers.add_parser(
        "classify",
        description="Classify FASTA file of marker sequences by species.",
        epilog="Each input file XXX.fasta will result in an output file "
        "named XXX.method.tsv in the specified output directory (default "
        "input dir).",
    )
    subcommand_parser.add_argument("-i", "--input", **ARG_INPUT_FASTA_OR_TSV)
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument(
        "-a",
        "--abundance",
        type=int,
        default=0,
        help="Minimum abundance applied to unique marker sequences in each "
        "FASTA sample file, default 0 (classify all).",
    )
    subcommand_parser.add_argument("-d", "--database", **ARG_DB_INPUT)
    subcommand_parser.add_argument("-k", "--marker", **ARG_MARKER_PICK)
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
    subcommand_parser.add_argument("-b", "--biom", **ARG_BIOM)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.add_argument("--cpu", **ARG_CPU)
    subcommand_parser.set_defaults(func=classify)
    del subcommand_parser  # To prevent accidentally adding more

    # assess-classification
    subcommand_parser = subparsers.add_parser(
        "assess",
        description="Assess accuracy of marker sequence classification.",
        epilog="Takes input sample-tally classification files (*.<method>.tsv), "
        "to be assessed against sample-tally classification files (*.<known>.tsv) "
        "or legacy per-sample expected classification files (<sample>.known.tsv). "
        "Produces a multi-species confusion matrix (output on request) and "
        "classifier performance metrics (to stdout by default). Assessment is "
        "at sample level (taking the union of species predicted by all "
        "sequences from each sample).",
    )
    subcommand_parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        nargs="+",
        help="TSV filenames or folders named *.<method>.tsv, and *.<known>.tsv "
        "where the filename suffixes can be set via the -m / --method and "
        "--known arguments.",
    )
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument(
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
    subcommand_parser.add_argument(
        # "-k",
        "--marker",
        type=str,
        default=None,
        help="Which amplicon marker from the DB wanted? Default all. "
        "Use if assessing a single marker to avoid extra false negatives.",
    )
    subcommand_parser.add_argument("-m", "--method", **ARG_METHOD_INPUT)
    subcommand_parser.add_argument(
        "-a",
        "--abundance",
        type=int,
        default="1",
        help="Minimum abundance to require before considering a classification. "
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
    del subcommand_parser  # To prevent accidentally adding more

    # summary
    subcommand_parser = subparsers.add_parser(
        "summary",
        description="Sample and sequence summary reports on classifier output.",
        epilog="Assumes you've run prepare-reads, sample-tally and classify "
        "giving tally-file XXX.tsv (with sequences counts per sample), and "
        "classifier output XXX.method.tsv (for the same sequences). "
        "The output is two sets of tables. The read tables have one row per "
        "unique sequence (can be thousands of rows) and one column per sample "
        "(often hundreds, typically 96 samples per plate). The sample tables "
        "have one row per sample, and one column per genus and species.",
    )
    subcommand_parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        nargs="+",
        help="Sample tally file (*.tally.tsv) and matching classifier output "
        "file (*.method.tsv), or folder name(s). "
        "The classifier method extension can be set via -m / --method.",
    )
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        metavar="STEM",
        help="Output filename stem for generating reports. Not a directory.",
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
        help="Minimum sample level abundance to require for the report. "
        f"Default {DEFAULT_MIN_ABUNDANCE} reflects default in prepare-reads. "
        "Rather than re-running the prepare or classifier steps with a stricter "
        "minimum abundance you can apply it here. Use zero or one to look at "
        "everything in your input files.",
    )
    subcommand_parser.add_argument("-t", "--metadata", **ARG_METADATA)
    subcommand_parser.add_argument("-e", "--metaencoding", **ARG_METAENCODING)
    subcommand_parser.add_argument("-c", "--metacols", **ARG_METACOLS)
    subcommand_parser.add_argument("-x", "--metaindex", **ARG_METAINDEX)
    subcommand_parser.add_argument("-g", "--metagroups", **ARG_METAGROUPS)
    subcommand_parser.add_argument("--metafields", **ARG_METAFIELDS)
    subcommand_parser.add_argument("-q", "--requiremeta", **ARG_REQUIREMETA)
    subcommand_parser.add_argument("-u", "--unsequenced", **ARG_UNSEQUENCED)
    subcommand_parser.add_argument("-b", "--biom", **ARG_BIOM)
    subcommand_parser.add_argument("-v", "--verbose", **ARG_VERBOSE)
    subcommand_parser.set_defaults(func=summary)
    del subcommand_parser  # To prevent accidentally adding more

    # edit-graph
    subcommand_parser = subparsers.add_parser(
        "edit-graph",
        description="Draw network graph of marker sequences using edit distance.",
        epilog="Takes a marker database and/or prepared sample-tally file as input. "
        "The output is a network graph (in a choice of format) with unique "
        "sequences as nodes (in the PDF labelled by the database taxonomy, "
        "colored by genus, size set by total abundance in the FASTA files), "
        "and short edit distances as edges between nodes. For Cytoscape "
        "we recommend generating XGMML output here, then start Cytoscape, "
        "menu 'File', 'Import', 'Import from file', and then run a layout. "
        "Both 'Perfuse Force Directed' and 'Edge-weighted Spring Embedded' "
        "work well.",
    )
    arg = subcommand_parser.add_argument("-d", "--database", **ARG_DB_INPUT)
    arg.help += " Used for labels and colors. Use '' to mean no DB."
    del arg
    subcommand_parser.add_argument(
        "-i",
        "--input",
        type=str,
        metavar="FILENAME",
        required=False,
        help="Sample-tally TSV file (*.tally.tsv) or classifier output "
        "(*.<method>.tsv)",
    )
    subcommand_parser.add_argument("--ignore-prefixes", **ARG_IGNORE_PREFIXES)
    subcommand_parser.add_argument(
        "-a",
        "--abundance",
        type=int,
        default=str(DEFAULT_MIN_ABUNDANCE),
        help="Minimum sample level abundance for FASTA sequences. "
        f"Default {DEFAULT_MIN_ABUNDANCE} reflects default in prepare-reads.",
    )
    subcommand_parser.add_argument(
        "-t",
        "--total",
        type=int,
        default="0",
        help="Minimum total abundance for FASTA sequences. "
        "Applied after per-sample level minimum (-a / --abundance). "
        "Offered as a way to simplify the final graph.",
    )
    subcommand_parser.add_argument(
        "--min-samples",
        type=int,
        default="0",
        help="Minimum number of samples a unique sequence must appear in. "
        "Offered as a way to simplify the final graph. Default 0.",
    )
    subcommand_parser.add_argument(
        "-k",
        "--marker",
        type=str,
        default=None,
        help="Show only entries from this marker. Required if only drawing DB entries.",
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
        choices=["graphml", "gexf", "gml", "xgmml", "pdf", "matrix"],
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
        metavar="FILENAME",
        help="File to write to (default '-' meaning stdout)",
    )
    subcommand_parser.add_argument("-t", "--metadata", **ARG_METADATA)
    subcommand_parser.add_argument("-e", "--metaencoding", **ARG_METAENCODING)
    subcommand_parser.add_argument("-c", "--metacols", **ARG_METACOLS)
    subcommand_parser.add_argument("-x", "--metaindex", **ARG_METAINDEX)
    subcommand_parser.add_argument("--metafields", **ARG_METAFIELDS)
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
