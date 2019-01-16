"""This defines the thapbi_pict command line tool.

This works via ``setup.py`` where under ``entry_points`` we define a
``console_scripts`` entry for ``thapbi_pict`` (executable name) pointing to
the ``main()`` function define in this Python file.
"""

import argparse
import os
import sys

from . import __version__
from .classify import method_classifier


def check_output_directory(out_dir):
    """Command line validation of output directory value."""
    if out_dir == "-" or os.path.isdir(out_dir):
        return True
    elif os.path.isfile(out_dir):
        sys.exit("ERROR: Output directory name is a file: %s\n" % out_dir)
    else:
        sys.exit("ERROR: Output directory does not exist")


def expand_database_argument(text):
    """Expand an SQLite3 filename to an SQLalchemy URL."""
    # TODO: Expand this to allow other DB prefixes later
    # Note we are not currently checking file exists,
    # as we might be about to create it.
    if not text:
        sys.exit("The database argument is required.\n")
    prefix = "sqlite:///"
    if text.startswith(prefix):
        return text
    return prefix + text


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
        fasta_file=args.fasta,
        db_url=expand_database_argument(args.database),
        name=args.name,
        validate_species=args.validate_species,
        debug=args.verbose,
    )


def legacy_import(args=None):
    """Subcommand to import a legacy ITS1 FASTA file into a database."""
    from .legacy import main

    return main(
        fasta_file=args.fasta,
        db_url=expand_database_argument(args.database),
        name=args.name,
        validate_species=args.validate_species,
        debug=args.verbose,
    )


def dump(args=None):
    """Subcommand to dump a database to a text file."""
    from .dump import main

    return main(
        db_url=expand_database_argument(args.database),
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
    return main(
        fastq=args.fastq,
        controls=args.controls,
        out_dir=args.output,
        min_abundance=args.abundance,
        debug=args.verbose,
        cpu=args.cpu,
    )


def classify(args=None):
    """Subcommand to classify ITS1 sequences using a database."""
    from .classify import main

    check_output_directory(args.output)
    return main(
        fasta=args.fasta,
        db_url=expand_database_argument(args.database),
        method=args.method,
        out_dir=args.output,
        debug=args.verbose,
        cpu=args.cpu,
    )


def assess_classification(args=None):
    """Subcommand to assess classification using known true taxonomy."""
    from .assess import main

    check_output_directory(args.output)
    return main(
        fasta=args.fasta,
        known=args.known,
        method=args.method,
        out_dir=args.output,
        debug=args.verbose,
    )


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
        epilog="e.g. run 'thapbi_pict dump -h' for the dump subcommand help.",
    )
    parser.add_argument(
        "-v", "--version", action="version", version="THAPBI PICT v%s" % __version__
    )
    subparsers = parser.add_subparsers(
        title="subcommands", help="Each subcommand has its own additional help"
    )

    # load-tax
    parser_load_tax = subparsers.add_parser(
        "load-tax", description="Load an NCBI taxonomy dump into an ITS1 database."
    )
    parser_load_tax.add_argument(
        "-t",
        "--tax",
        type=str,
        required=True,
        help="Folder containing NCBI taxonomy dump files names.dmp etc.",
    )
    parser_load_tax.add_argument(
        "-d",
        "--database",
        type=str,
        required=True,
        help="Which database to write to (or create)",
    )
    parser_load_tax.add_argument(
        "-a",
        "--ancestors",
        type=str,
        default="4783",
        help="Comma separated lists of taxids at genus level or higher, "
        "default 4783 for Phytophthora",
    ),
    parser_load_tax.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose logging"
    )
    parser_load_tax.set_defaults(func=load_tax)

    # ncbi-import
    parser_ncbi_import = subparsers.add_parser(
        "ncbi-import",
        description="Load an NCBI format ITS1 FASTA file into a database.",
    )
    parser_ncbi_import.add_argument("fasta", type=str, help="One ITS1 fasta filename.")
    parser_ncbi_import.add_argument(
        "-d",
        "--database",
        type=str,
        required=True,
        help="Which database to write to (or create)",
    )
    parser_ncbi_import.add_argument(
        "-n",
        "--name",
        type=str,
        default="",
        help="Data source name (string, ideally avoiding spaces etc)",
    )
    parser_ncbi_import.add_argument(
        "-s",
        "--validate_species",
        default=False,
        action="store_true",
        help="Only load ITS1 entries matching a known species name",
    )
    parser_ncbi_import.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose logging"
    )
    parser_ncbi_import.set_defaults(func=ncbi_import)

    # legacy-import
    parser_legacy_import = subparsers.add_parser(
        "legacy-import",
        description="Load one of our legacy ITS1 FASTA files into a database.",
    )
    parser_legacy_import.add_argument(
        "fasta", type=str, help="One ITS1 fasta filename."
    )
    parser_legacy_import.add_argument(
        "-d",
        "--database",
        type=str,
        required=True,
        help="Which database to write to (or create)",
    )
    parser_legacy_import.add_argument(
        "-n",
        "--name",
        type=str,
        default="",
        help="Data source name (string, ideally avoiding spaces etc)",
    )
    parser_legacy_import.add_argument(
        "-s",
        "--validate_species",
        default=False,
        action="store_true",
        help="Only load ITS1 entries matching a known species name",
    )
    parser_legacy_import.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose logging"
    )
    parser_legacy_import.set_defaults(func=legacy_import)

    # dump
    parser_dump = subparsers.add_parser(
        "dump",
        description="Export an ITS1 database to a text file.",
        epilog="""Examples:

$ thapbi_pict dump -d ... -c 8a,8b -o clade_8a_8b.txt

$ thapbi_pict dump -d ... -g Phytophthora -s "ilicis, sp. aff. meadii" -o Phytophthora.txt

Note that the -s or --species option allows spaces after the
comma.
""",  # noqa: E501
    )
    parser_dump.add_argument(
        "-d",
        "--database",
        type=str,
        required=True,
        help="Which database to export from",
    )
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
        help="Which genus (or genera) export (comma separated list). "
        "Default is not to filter by genus.",
    )
    parser_dump.add_argument(
        "-s",
        "--species",
        type=str,
        default="",
        help="Which species to export (comma separated list). "
        "Requires a single genus argument be given. "
        "Default is not to filter by species.",
    )
    parser_dump.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose logging"
    )
    parser_dump.set_defaults(func=dump)

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
    parser_prepare_reads.add_argument(
        "fastq",
        type=str,
        nargs="+",
        help="One or more ITS1 FASTQ filenames or folder names "
        "(containing files named *.fastq or *.fastq.gz).",
    )
    parser_prepare_reads.add_argument(
        "-c",
        "--controls",
        type=str,
        nargs="+",
        help="One or more control FASTQ filenames or folder names "
        "(which can also be included in the FASTQ argument). "
        "The paired FASTQ controls reads are processed in order to"
        "determine the minimal abundance threshold automatically.",
    )
    parser_prepare_reads.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="DIRNAME",
        help="Directory to write output FASTA files to, "
        "default is next to each input file.",
    )
    parser_prepare_reads.add_argument(
        "-a",
        "--abundance",
        type=int,
        default="100",
        help="Mininum abundance to apply to final candidate ITS1 "
        "sequences in the output FASTA file (default 100). "
        "This may be increased based on any FASTQ controls.",
    )
    parser_prepare_reads.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose logging"
    )
    parser_prepare_reads.add_argument(
        "--cpu",
        type=int,
        default=0,
        help="Number of parallel threads to use in called tools.",
    )
    parser_prepare_reads.set_defaults(func=prepare_reads)

    # classify
    parser_classify = subparsers.add_parser(
        "classify",
        description="Classify FASTA file of ITS1 sequences by species.",
        epilog="Each input file XXX.fasta will result in output files "
        "namesd XXX.method-reads.tsv and XXX.method-tax.tsv in "
        "the specified output directory (default input dir).",
    )
    parser_classify.add_argument(
        "fasta",
        type=str,
        nargs="+",
        help="One or more ITS1 FASTA filenames or folder names "
        "(containing files named *.fasta).",
    )
    parser_classify.add_argument(
        "-d",
        "--database",
        type=str,
        required=True,
        help="Which ITS1 database to use for species classification.",
    )
    parser_classify.add_argument(
        "-m",
        "--method",
        type=str,
        default="identity",
        choices=list(method_classifier),
        help="Method to use, default uses simple identity.",
    )
    parser_classify.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="DIRNAME",
        help="Directory to write output reports to, default "
        "is next to each input file.",
    )
    parser_classify.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose logging"
    )
    parser_classify.add_argument(
        "--cpu",
        type=int,
        default=0,
        help="Max number of parallel threads to use in called tools.",
    )
    parser_classify.set_defaults(func=classify)

    # assess-classification
    parser_assess = subparsers.add_parser(
        "assess",
        description="Assess accuracy of ITS1 read classification.",
        epilog="Takes as input XXX.known-reads.tsv and matching "
        "predictions in XXX.method-reads.tsv (in same directory) "
        "to produce a multi-species confusion matrix named "
        "XXX.method-vs-known.tsv, and a summary to stdout. You can "
        "deliberately compare to prediction methods to each other "
        "using this.",
    )
    parser_assess.add_argument(
        "fasta",
        type=str,
        nargs="+",
        help="One or more FASTA file or folder names. Next to each "
        "FASTA file expects matching files *.method-reads.tsv to be "
        "assessed against *.known-reads.tsv, where these filenames "
        "can be set via -m / --method and -k / --known arguments. ",
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
    parser_assess.add_argument(
        "-m",
        "--method",
        type=str,
        default="identity",
        choices=list(method_classifier),
        help="Method to assess (used to infer filenames), default is identity.",
    )
    parser_assess.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        metavar="DIRNAME",
        help="Directory to write output reports to, default "
        "is next to each input file.",
    )
    parser_assess.add_argument(
        "-v", "--verbose", action="store_true", help="Verbose logging"
    )
    parser_assess.set_defaults(func=assess_classification)

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
