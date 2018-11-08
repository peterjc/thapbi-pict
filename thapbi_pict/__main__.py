"""This defines the thapbi_pict command line tool.

This works via ``setup.py`` where under ``entry_points`` we define a
``console_scripts`` entry for ``thapbi_pict`` (executable name) pointing to
the ``main()`` function define in this Python file.
"""

import argparse
import sys

from . import __version__


def expand_database_argument(text):
    """Expand an SQLite3 filename to an SQLalchemy URL."""
    # TODO: Expand this to allow other DB prefixes later
    # Note we are not currently checking file exists,
    # as we might be about to create it.
    prefix = "sqlite:///"
    if text.startswith(prefix):
        return text
    return prefix + text


def legacy_import(args=None):
    """Subcommand to import a legacy ITS1 FASTA file into a database."""
    from .legacy import main
    return main(
        fasta_files=args.fasta,
        db_url=expand_database_argument(args.database),
        debug=args.verbose)


def dump(args=None):
    """Subcommand to dump a database to a text file."""
    from .dump import main
    return main(db_url=expand_database_argument(args.database),
                output_txt=args.output,
                clade=args.clade,
                debug=args.verbose)


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
        description=("THAPBI Phytophthora ITS1 Classifier Tool (PICT), v%s."
                     % __version__),
        epilog="e.g. run 'thapbi_pict dump -h' for the dump subcommand help.")
    parser.add_argument("-v", "--version", action="version",
                        version="THAPBI PICT v%s" % __version__)
    subparsers = parser.add_subparsers(title="subcommands",
                                       description="valid subcommands",
                                       help="additional help")

    # legacy-import
    parser_import = subparsers.add_parser(
        "legacy-import",
        description="Load one of our legacy ITS1 FASTA files into a database.")
    parser_import.add_argument(
        'fasta', nargs='+',
        help='One or more ITS1 fasta filenames.')
    parser_import.add_argument(
        "-d", "--database", type=str,
        help="Which database to write to (or create)")
    parser_import.add_argument(
        "-v", "--verbose", action='store_true',
        help="Verbose logging")
    parser_import.add_argument(
        "-t", "--taxid", type=int, default=0,
        help="Default NCBI taxid (integer)")
    parser_import.set_defaults(func=legacy_import)

    # dump
    parser_dump = subparsers.add_parser(
        "dump",
        description="Export an ITS1 database to a text file.",
        epilog="e.g. 'thapbi_pict dump -d ... -c 8a,8b -o clade_8a_8b.txt'")
    parser_dump.add_argument(
        "-d", "--database", type=str, required=True,
        help="Which database to export from")
    parser_dump.add_argument(
        "-o", "--output", type=str, default="-",
        help="File to write to (default '-' meaning stdout)")
    parser_dump.add_argument(
        "-c", "--clade", type=str, default="",
        help="Which clade(s) to export (comma separated list, "
             "with '-' meaning no clade defined). "
             "Default is not to filter by clade.")
    parser_dump.add_argument(
        "-v", "--verbose", action='store_true',
        help="Verbose logging")
    parser_dump.set_defaults(func=dump)

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
