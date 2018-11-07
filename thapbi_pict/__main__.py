"""This defines the thapbi_pict command line tool.

This works via ``setup.py`` where under ``entry_points`` we define a
``console_scripts`` entry for ``thapbi_pict`` (executable name) pointing to
the ``main()`` function define in this Python file.
"""

import sys

from . import __version__


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

    if "-v" in args or "--version" in args:
        print("THAPBI PICT v%s" % __version__)
        sys.exit(0)

    # TODO - Add argument parsing here...
    print("This is THAPBI PICT v%s" % __version__)


if __name__ == "__main__":
    main()
