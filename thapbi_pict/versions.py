# Copyright 2019-2024 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Helper code to get command line tool versions.

Defines various functions to check a tool is on the ``$PATH`` and if so,
return the tool version as a short string (sometimes including a date).

These functions are called from various THAPBI-PICT subcommands which call
external tools to ensure a clear missing dependency message, and to log the
version of the external tool used.

If the tool is not on the path, the commands all return None.

If we cannot parse the output, again the commands return None - which is
likely an indication of a major version change, meaning the tool ought to be
re-evaluated for use with THAPBI-PICT.
"""

from __future__ import annotations

import sys
from subprocess import getoutput


def check_rapidfuzz() -> str:
    """Check can import rapidfuzz and confirm recent enough."""
    try:
        import rapidfuzz
    except ImportError:
        sys.exit("ERROR: Missing Python library rapidfuzz")
    try:
        version: str = rapidfuzz.__version__
    except AttributeError:
        sys.exit("ERROR: Could not check rapidfuzz.__version__")
    try:
        major_minor = tuple(int(_) for _ in version.split(".")[:2])
    except ValueError:
        sys.exit(
            f"ERROR: Could not parse rapidfuzz.__version__ {rapidfuzz.__version__}"
        )
    if major_minor < (2, 4):
        sys.exit(f"ERROR: Need at least rapidfuzz v2.4, have {rapidfuzz.__version__}")
    return version


def check_tools(names: list[str], debug: bool) -> list[str]:
    """Verify the named tools are present, log versions if debug=True.

    Argument names should be an iterable of tool binary names.

    If all the tools are present, returns a list of version strings.

    If any tools are missing (or have a version we could not parse),
    aborts.
    """
    easy = {
        "blastn": version_blast,
        "cutadapt": version_cutadapt,
        "flash": version_flash,
        "fdp": version_graphviz_fdp,
        "makeblastdb": version_blast,
        "usearch": version_usearch,
        "vsearch": version_vsearch,
    }
    missing = []
    versions = []
    for name in names:
        if name in easy:
            # Just call the associated function, with the binary name
            version = easy[name](name)
            if not version:
                missing.append(name)
            else:
                versions.append(name)
                if debug:
                    sys.stderr.write(f"DEBUG: version of {name}: {version}\n")
        else:
            sys.exit(f"ERROR: Unsupported external tool name, {name}")
    if len(missing) == 1:
        sys.exit("ERROR: Missing external tool: " + missing[0])
    elif missing:
        sys.exit("ERROR: Missing external tool(s): " + ",".join(missing))
    else:
        return versions


def version_graphviz_fdp(cmd: str = "fdp") -> str | None:
    """Return the version of the GraphViz tool fdp (as a short string).

    Depends on the -V switch::

        $ fdp -V
        fdp - graphviz version 9.0.0 (0)

    In the above example, it would behave as follows:

    >>> version_graphviz_fdp()
    '9.0.0'

    If the command is not on the path, returns None.
    """
    text = getoutput(cmd + " -V")
    for line in text.splitlines():
        if line.startswith("fdp - graphviz version "):
            words = line.strip().split()
            return words[4]
    return None


def version_blast(cmd: str = "blastn") -> str | None:
    """Return the version of the NCBI BLAST+ suite's blastn (as a short string).

    In the absence of a built in version switch like ``-v``, this works by
    parsing the short help output with ``-h`` (which does vary between the
    tools in the suite)::

        $ makeblastdb -h | grep BLAST
        Application to create BLAST databases, version 2.16.0+

        $ blastn -h | grep BLAST
        Nucleotide-Nucleotide BLAST 2.16.0+

    In the above examples, it would behave as follows:

    >>> version_blast("makeblastdb")
    '2.16.0+'
    >>> version_blast("blastn")
    '2.16.0+'

    If the command is not on the path, returns None.
    """
    text = getoutput(cmd + " -h")
    for line in text.splitlines():
        if line.strip().endswith("+"):
            words = line.strip().split()
            if "BLAST" in words or "version" in words:
                return words[-1]
    return None


def version_cutadapt(cmd: str = "cutadapt") -> str | None:
    """Return the version of cutadapt (as a short string).

    Uses the output with ``--version``::

        $ cutadapt --version
        4.7

    It would capture this:

    >>> version_cutadapt()
    '4.7'

    If the command is not on the path, returns None.
    """
    text = getoutput(cmd + " --version").strip().split("\n", 1)[0]
    if "." in text:
        return text
    return None


def version_flash(cmd: str = "flash") -> str | None:
    """Return the version of flash (as a short string).

    Parses the output with ``-v``::

        $ flash -v | head -n 1
        FLASH v1.2.11

    It would capture the version from the first line as follows:

    >>> version_flash()
    '1.2.11'

    If the command is not on the path, returns None.
    """
    text = getoutput(cmd + " -v")
    ver = text.split("\n", 1)[0]
    if ver.upper().startswith("FLASH V"):
        return ver[7:]
    return None


def version_usearch(cmd: str = "usearch") -> str | None:
    """Return the version of usearch (as a short string).

    Uses the output with ``--version``::

        $ usearch --version
        usearch v11.0.667_i86linux32

    It would capture this:

    >>> version_usearch()  # doctest: +SKIP
    '11.0.667'

    If the command is not on the path, returns None.
    """
    text = getoutput(cmd + " --version").strip().split("\n", 1)[0]
    ver = text.split("_", 1)[0]
    if ver.lower().startswith("usearch v"):
        return ver[9:]
    return None


def version_vsearch(cmd: str = "vsearch") -> str | None:
    """Return the version of vsearch (as a short string).

    Uses the output with ``--version``::

        $ vsearch --version
        vsearch v2.28.1_macos_x86_64, 8.0GB RAM, 8 cores
        ...

    It would capture this:

    >>> version_vsearch()
    '2.28.1'

    If the command is not on the path, returns None.
    """
    # Seems to be first line on Linux/Mac, but not on Windows
    # Could be due to stdout vs stderr sorting/buffering?
    for text in getoutput(cmd + " --version").strip().split("\n"):
        ver = text.split(",", 1)[0].split("_", 1)[0]
        if ver.lower().startswith("vsearch v"):
            return ver[9:]
    return None
