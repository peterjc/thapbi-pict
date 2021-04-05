# Copyright 2019-2020 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Helper code to get command line tool versions.

Defines various functions to check a tool is on the $PATH and if so,
return the tool version as a short string (sometimes including a date).

These functions are called from various THAPBI-PICT subcommands which
call external tools to ensure a clear missing dependency message, and
to log the version of the external tool used.

If the tool is not on the path, the commands all return None.

If we cannot parse the output, again the commands return None - which
is likely an indication of a major version change, meaning the tool
ought to be re-evaluated for use with THAPBI-PICT.
"""
import sys
from subprocess import getoutput


def check_tools(names, debug):
    """Verify the named tools are present, log versions if debug=True.

    Argument names should be an interable of tool binary names.

    If all the tools are present, returns a list of version strings.

    If any tools are missing (or have a version we would not parse),
    aborts.
    """
    easy = {
        "blastn": version_blast,
        "cutadapt": version_cutadapt,
        "flash": version_flash,
        "hmmscan": version_hmmer,
        "makeblastdb": version_blast,
        "swarm": version_swarm,
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
    if missing:
        sys.exit("ERROR: Missing external tool(s): " + ",".join(missing))
    else:
        return versions


def version_blast(cmd="blastn"):
    """Return the version of the NCBI BLAST+ suite's blastn (as a short string).

    In the absense of a version switch, works by parsing the short
    help output with ``-h`` (which does vary between the tools in
    the suite)::

        $ makeblastdb -h | grep BLAST
        Application to create BLAST databases, version 2.7.1+

        $ blastn -h | grep BLAST
        Nucleotide-Nucleotide BLAST 2.7.1+

    In the above examples, it would behave as follows:

    >>> version_blast("makeblastdb")
    '2.7.1+'
    >>> version_blast("blastn")
    '2.7.1+'

    If the command is not on the path, returns None.
    """
    text = getoutput(cmd + " -h")
    for line in text.split("\n"):
        if line.strip().endswith("+"):
            words = line.strip().split()
            if "BLAST" in words or "version" in words:
                return words[-1]


def version_cutadapt(cmd="cutadapt"):
    """Return the version of cutadapt (as a short string).

    Uses the output with ``--version``::

        $ cutadapt --version
        1.18

    It would capture this:

    >>> version_cutadapt()
    '1.18'

    If the command is not on the path, returns None.
    """
    text = getoutput(cmd + " --version").strip().split("\n", 1)[0]
    if "." in text:
        return text


def version_hmmer(cmd="hmmscan"):
    """Return the version of hmmer (as a short string).

    Parses the output with ``-h``::

        $ hmmscan -h | grep HMMER
        # HMMER 3.2.1 (June 2018); http://hmmer.org/

    In this example, it would behave as follows:

    >>> version_hmmer()
    '3.2.1 (June 2018)'

    If the command is not on the path, returns None.
    """
    text = getoutput(cmd + " -h")
    for line in text.split("\n"):
        if line.startswith("# HMMER"):
            if line.endswith("; http://hmmer.org/"):
                line = line.rsplit(";", 1)[0]
            return line[7:].strip()


def version_flash(cmd="flash"):
    """Return the version of flash (as a short string).

    Parses the output with ``-v``::

        $ flash -v | head -n 1
        FLASH v1.2.11

    It would capture the version from the first line as follows:

    >>> version_flash()
    'v1.2.11'

    If the command is not on the path, returns None.
    """
    text = getoutput(cmd + " -v")
    ver = text.split("\n", 1)[0]
    if ver.upper().startswith("FLASH V"):
        return ver[7:]


def version_swarm(cmd="swarm"):
    """Return the version of swarm (as a short string).

    e.g. Given the following output on stderr::

        $ swarm -v
        Swarm 2.2.2 [Jun 21 2018 17:48:50]
        Copyright (C) 2012-2017 Torbjorn Rognes and Frederic Mahe
        https://github.com/torognes/swarm

        Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2014)
        Swarm: robust and fast clustering method for amplicon-based studies
        PeerJ 2:e593 https://doi.org/10.7717/peerj.593

        Mahe F, Rognes T, Quince C, de Vargas C, Dunthorn M (2015)
        Swarm v2: highly-scalable and high-resolution amplicon clustering
        PeerJ 3:e1420 https://doi.org/10.7717/peerj.1420

    It would capture the version from the first line, as follows:

    >>> version_swarm()
    '2.2.2 [Jun 21 2018 17:48:50]'

    If the command is not on the path, returns None.
    """
    text = getoutput(cmd + " -v")
    ver = text.split("\n", 1)[0].strip()
    if ver.lower().startswith("swarm "):
        return ver[6:].strip()
