"""Helper functions for THAPB-PICT code."""

import subprocess
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser


def cmd_as_string(cmd):
    """Express a list command as a suitably quoted string.

    Intended for using in debugging or error messages.
    """
    if isinstance(cmd, list):
        # Quote any entries with spaces
        return ' '.join('"%s"' % _ if ' ' in _ else _ for _ in cmd)
    else:
        return cmd


def run(cmd, debug=False):
    """Run a command via subprocess, abort if fails."""
    if debug:
        sys.stderr.write("Calling command: %s\n" % cmd_as_string(cmd))
    try:
        # On Python 3.7 onwards, could use capture_output=True
        # rather than stdout=PIPE and stderr=PIPE
        if isinstance(cmd, list):
            return subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                check=True)
        else:
            return subprocess.run(
                cmd, shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                check=True)
    except subprocess.CalledProcessError as e:
        if debug:
            # Used universal_newlines=True above so that this just works
            # (equivalent to text=True in Python 3.7 onwards):
            sys.stdout.write(e.stdout)
            sys.stderr.write(e.stderr)
        sys.exit("This command failed with return code %i:\n%s\n"
                 % (e.returncode, cmd_as_string(cmd)))


def abundance_from_read_name(text, debug=False):
    """Extract abundance from SWARM style read name.

    >>> abundance_from_read_name(">9e8f051c64c2b9cc3b6fcb27559418ca_988")
    988

    If fails, will return one.
    """
    try:
        return int(text.rsplit("_", 1)[1])
    except (ValueError, IndexError):
        if debug:
            sys.stderr.write(
                "WARNING: No abundance suffix in %r\n" % text)
        return 1


def abundance_values_in_fasta(fasta_file):
    """Return total and maximum abundance encoded in read names."""
    total_a = max_a = 0
    with open(fasta_file) as handle:
        for title, seq in SimpleFastaParser(handle):
            a = abundance_from_read_name(title.split(None, 1)[0])
            max_a = max(a, max_a)
            total_a += a
    return total_a, max_a
