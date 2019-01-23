"""Helper functions for THAPB-PICT code."""

import os
import subprocess
import sys

from Bio.SeqIO.FastaIO import SimpleFastaParser


def cmd_as_string(cmd):
    """Express a list command as a suitably quoted string.

    Intended for using in debugging or error messages.
    """
    if isinstance(cmd, list):
        # Quote any entries with spaces
        return " ".join('"%s"' % _ if " " in _ else _ for _ in cmd)
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
                check=True,
            )
        else:
            return subprocess.run(
                cmd,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                check=True,
            )
    except subprocess.CalledProcessError as e:
        if debug:
            # Used universal_newlines=True above so that this just works
            # (equivalent to text=True in Python 3.7 onwards):
            sys.stdout.write(e.stdout)
            sys.stderr.write(e.stderr)
        sys.exit(
            "This command failed with return code %i:\n%s\n"
            % (e.returncode, cmd_as_string(cmd))
        )


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
            sys.stderr.write("WARNING: No abundance suffix in %r\n" % text)
        return 1


def abundance_values_in_fasta(fasta_file):
    """Return total and maximum abundance encoded in read names."""
    total_a = max_a = 0
    with open(fasta_file) as handle:
        for title, _ in SimpleFastaParser(handle):
            a = abundance_from_read_name(title.split(None, 1)[0])
            max_a = max(a, max_a)
            total_a += a
    return total_a, max_a


def find_requested_files(filenames_or_folders, ext=".fasta", debug=False):
    """Interpret a list of filenames and/or foldernames.

    The extensions argument can be a tuple.
    """
    answer = []
    for x in filenames_or_folders:
        if os.path.isdir(x):
            if debug:
                sys.stderr.write("Walking directory %r\n" % x)
            for root, _, files in os.walk(x):
                for f in files:
                    if f.endswith(ext):
                        # Check not a directory?
                        answer.append(os.path.join(root, f))
        elif os.path.isfile(x):
            if x.endswith(ext):
                answer.append(x)
            else:
                sys.exit(
                    "Specified filename %r does not have expected extension %r."
                    % (x, ext)
                )
        else:
            sys.exit("ERROR: %r is not a file or a directory\n" % x)
    # Warn if there were duplicates?
    return sorted(set(answer))


def find_paired_files(filenames_or_folders, ext1, ext2, debug=False):
    """Interpret a list of filenames and/or foldernames to find pairs.

    Looks for paired files named XXX.ext1 and XXX.ext2 which can be
    in different directories - duplicated filenames (in different
    directories) are considered to be an error.

    Having XXX.ext1 without XXX.ext2 is treated as a warning.

    Having XXX.ext2 without XXX.ext1 is ignored without warning.

    The arguments ext1 and ext2 should include the leading dot.
    """
    file_list = find_requested_files(
        filenames_or_folders, ext=(ext1, ext2), debug=debug
    )
    ext1_list = [_ for _ in file_list if _.endswith(ext1)]
    ext2_list = [_ for _ in file_list if _.endswith(ext2)]
    del file_list

    # Dicts mapping stem to filename
    ext1_dict = dict((os.path.basename(_).rsplit(".", 2)[0], _) for _ in ext1_list)
    ext2_dict = dict((os.path.basename(_).rsplit(".", 2)[0], _) for _ in ext2_list)

    # This could happen if have same filename used in different folders:
    if len(ext1_dict) < len(ext1_list):
        sys.exit("ERROR: Duplicate *.%s file names" % ext1)
    if len(ext2_dict) < len(ext2_list):
        sys.exit("ERROR: Duplicate *.%s file names" % ext2)
    del ext1_list, ext2_list

    input_list = []
    for stem in ext1_dict:
        if stem in ext2_dict:
            input_list.append((ext1_dict[stem], ext2_dict[stem]))
        else:
            # Acceptable in motivating use case where on a given plate
            # only some of the samples would be known positive controls:
            sys.stderr.write(
                "WARNING: Have %s but missing %s%s\n" % (ext1_dict[stem], stem, ext2)
            )
    # TODO: Check for XXX.ext2 without XXX.ext1 here?
    del ext1_dict, ext2_dict

    return input_list
