# Copyright 2018-2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

"""Helper functions for THAPB-PICT code."""

import os
import hashlib
import subprocess
import sys
import time

from collections import Counter

from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.SeqIO.FastaIO import SimpleFastaParser


def sample_sort(sample_names):
    """Sort sample names like a human.

    Our samples have occasionally used underscores and minus signs inconsistently
    in sample names as field separators or space substitutions. Therefore simple
    ASCII sorting can give surprises such as not grouping by prefix (since the
    underscore is sorted after the digits and letter):

        >>> sorted(["N01-a", "N01_b", "N01 c", "N011-a"])
        ['N01 c', 'N01-a', 'N011-a', 'N01_b']

    We specifically want "_" (ASCII , after the letters) to sort like " "
    or "-" (ASCII 32 or 45, both before the digits and letters). In case
    any samples are using plus/minus, will map underscore and space to
    the minus sign for sorting.

        >>> sample_sort(["N01-a", "N01_b", "N01 c", "N011-d"])
        ['N01-a', 'N01_b', 'N01 c', 'N011-a']

    """
    # TODO: Clever sorting of e.g. A1, A2, ..., A10, A11, A99 or
    # e.g. "Sample 1", "Sample 2", ... "Sample 10", "Sample 11", ...
    # Just sort on the munged sample name
    return sorted(sample_names, key=lambda _: _.replace("_", "-").replace(" ", "-"))


def genus_species_name(genus, species):
    """Return name, genus with species if present.

    Copes with species being None (or empty string).
    """
    # This is a simple function, centralising it for consistency
    assert genus and genus == genus.strip()
    if species:
        assert species == species.strip()
        return "%s %s" % (genus, species)
    else:
        return genus


def species_level(prediction):
    """Is this prediction at species level.

    Returns True for a binomial name (at least one space), False for genus
    only or no prediction.
    """
    assert ";" not in prediction, prediction
    return prediction and " " in prediction


def onebp_substitutions(seq):
    """Generate all 1bp substitutions of the sequence.

    Assumes unambiguous IUPAC codes A, C, G, T only.
    """
    seq = seq.upper()
    variants = set()
    for i in range(len(seq)):
        for s in "ACGT":
            # One base substitions
            variants.add(seq[:i] + s + seq[i + 1 :])
    variants.remove(seq)
    return variants


def onebp_deletions(seq):
    """Generate all variants of sequence with 1bp deletion.

    Assumes unambiguous IUPAC codes A, C, G, T only.
    """
    seq = seq.upper()
    variants = set()
    for i in range(len(seq)):
        # One base deletion
        variants.add(seq[:i] + seq[i + 1 :])
    return variants


def onebp_inserts(seq):
    """Generate all variants of sequence with 1bp insert.

    Assumes unambiguous IUPAC codes A, C, G, T only.
    """
    seq = seq.upper()
    variants = set()
    for i in range(len(seq)):
        for s in "ACGT":
            # One base insertions
            variants.add(seq[:i] + s + seq[i:])
    for s in "ACGT":
        # One base "insertion" at the end
        variants.add(seq + s)
    return variants


def onebp_variants(seq):
    """Generate all 1bp variants of the sequence (substitution, deletion or insertion).

    Assumes unambiguous IUPAC codes A, C, G, T only.
    """
    seq = seq.upper()
    variants = set()
    for i in range(len(seq)):
        # One base deletion
        variants.add(seq[:i] + seq[i + 1 :])
        for s in "ACGT":
            # One base substitions
            variants.add(seq[:i] + s + seq[i + 1 :])
            # One base insertions
            variants.add(seq[:i] + s + seq[i:])
    for s in "ACGT":
        # One base "insertion" at the end
        variants.add(seq + s)
    variants.remove(seq)
    return variants


assert set(onebp_variants("A")) == {
    "",
    "C",
    "G",
    "T",
    "AA",
    "CA",
    "GA",
    "TA",
    "AC",
    "AG",
    "AT",
}

assert set(onebp_variants("AA")) == {
    "A",
    "CA",
    "GA",
    "TA",
    "AC",
    "AG",
    "AT",
    "AAA",
    "CAA",
    "GAA",
    "TAA",
    "ACA",
    "AGA",
    "ATA",
    "AAC",
    "AAG",
    "AAT",
}


def expand_IUPAC_ambiguity_codes(seq):
    """Convert to upper case and iterate over possible unabmigous interpretations.

    This is a crude recursive implementation, intended for use on sequences with
    just a few ambiguity codes in them - it may not scale very well!
    """
    seq = seq.upper()
    ambig_letters = set(seq).difference(["A", "C", "G", "T"])
    for nuc in ambig_letters:
        if nuc not in ambiguous_dna_values:
            raise ValueError("%r is not an IUPAC ambiguous DNA letter" % nuc)
    if not ambig_letters:
        yield seq
    else:
        # Recursive!
        nuc = list(ambig_letters)[0]
        i = seq.index(nuc)  # first appearance of nuc
        before = seq[:i]
        after = seq[i + 1 :]
        for alt in ambiguous_dna_values[nuc]:
            yield from expand_IUPAC_ambiguity_codes(before + alt + after)


# Example has 4 ambiguity codes, each has 2 possible values,
# thus expect 2**4 = 16 interpretations:
assert sorted(expand_IUPAC_ambiguity_codes("GYRGGGACGAAAGTCYYTGC")) == [
    "GCAGGGACGAAAGTCCCTGC",
    "GCAGGGACGAAAGTCCTTGC",
    "GCAGGGACGAAAGTCTCTGC",
    "GCAGGGACGAAAGTCTTTGC",
    "GCGGGGACGAAAGTCCCTGC",
    "GCGGGGACGAAAGTCCTTGC",
    "GCGGGGACGAAAGTCTCTGC",
    "GCGGGGACGAAAGTCTTTGC",
    "GTAGGGACGAAAGTCCCTGC",
    "GTAGGGACGAAAGTCCTTGC",
    "GTAGGGACGAAAGTCTCTGC",
    "GTAGGGACGAAAGTCTTTGC",
    "GTGGGGACGAAAGTCCCTGC",
    "GTGGGGACGAAAGTCCTTGC",
    "GTGGGGACGAAAGTCTCTGC",
    "GTGGGGACGAAAGTCTTTGC",
]


def md5seq(seq):
    """Return MD5 hash of the (upper case) sequence."""
    return hashlib.md5(seq.upper().encode("ascii")).hexdigest()


def cmd_as_string(cmd):
    """Express a list command as a suitably quoted string.

    Intended for using in debugging or error messages.
    """
    if isinstance(cmd, list):
        # Quote any entries with spaces
        return " ".join('"%s"' % _ if " " in _ else _ for _ in cmd)
    else:
        return cmd


def run(cmd, debug=False, attempts=1):
    """Run a command via subprocess, abort if fails."""
    for i in range(attempts):
        if debug:
            if attempts:
                sys.stderr.write(
                    "Attempt %i of %i calling command: %s\n"
                    % (i + 1, attempts, cmd_as_string(cmd))
                )
            else:
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
            if i + 1 < attempts:
                sys.stderr.write(
                    "WARNING: Attempt %i of %i failed with return code %i, cmd:\n%s\n"
                    % (i + 1, attempts, e.returncode, cmd_as_string(cmd))
                )
                time.sleep(min(5, i + 1))
            else:
                if debug:
                    # Used universal_newlines=True above so that this just works
                    # (equivalent to text=True in Python 3.7 onwards):
                    sys.stdout.flush()
                    sys.stdout.write(e.stdout)
                    sys.stdout.flush()
                    sys.stderr.write(e.stderr)
                sys.exit(
                    "ERROR: Attempt %i of %i failed with return code %i, cmd:\n%s\n"
                    % (i + 1, attempts, e.returncode, cmd_as_string(cmd))
                )


def abundance_from_read_name(text, debug=False):
    """Extract abundance from SWARM style read name.

    >>> abundance_from_read_name("9e8f051c64c2b9cc3b6fcb27559418ca_988")
    988

    If fails, will return one.
    """
    try:
        return int(text.rsplit("_", 1)[1])
    except (ValueError, IndexError):
        if debug:
            sys.stderr.write("WARNING: No abundance suffix in %r\n" % text)
        return 1


def split_read_name_abundance(text, debug=False):
    """Split SWARM style read name into prefix and abundance.

    >>> abundance_from_read_name("9e8f051c64c2b9cc3b6fcb27559418ca_988")
    '9e8f051c64c2b9cc3b6fcb27559418ca', 988

    If fails to detect the abundance, will return the original text
    as the prefix with an abundance of 1.
    """
    try:
        prefix, abundance = text.rsplit("_", 1)
        return prefix, int(abundance)
    except (ValueError, IndexError):
        if debug:
            sys.stderr.write("WARNING: No abundance suffix in %r\n" % text)
        return text, 1


def abundance_values_in_fasta(fasta_file):
    """Return total and maximum abundance encoded in read names."""
    total_a = 0
    max_a = Counter()
    with open(fasta_file) as handle:
        for title, _ in SimpleFastaParser(handle):
            a = abundance_from_read_name(title.split(None, 1)[0])
            hmm = title.split(None, 1)[1].strip()
            max_a[hmm] = max(a, max_a[hmm])
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
                    "ERROR: Specified filename %r does not have expected extension %r."
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
    ext1_dict = {os.path.basename(_)[: -len(ext1)]: _ for _ in ext1_list}
    ext2_dict = {os.path.basename(_)[: -len(ext2)]: _ for _ in ext2_list}

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
        elif debug:
            # Acceptable in motivating use case where on a given plate
            # only some of the samples would be known positive controls:
            sys.stderr.write(
                "WARNING: Have %s but missing %s%s\n" % (ext1_dict[stem], stem, ext2)
            )
    # TODO: Check for XXX.ext2 without XXX.ext1 here?
    del ext1_dict, ext2_dict

    return input_list


def parse_species_list_from_tsv(tabular_file):
    """Extract species list from TSV header line (ignores genus only)."""
    with open(tabular_file) as handle:
        line = handle.readline()
    if not line.startswith("#") or line.count("\t") != 3:
        sys.exit(
            "ERROR: %s does not have 4 column TSV header:\n%s" % (tabular_file, line)
        )
    parts = line.rstrip("\n").split("\t")
    if (
        parts[0] != "#sequence-name"
        or parts[1] != "taxid"
        or not parts[2].startswith("genus-species")
        or parts[3] != "note"
    ):
        sys.stderr.write("%r\n" % parts)
        sys.exit(
            "ERROR: %s does not have expected 4 column TSV headers "
            "(sequence-name, taxid, genus-species:..., note):\n%s"
            % (tabular_file, line)
        )
    if not parts[2].startswith("genus-species:"):
        sys.exit(
            "ERROR: %s does not have species list in genus-species column header:\n%s"
            % (tabular_file, line)
        )
    return [_ for _ in parts[2][14:].split(";") if species_level(_)]


def parse_species_tsv(tabular_file, min_abundance=0, req_species_level=False):
    """Parse file of species assignments/predictions by sequence."""
    with open(tabular_file) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            if line.count("\t") != 3:
                sys.exit(
                    "ERROR: %s is not 4 column TSV "
                    "(name, taxid, genus-species, note):\n%s" % (tabular_file, line)
                )
            name, taxid, genus_species, _ = line.split("\t", 3)
            if name == "*":
                raise ValueError("Wildcard species name found")
            if min_abundance > 1 and abundance_from_read_name(name) < min_abundance:
                continue
            if req_species_level:
                assert taxid.count(";") == genus_species.count(";"), line
                wanted = [
                    (t, s)
                    for (t, s) in zip(taxid.split(";"), genus_species.split(";"))
                    if species_level(s)
                ]
                taxid = ";".join(t for (t, s) in wanted)
                genus_species = ";".join(s for (t, s) in wanted)
            yield name, taxid, genus_species


def load_metadata(
    metadata_file,
    metadata_cols,
    metadata_name_row=1,
    metadata_index=0,
    metadata_index_sep=";",
    sequenced_samples=None,
    metadata_sort=True,
    debug=False,
):
    """Load specified metadata as several lists.

    The columns argument should be a string like "1,3,5" - a comma
    separated list of columns to output. The column numbers are assumed
    to be one-based as provided by the command line user.

    The name row indicates which row in the table contains the names
    or descriptions of the metadata columns (one-based).

    The index column is assumed to contain one or more sequenced sample
    names separated by the character specified (default is semi-colon).
    This one-to-many mapping reflecting that a single field sample could
    be sequenced more than once (e.g. technical replicates). These sample
    names are matched against the file name stems, see function find_metadata.

    if metadata_sort=True, then the table rows are sorted based on the
    requested metadata - otherwise the table row order is preserved.

    Returns:
     - list of field sample metadata values (each a list of N values)
     - matching list of associated sequenced sample(s),
     - list of the N field names,
     - matching default value (list of N empty strings).
     - samples without metadata

    Without labeling one of the metadata values as a field-sample ID,
    we cannot index the return values - instead expect caller to use
    the row number.

    Optional argument sequenced_samples should be a set or list of
    sample names which will be cross-checked against the metadata_index
    column. Samples not in the metadata file are one of the return values
    and will generate warning message. The other way round gives a
    missing file warning too (but they are left in the returned data).

    """
    # TODO - Accept Excel style A, ..., Z, AA, ... column names?

    if debug:
        sys.stderr.write(
            "DEBUG: Loading metadata from %r, column specification %r, "
            "field names from row %r\n"
            % (metadata_file, metadata_cols, metadata_name_row)
        )
        if sequenced_samples is not None:
            sys.stderr.write(
                "DEBUG: Have %i sequenced samples\n" % len(sequenced_samples)
            )

    if not metadata_file or not metadata_cols:
        if debug:
            sys.stderr.write("DEBUG: Not loading any metadata\n")
        return [], [], [], [], sequenced_samples

    try:
        value_cols = [int(_) - 1 for _ in metadata_cols.split(",")]
    except ValueError:
        sys.exit(
            "ERROR: Output metadata columns should be a comma separated list "
            "of positive integers, not %r." % metadata_cols
        )
    if min(value_cols) < 0:
        sys.exit("ERROR: Invalid metadata output column, should all be positive.")
    if metadata_index:
        sample_col = int(metadata_index) - 1
        if sample_col < 0:
            sys.exit(
                "ERROR: Invalid metadata index column, should be positive, not %r."
                % metadata_index
            )
    else:
        sample_col = value_cols[0]  # Default is first output column
    if debug:
        sys.stderr.write(
            "DEBUG: Matching sample names to metadata column %i\n" % (sample_col + 1)
        )

    names = [""] * len(value_cols)  # default
    meta = {}
    default = [""] * len(value_cols)
    try:
        # Try with default encoding first...
        with open(metadata_file) as handle:
            lines = list(handle)
    except UnicodeDecodeError:
        # Automatically try latin1, which seems to be
        # default when save tab-delimited from macOS Excel
        with open(metadata_file, encoding="latin1") as handle:
            lines = list(handle)

    if metadata_name_row:
        line = lines[metadata_name_row - 1]
        if line.startswith("#"):
            line = line[1:]
        parts = line.rstrip("\n").split("\t")
        names = [parts[_].strip() for _ in value_cols]
        if debug:
            sys.stderr.write(
                "DEBUG: Row %i gave metadata field names: %r\n"
                % (metadata_name_row, names)
            )

    # Remove header lines,
    lines = [_ for _ in lines[metadata_name_row:] if not _.startswith("#")]

    # Break up line into fields
    lines = [_.rstrip("\n").split("\t") for _ in lines]

    # Select columns of interest
    meta_plus_idx = [[_[i].strip() for i in value_cols + [sample_col]] for _ in lines]

    # Remove blanks
    meta_plus_idx = [_ for _ in meta_plus_idx if any(_)]
    del lines

    # Sort on metadata if requested
    if metadata_sort:
        meta_plus_idx.sort()

    # Select desired columns,
    meta = [_[:-1] for _ in meta_plus_idx]
    index = [[s.strip() for s in _[-1].split(";") if s.strip()] for _ in meta_plus_idx]
    del meta_plus_idx

    back = {}
    for i, samples in enumerate(index):
        for sample in samples:
            try:
                back[sample].add(i)
            except KeyError:
                back[sample] = {i}
    if debug:
        sys.stderr.write(
            "DEBUG: Loaded metadata for %i field samples, %i sequenced samples\n"
            % (len(meta), len(back))
        )

    missing_meta = []
    if sequenced_samples is not None and set(back) != set(sequenced_samples):
        # Using list comprehensions to preserve any meaningful order
        missing_files = [_ for _ in back if _ not in sequenced_samples]
        if missing_files:
            sys.stderr.write(
                "WARNING: Missing %i sequenced samples listed in metadata, %s (etc)\n"
                % (len(missing_files), missing_files[0])
            )
        missing_meta = [_ for _ in sequenced_samples if _ not in back]
        if missing_meta:
            sys.stderr.write(
                "WARNING: %s sequenced samples without metadata, %s (etc)\n"
                % (len(missing_meta), missing_meta[0])
            )

    bad = sample_sort(sample for sample in back if len(back[sample]) > 1)
    if bad:
        print(bad)
        sys.exit(
            "ERROR: Duplicated metadata for %i samples, %s (etc)" % (len(bad), bad[1])
        )
    del back

    return meta, index, names, default, missing_meta
