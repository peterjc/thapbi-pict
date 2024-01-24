# Copyright 2018-2023 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Helper functions for THAPB-PICT code."""
import gzip
import hashlib
import os
import subprocess
import sys
import time
from collections import Counter
from keyword import iskeyword

from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.SeqIO.FastaIO import SimpleFastaParser


KMER_LENGTH = 31


def valid_marker_name(text):
    """Check the proposed string valid for use as a marker name.

    Want to be able to use the string for file or directory names, and also
    column names etc in reports. At very least should reject whitespace, line
    breaks, and slashes.

    Also rejecting all digits, as might want to accept integers as argument
    (e.g. cluster array job mapping task numbers to marker numbers).

    Also rejecting the underscore, as may want to use it as a field separator
    in sequence names (e.g. ``marker_md5_abundance``), and full stop as may
    use it as a field separator in filenames.

    May want to relax this later, thus defining this central function.
    """
    return (
        text.replace("-", "").isalnum()
        and not iskeyword(text)
        and not text.isdigit()
        and text != "pooled"
    )


def primer_clean(primer):
    """Handle non-IUPAC entries in primers, maps I for inosine to N.

    >>> primer_clean("I")
    'N'

    Inosine is found naturally at the wobble position of tRNA, and can match
    any base. Structurally similar to guanine (G), it preferentially binds
    cytosine (C). It sometimes used in primer design (Ben-Dov et al, 2006),
    where degeneracy N would give similar results.
    """
    return primer.upper().replace("I", "N")


def reject_species_name(species):
    """Reject species names like 'environmental samples' or 'uncultured ...'."""
    return (
        not species
        or species.split(None, 1)[0] in ("unclassified", "uncultured", "unidentified")
        or species == "environmental samples"
    )


def genus_species_name(genus, species):
    """Return name, genus with species if present.

    Copes with species being None (or empty string).
    """
    # This is a simple function, centralising it for consistency
    assert genus and genus == genus.strip(), repr(genus)
    if species:
        assert species == species.strip(), repr(species)
        return f"{genus} {species}"
    else:
        return genus


def genus_species_split(name):
    """Return (genus, species) splitting on first space.

    If there are no spaces, returns (name, '') instead.
    """
    if " " in name:
        return name.split(" ", 1)
    else:
        return name, ""


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
            # One base substitutions
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


def expand_IUPAC_ambiguity_codes(seq):
    """Convert to upper case and iterate over possible unabmigous interpretations.

    This is a crude recursive implementation, intended for use on sequences with
    just a few ambiguity codes in them - it may not scale very well!
    """
    seq = seq.upper()
    ambig_letters = set(seq).difference(["A", "C", "G", "T"])
    for nuc in ambig_letters:
        if nuc not in ambiguous_dna_values:
            msg = f"{nuc!r} is not an IUPAC ambiguous DNA letter"
            raise ValueError(msg)
    if not ambig_letters:
        yield seq
    else:
        # Recursive!
        nuc = next(iter(ambig_letters))
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
    """Return MD5 32-letter hex digest of the (upper case) sequence.

    >>> md5seq("ACGT")
    'f1f8f4bf413b16ad135722aa4591043e'

    """
    return hashlib.md5(seq.upper().encode("ascii")).hexdigest()


def md5_hexdigest(filename, chunk_size=1024):
    """Return the MD5 hex-digest of the given file."""
    hash_md5 = hashlib.md5()
    with open(filename, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                # EOF
                break
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def cmd_as_string(cmd):
    """Express a list command as a suitably quoted string.

    Intended for using in debugging or error messages.
    """
    if isinstance(cmd, list):
        # Quote any entries with spaces or semi-colons
        return " ".join('"%s"' % _ if (" " in _ or ";" in _) else _ for _ in cmd)
    else:
        return cmd


def run(cmd, debug=False, attempts=1):
    """Run a command via subprocess, abort if fails.

    Returns a subprocess.CompletedProcess object.
    """
    for i in range(attempts):
        if debug:
            if attempts > 1:
                sys.stderr.write(
                    f"Attempt {i + 1} of {attempts} calling command:"
                    f" {cmd_as_string(cmd)}\n"
                )
            else:
                sys.stderr.write(f"Calling command: {cmd_as_string(cmd)}\n")
        try:
            # On Python 3.7 onwards, could use capture_output=True
            # rather than stdout=PIPE and stderr=PIPE
            if isinstance(cmd, list):
                return subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=True,
                )
            else:
                return subprocess.run(
                    cmd,
                    shell=True,
                    capture_output=True,
                    text=True,
                    check=True,
                )
        except subprocess.CalledProcessError as e:
            if i + 1 < attempts:
                sys.stderr.write(
                    f"WARNING: Attempt {i + 1} of {attempts} failed"
                    f" with return code {e.returncode}, cmd:\n{cmd_as_string(cmd)}\n"
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
                    f"ERROR: Attempt {i + 1} of {attempts} failed"
                    f" with return code {e.returncode}, cmd:\n{cmd_as_string(cmd)}\n"
                )
    return None


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
            sys.stderr.write(f"WARNING: No abundance suffix in {text!r}\n")
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
            sys.stderr.write(f"WARNING: No abundance suffix in {text!r}\n")
        return text, 1


def abundance_values_in_fasta(fasta_file, gzipped=False):
    """Return unique count, total abundance, and maximum abundances by spike-in."""
    unique = 0
    total = 0
    max_a = Counter()
    if gzipped:
        handle = gzip.open(fasta_file, "rt")
    else:
        handle = open(fasta_file)
    with handle:
        for title, _ in SimpleFastaParser(handle):
            unique += 1
            a = abundance_from_read_name(title.split(None, 1)[0])
            total += a
            try:
                hmm = title.split(None, 1)[1].strip()
            except IndexError:
                hmm = ""  # prepared with no HMM
            max_a[hmm] = max(a, max_a[hmm])
    return unique, total, max_a


def abundance_filter_fasta(input_fasta, output_fasta, min_abundance):
    """Apply a minimum abundance filter to a FASTA file."""
    with open(input_fasta) as in_handle:
        with open(output_fasta, "w") as out_handle:
            for title, seq in SimpleFastaParser(in_handle):
                a = abundance_from_read_name(title.split(None, 1)[0])
                if a < min_abundance:
                    continue
                out_handle.write(f">{title}\n{seq}\n")


def file_to_sample_name(filename):
    """Given filename (with or without a directory name), return sample name only.

    i.e. XXX.fasta, XXX.fastq.gz, XXX.method.tsv --> XXX
    """
    if filename.endswith(".fasta"):
        return os.path.basename(filename).rsplit(".", 1)[0]
    elif filename.endswith((".fasta.gz", ".tsv", ".fastq.gz")):
        return os.path.basename(filename).rsplit(".", 2)[0]
    else:
        msg = f"Invalid file_to_sample_name arg: {filename}"
        raise ValueError(msg)


def find_requested_files(
    filenames_or_folders, ext=".fasta", ignore_prefixes=None, debug=False
):
    """Interpret a list of filenames and/or foldernames.

    The extensions argument can be a tuple.
    """
    assert ignore_prefixes  # DEBUG - TODO - remove this
    answer = []
    for x in filenames_or_folders:
        if os.path.isdir(x):
            if debug:
                sys.stderr.write(f"DEBUG: Walking directory {x}\n")
            for root, _, files in os.walk(x, followlinks=True):
                for f in files:
                    if f.endswith(ext):
                        # Check not a directory?
                        if ignore_prefixes and f.startswith(ignore_prefixes):
                            if debug:
                                sys.stderr.write(
                                    f"DEBUG: Ignoring {os.path.join(root, f)}\n"
                                )
                            continue
                        answer.append(os.path.join(root, f))
        elif os.path.isfile(x):
            if x.endswith(ext):
                if ignore_prefixes and x.startswith(ignore_prefixes):
                    if debug:
                        sys.stderr.write(f"DEBUG: Due to prefix ignoring {x}\n")
                    continue
                answer.append(x)
            elif debug:
                sys.stderr.write(
                    f"DEBUG: Looking for extension {ext} so ignoring {x}\n"
                )
        else:
            sys.exit(f"ERROR: {x!r} is not a file or a directory\n")
    # Warn if there were duplicates?
    return sorted(set(answer))


def find_paired_files(
    filenames_or_folders,
    ext1,
    ext2,
    ignore_prefixes=None,
    debug=False,
    strict=False,
):
    """Interpret a list of filenames and/or foldernames to find pairs.

    Looks for paired files named XXX.ext1 and XXX.ext2 which can be
    in different directories - duplicated filenames (in different
    directories) are considered to be an error.

    Having XXX.ext1 without XXX.ext2 is an error in strict mode, or a warning
    in debug mode, otherwise silently ignored.

    Having XXX.ext2 without XXX.ext1 is silently ignored.

    The arguments ext1 and ext2 should include the leading dot.
    """
    assert ext1.startswith("."), ext1
    assert ext2.startswith("."), ext2
    file_list = find_requested_files(
        filenames_or_folders,
        ext=(ext1, ext2),
        ignore_prefixes=ignore_prefixes,
        debug=debug,
    )
    ext1_list = [_ for _ in file_list if _.endswith(ext1)]
    ext2_list = [_ for _ in file_list if _.endswith(ext2)]
    del file_list

    # Dicts mapping stem to filename
    ext1_dict = {os.path.basename(_)[: -len(ext1)]: _ for _ in ext1_list}
    ext2_dict = {os.path.basename(_)[: -len(ext2)]: _ for _ in ext2_list}

    # This could happen if have same filename used in different folders:
    if len(ext1_dict) < len(ext1_list):
        sys.exit(f"ERROR: Duplicate *{ext1} file names")
    if len(ext2_dict) < len(ext2_list):
        sys.exit(f"ERROR: Duplicate *{ext2} file names")
    del ext1_list, ext2_list

    input_list = []
    for stem in ext1_dict:
        if stem in ext2_dict:
            input_list.append((ext1_dict[stem], ext2_dict[stem]))
        elif strict:
            sys.exit(f"ERROR: {ext1_dict[stem]} without {stem}{ext2}")
        elif debug:
            # Acceptable in motivating use case where on a given plate
            # only some of the samples would be known positive controls:
            sys.stderr.write(
                f"WARNING: Ignoring {ext1_dict[stem]} as missing {stem}{ext2}\n"
            )
    # TODO: Check for XXX.ext2 without XXX.ext1 here?
    del ext1_dict, ext2_dict

    return input_list


def export_sample_biom(output_file, seqs, seq_meta, sample_meta, counts, gzipped=True):
    """Export a sequence vs samples counts BIOM table, with metadata.

    Similar to the export_sample_tsv file (our TSV output), expects same
    arguments as loaded from one of our TSV files via he parse_sample_tsv
    function.

    Will save a BIOM v2 HDF5 file if possible and return True. If output
    fails, returns False.
    """
    from thapbi_pict import __version__

    try:
        from biom.table import Table
        from biom.util import biom_open
    except ImportError:
        return False

    biom_table = Table(
        # BIOM want counts indexed by sequence and sample integer index:
        {
            (list(seq_meta).index((marker, md5)), list(sample_meta).index(sample)): v
            for (marker, md5, sample), v in counts.items()
        },
        # BIOM wants single string names for sequences
        [f"{marker}/{md5}" for (marker, md5) in seq_meta],
        list(sample_meta),
        # Add the sequence itself to the metadata dict for BIOM export:
        [dict([("Sequence", seqs[k]), *list(v.items())]) for k, v in seq_meta.items()],
        sample_meta.values(),
        # Required attribute in BIOM format:
        type="OTU table",
    )
    del seqs, seq_meta, sample_meta, counts

    tag = "THAPBI PICT " + __version__
    # TODO - override refault date of now for reproducibility?

    # Is JSON output useful? BIOM v1 format:
    # with open(output_file, "w") as handle:
    #     biom_table.to_json(generated_by=tag, direct_io=handle)

    with biom_open(output_file, "w") as handle:
        biom_table.to_hdf5(handle, generated_by=tag, compress=gzipped)
    return True


def export_sample_tsv(output_file, seqs, seq_meta, sample_meta, counts, gzipped=False):
    """Export a sequence vs sample counts TSV table, with metadata.

    The TSV file ought to be readable by the parse_sample_tsv function, and
    is first generated in our pipeline by the sample-tally command, and then
    extended by the classify command to add taxonomic sequence metadata.

    If the output tabular file argument is "-", it writes to stdout (not
    supported with gzipped mode).

    With no sequence metadata this should be accepted as a TSV BIOM file.
    """
    if (
        not isinstance(seqs, dict)
        or not isinstance(seq_meta, dict)
        or not isinstance(sample_meta, dict)
        or not isinstance(counts, dict)
    ):
        msg = "Expect dictionaries are arguments"
        raise TypeError(msg)
    if output_file == "-":
        if gzipped:
            msg = "Does not support gzipped output to stdout."
            raise ValueError(msg)
        out_handle = sys.stdout
    elif gzipped:
        out_handle = gzip.open(output_file, "wt")
    else:
        out_handle = open(output_file, "w")

    # Write a header row for any sample metadata
    if sample_meta:
        sample_fields = None
        for sample, values in sample_meta.items():
            if sample_fields is None:
                sample_fields = list(values.keys())
            elif sample_fields != list(values.keys()):
                msg = f"Inconsistent sample metadata keys in {sample}"
                raise ValueError(msg)
    else:
        sample_fields = []
    if seq_meta:
        seq_fields = None
        for key, values in seq_meta.items():
            if seq_fields is None:
                seq_fields = list(values.keys())
            elif seq_fields != list(values.keys()):
                msg = f"Inconsistent seq metadata keys in {key}"
                raise ValueError(msg)
    else:
        seq_fields = []

    for stat in sample_fields:
        # Using "-" as missing value to match default in summary reports
        values = [sample_meta[sample][stat] for sample in sample_meta]
        out_handle.write(
            "\t".join(
                ["#" + stat]
                + ["-" if _ is None else str(_) for _ in values]
                + [""] * len(seq_fields)
                + ["\n"]
            )
        )

    samples = list(sample_meta)
    out_handle.write(
        "\t".join(["#Marker/MD5_abundance", *samples, "Sequence", *seq_fields]) + "\n"
    )

    for (marker, idn), seq in seqs.items():
        values = [counts.get((marker, idn, sample), 0) for sample in samples]
        out_handle.write(
            "\t".join(
                [f"{marker}/{idn}_{sum(values)}"]
                + [str(_) for _ in values]
                + [seq]
                + [str(seq_meta[marker, idn][_]) for _ in seq_fields]
            )
            + "\n"
        )

    if output_file != "-":
        out_handle.close()


def parse_sample_tsv(tabular_file, min_abundance=0, debug=False, force_upper=True):
    """Parse file of sample abundances and sequence (etc).

    Optional argument min_abundance is applied to the per sequence per sample
    values (i.e. the matrix elements, not the row/column totals).

    Columns are:
    * Sequence label, <marker>/<identifier>_<abundance>
    * Column per sample giving the sequence count
    * Sequence itself
    * Optional additional columns for sequence metadata (e.g. chimera flags)

    Supports optional sample metadata header too as # prefixed header lines.

    Returns dictionaries of:
    * Sequence keyed on [<marker>, <identitifer>], string
    * Sequence metadata keyed [<marker>, <identitifer>], dict of key:value pairs
    * Sample metadata keyed on [<sample>], dict of key:value pairs
    * Counts keyed on 3-tuple [<marker>, <identifier>, <sample>], integer
    """
    header_lines = []
    samples = []
    counts = {}
    seqs = {}
    seq_meta_keys = None
    seq_meta = {}
    seq_col = None
    with open(tabular_file) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == "#Marker/MD5_abundance":
                seq_col = parts.index("Sequence")
                samples = parts[1:seq_col]
                seq_meta_keys = parts[seq_col + 1 :]
                if debug:
                    sys.stderr.write(
                        f"DEBUG: {len(samples)} samples in {tabular_file}\n"
                    )
            elif line.startswith("#"):
                header_lines.append(parts)
                if debug:
                    sys.stderr.write(f"DEBUG: Header {parts[0]} in {tabular_file}\n")
            elif samples:
                marker, idn = parts[0].split("/")
                idn = idn.rsplit("_")[0]  # drop the total count
                above_threshold = False
                assert seq_col, "Error: Did not find 'Sequence' column"
                for sample, value in zip(samples, parts[1:seq_col]):
                    if value == "0":  # Don't bother with int("0")
                        continue  # Don't store blank values!
                    try:
                        value = int(value)
                    except ValueError:
                        msg = f"ERROR: Non-integer count for {parts[0]} vs {sample}"
                        raise ValueError(msg) from None
                    if min_abundance <= value:
                        counts[marker, idn, sample] = value
                        above_threshold = True
                if above_threshold:
                    if force_upper:
                        seqs[marker, idn] = parts[seq_col].upper()
                    else:
                        seqs[marker, idn] = parts[seq_col]
                    if seq_meta_keys:
                        seq_meta[marker, idn] = dict(
                            zip(seq_meta_keys, parts[seq_col + 1 :])
                        )
            else:
                msg = (
                    r"ERROR: Missing #Marker/MD5_abundance(tab)...(tab)Sequence\n"
                    f" line in {tabular_file}"
                )
                raise ValueError(msg)
    # TODO: Turn counts into an array?
    sample_headers = {sample: {} for sample in samples}
    for parts in header_lines:
        name = parts[0][1:]  # Drop the leading "#"
        for sample, value in zip(samples, parts[1:seq_col]):
            sample_headers[sample][name] = value
    return seqs, seq_meta, sample_headers, counts


def parse_species_tsv(
    tabular_file, min_abundance=0, req_species_level=False, allow_wildcard=False
):
    """Parse file of species assignments/predictions by sequence.

    Yields tuples of marker name (from the file header line), sequence name,
    taxid, and genus_species.
    """
    marker = None
    tally_mode = False
    with open(tabular_file) as handle:
        for line in handle:
            if line.startswith("#"):
                parts = line[1:].strip("\n").split("\t")
                if parts[0].endswith("/sequence-name") and parts[1] == "taxid":
                    marker = parts[0][:-14]
                if (
                    "Sequence" in parts
                    and "taxid" in parts
                    and "genus-species" in parts
                ):
                    tally_mode = (parts.index("taxid"), parts.index("genus-species"))
                    allow_wildcard = False
                continue
            if tally_mode:
                parts = line.strip("\n").split("\t")
                marker2, name = parts[0].split("/")
                if marker is None:
                    marker = marker2
                elif marker != marker2:
                    msg = "Mixed markers in {tabular_file}"
                    raise ValueError(msg)
                taxid, genus_species = (parts[_] for _ in tally_mode)
            elif line.count("\t") == 2:
                name, taxid, genus_species = line.rstrip("\n").split("\t", 3)
            elif line.count("\t") == 3:
                name, taxid, genus_species, _ = line.split("\t", 3)
            else:
                sys.exit(
                    f"ERROR: {tabular_file} is not 3 or 4 column TSV"
                    f" (name, taxid, genus-species, optional note):\n{line}"
                )
            if name == "*" and not allow_wildcard:
                msg = "Wildcard species name found"
                raise ValueError(msg)
            if (
                min_abundance > 1
                and abundance_from_read_name(name) < min_abundance
                and name != "*"
            ):
                continue
            if req_species_level:
                if taxid and taxid != "0":
                    assert taxid.count(";") == genus_species.count(";"), line
                    wanted = [
                        (t, s)
                        for (t, s) in zip(taxid.split(";"), genus_species.split(";"))
                        if species_level(s)
                    ]
                    taxid = ";".join(t for (t, s) in wanted)
                    genus_species = ";".join(s for (t, s) in wanted)
                else:
                    taxid = ""
                    genus_species = ";".join(
                        s for s in genus_species.split(";") if species_level(s)
                    )
            yield marker, name, taxid, genus_species


def load_metadata(
    metadata_file,
    metadata_encoding,
    metadata_cols,
    metadata_groups=None,
    metadata_name_row=1,
    metadata_index=0,
    metadata_index_sep=";",
    ignore_prefixes=("Undetermined",),
    debug=False,
):
    """Load specified metadata as several lists.

    The encoding argument can be None or "", meaning use the default.

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

    The metadata table rows are sorted based on the requested colunms.

    Return values:

    - Dict mapping FASTQ stems to metadata tuples
    - Ordered dict mapping metadata tuples to lists of FASTQ stems
    - list of the N field names
    - Color grouping offset into the N values
    """
    # TODO - Accept Excel style A, ..., Z, AA, ... column names?

    if not metadata_file or not metadata_cols:
        if debug:
            sys.stderr.write("DEBUG: Not loading any metadata\n")
        return {}, {}, [], 0

    if metadata_groups and not metadata_cols:
        sys.exit("ERROR: Using -g / --metagroups requires -c / --metacols")

    if debug:
        sys.stderr.write(
            f"DEBUG: Loading metadata from {metadata_file!r},"
            f" column specification {metadata_cols!r},"
            f" field names from row {metadata_name_row!r}\n"
        )

    try:
        value_cols = [int(_) - 1 for _ in metadata_cols.split(",")]
    except ValueError:
        sys.exit(
            "ERROR: Output metadata columns should be a comma separated list"
            f" of positive integers, not {metadata_cols!r}."
        )
    if min(value_cols) < 0:
        sys.exit("ERROR: Invalid metadata output column, should all be positive.")
    if metadata_index:
        sample_col = int(metadata_index) - 1
        if sample_col < 0:
            sys.exit(
                "ERROR: Invalid metadata index column, should be positive,"
                f" not {metadata_index!r}."
            )
    else:
        sample_col = value_cols[0]  # Default is first output column
    if debug:
        sys.stderr.write(
            f"DEBUG: Matching sample names to metadata column {sample_col + 1}\n"
        )

    if metadata_groups:
        try:
            group_col = int(metadata_groups) - 1
        except ValueError:
            sys.exit(
                "ERROR: Invalid metadata group column, should be positive or 0,"
                f" not {metadata_groups!r}."
            )
        if group_col not in value_cols:
            sys.exit(
                f"ERROR: Metadata group column {metadata_groups}"
                " not included in reported metadata.\n"
            )
        group_col = value_cols.index(group_col)  # i.e. which of requested columns
    else:
        group_col = 0  # Default is the first requested column

    def dequote_line(fields_in_line):
        return [
            _[1:-1] if _.count('"') == 2 and _[0] == '"' and _[-1] == '"' else _
            for _ in fields_in_line
        ]

    assert dequote_line(["A", "A is B", '"A is also C"', "Not D"]) == [
        "A",
        "A is B",
        "A is also C",
        "Not D",
    ]

    names = [""] * len(value_cols)  # default
    meta = {}
    if not metadata_encoding:
        metadata_encoding = None
    try:
        # Use metadata_encoding=None to try with default encoding:
        with open(metadata_file, encoding=metadata_encoding) as handle:
            lines = list(handle)
    except UnicodeDecodeError:
        if metadata_encoding:
            sys.exit(
                f"ERROR: Specified encoding {metadata_encoding!r}"
                f" incompatible with {metadata_file}"
            )
        else:
            sys.exit(f"ERROR: Default encoding incompatible with {metadata_file}")

    if metadata_name_row:
        line = lines[metadata_name_row - 1]
        if line.startswith("#"):
            line = line[1:]
        parts = dequote_line(line.rstrip("\n").split("\t"))
        if len(parts) < max(value_cols) + 1:
            sys.exit("ERROR: Not enough columns in metadata name row")
        names = [parts[_].strip() for _ in value_cols]
        if debug:
            sys.stderr.write(
                f"DEBUG: Row {metadata_name_row} gave metadata field names:"
                f" {names!r}\n"
            )
            sys.stderr.write(
                f"DEBUG: Grouping on {names[group_col]} for colour bands\n"
            )

    # Remove header lines, and empty lines
    lines = [
        _ for _ in lines[metadata_name_row:] if _.strip() and not _.startswith("#")
    ]

    # Break up line into fields
    lines = [dequote_line(_.rstrip("\n").split("\t")) for _ in lines]

    for _ in lines:
        if len(_) <= sample_col:
            sys.exit(f"ERROR: Missing sample column {sample_col+1} in metadata")
        if len(_) <= max(value_cols):
            sys.exit(
                f"ERROR: Missing column {max(value_cols)+1} for sample {_[sample_col]}"
            )

    # Select columns of interest
    meta_plus_idx = [[_[i].strip() for i in [*value_cols, sample_col]] for _ in lines]

    # Remove blanks
    meta_plus_idx = [_ for _ in meta_plus_idx if any(_)]
    del lines

    # Sort on metadata
    meta_plus_idx.sort()

    bad = set()
    meta_to_stem = {}
    stem_to_meta = {}
    for meta_and_index in meta_plus_idx:
        meta = tuple(meta_and_index[:-1])
        stems = [_.strip() for _ in meta_and_index[-1].split(metadata_index_sep)]
        stems = [_ for _ in stems if _]  # drop any blanks
        if ignore_prefixes:
            stems = [_ for _ in stems if not _.startswith(ignore_prefixes)]
        for stem in stems:
            if stem in stem_to_meta:
                bad.add(stem)
            stem_to_meta[stem] = meta
        if meta in meta_to_stem:
            # combine this row with previous lines with same metadata columns
            meta_to_stem[meta].extend(stems)
        else:
            meta_to_stem[meta] = stems

    if debug:
        sys.stderr.write(
            f"DEBUG: Loaded {len(meta_to_stem)} metadata entries, and "
            f" {len(stem_to_meta)} sequenced samples\n"
        )
    if bad:
        sys.exit(
            f"ERROR: Duplicated metadata for {len(bad)} samples,"
            f" {sorted(bad)[0]} (etc)"
        )

    return stem_to_meta, meta_to_stem, names, group_col


def color_bands(meta_groups, sample_color_bands, debug=False):
    """Return a list for formats, one for each sample."""
    default = [None] * len(meta_groups)

    min_groups = 3
    max_groups = 0.8 * len(meta_groups)
    if min_groups > max_groups:
        if debug:
            sys.stderr.write("DEBUG: Not enough columns to bother with coloring\n")
        return default

    if len(set(meta_groups)) == 1:
        # All the same not helpful for banding
        if debug:
            sys.stderr.write(
                "DEBUG: All samples had same metadata in color grouping field:"
                f" {meta_groups[0]!r}\n"
            )
        return default

    if max_groups < len(set(meta_groups)):
        if debug:
            sys.stderr.write(
                "DEBUG: Too many coloring groups, trying first word only\n"
            )
        # Attempting heuristic, taking the first word/field will work on schemes like
        # SITE_DATE_NUMBER or SPECIES-SAMPLE etc.
        meta_groups = [
            _.replace("-", " ").replace("_", " ").split(None, 1)[0] if _ else ""
            for _ in meta_groups
        ]
        if len(set(meta_groups)) == 1:
            # That didn't work.
            if debug:
                sys.stderr.write(
                    "DEBUG: Too many coloring groups, but first word only was unique\n"
                )
            return default

    if len(set(meta_groups)) < min_groups or max_groups < len(set(meta_groups)):
        # (Almost) all same or (almost) unique not helpful
        if debug:
            sys.stderr.write(
                f"DEBUG: {len(set(meta_groups))} groups not suitable for coloring:\n"
            )
            sys.stderr.write(";".join(sorted(set(meta_groups))) + "\n")
        return default

    bands = []
    debug_msg = []
    for s, value in enumerate(meta_groups):
        if s == 0:
            bands.append(0)
            debug_msg.append(value)
        elif value == meta_groups[s - 1]:
            # Same
            bands.append(bands[-1])
        else:
            # Different
            if value in meta_groups[:s]:
                # Metadata values are not grouped, can't use for banding
                # (might be able to use with a color key?)
                sys.stderr.write(
                    "WARNING: Metadata not grouped nicely for coloring:"
                    f" {', '.join(debug_msg)} and then {value} (again).\n"
                )
                return default
            bands.append(max(bands) + 1)
            debug_msg.append(value)
    assert len(set(bands)) == max(bands) + 1, bands
    return [sample_color_bands[_ % len(sample_color_bands)] for _ in bands]


def load_fasta_header(fasta_file, gzipped=False):
    """Parse our FASTA hash-comment line header as a dict."""
    answer = {}
    if gzipped:
        handle = gzip.open(fasta_file, "rt")
    else:
        handle = open(fasta_file)
    for line in handle:
        if line.startswith("#") and ":" in line:
            tag, value = line[1:].strip().split(":", 1)
            try:
                value = int(value)
            except ValueError:
                pass
            answer[tag] = value
        elif line.startswith(">"):
            break
        elif not line.strip():
            pass
        else:
            sys.exit(f"ERROR: Unexpected line in headered FASTA file {fasta_file}")
    return answer


def kmers(sequence, k=KMER_LENGTH):
    """Make set of all kmers in the given sequence."""
    return {sequence[i : i + k] for i in range(len(sequence) - k + 1)}


def is_spike_in(sequence, spikes):
    """Return spike-in name if sequence matches, else empty string."""
    for spike_name, spike_seq, spike_kmers in spikes:
        if sequence == spike_seq:
            return spike_name
        # This will not work when len(spike) <~ kmer length
        # (fail gracefully with an impossible to meet value of 10)
        threshold = max((len(spike_seq) - KMER_LENGTH) / 3, 10)
        count = 0
        for i in range(len(sequence) - KMER_LENGTH + 1):
            if sequence[i : i + KMER_LENGTH] in spike_kmers:
                count += 1
                if count >= threshold:
                    return spike_name
    return ""
