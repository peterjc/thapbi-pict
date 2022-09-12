# Copyright 2018-2022 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""Classifying prepared marker sequences using a marker database.

This implements the ``thapbi_pict classify ...`` command.
"""
import os
import shutil
import sys
import tempfile
from collections import Counter
from collections import OrderedDict

from Bio.SeqIO.FastaIO import SimpleFastaParser
from numpy import int8
from rapidfuzz.distance import Levenshtein
from rapidfuzz.process import cdist

from .db_orm import MarkerDef
from .db_orm import MarkerSeq
from .db_orm import SeqSource
from .db_orm import Taxonomy
from .utils import abundance_filter_fasta
from .utils import abundance_from_read_name
from .utils import find_requested_files
from .utils import genus_species_name
from .utils import load_fasta_header
from .utils import run
from .utils import species_level
from .versions import check_rapidfuzz
from .versions import check_tools

MIN_BLAST_COVERAGE = 0.85  # percentage of query length

db_seqs = None  # global dict (seq to genus) for onebp and 1s?g distance classifiers
max_dist_genus = None  # global variable for 1s?g distance classifiers
genus_taxid = {}  # global variable to cache taxids for genus names


def unique_or_separated(values, sep=";"):
    """Return sole element, or a string joining all elements using the separator."""
    if len(set(values)) == 1:
        return values[0]
    else:
        return sep.join(str(_) for _ in values)


def consoliate_and_sort_taxonomy(genus_species_taxid):
    """Remove any redundant entries, returns new sorted list.

    Drops zero taxid entries if has matching non-zero entry.

    Drops genus only entries if have species level entries.
    Note ignoring the TaxID here - would need to know the parent/child
    relationship to confirm the genus we're removing does have species
    level children in the prediction set.
    """
    # First drop any redundant genus-only entries...
    answer = []
    genera_with_species = set()
    with_taxid = set()
    # Reverse sort to put genus only before genus with species
    for genus, species, taxid in sorted(genus_species_taxid, reverse=True):
        assert genus
        if species:
            genera_with_species.add(genus)
        if not species and genus in genera_with_species:
            # Drop this genus only entry as have an entry with species
            continue
        if taxid:
            with_taxid.add((genus, species))
        answer.append((genus, species, taxid))
    # Now drop any redundant zero-taxid entries...
    if with_taxid:
        answer = [
            (genus, species, taxid)
            for (genus, species, taxid) in answer
            if taxid or (genus, species) not in with_taxid
        ]
    # Done, just need to reverse the answer...
    return sorted(answer)


assert consoliate_and_sort_taxonomy([("Genie", "alpha", 101), ("Genie", "", 100)]) == [
    ("Genie", "alpha", 101)
]
assert consoliate_and_sort_taxonomy([("Genie", "", 100), ("Genie", "", 0)]) == [
    ("Genie", "", 100)
]
assert consoliate_and_sort_taxonomy(
    [
        ("Genie", "alpha", 0),
        ("Genie", "", 0),
        ("Genie", "", 100),
        ("Genie", "alpha", 101),
    ]
) == [("Genie", "alpha", 101)]


def taxid_and_sp_lists(taxon_entries):
    """Return semi-colon separated summary of the taxonomy objects from DB.

    Will discard genus level predictions (e.g. 'Phytophthora') if there is a
    species level prediciton within that genus (e.g. 'Phytophthora infestans').

    If there is a single result, returns a tuple of taxid (integer), genus-species,
    and debugging comment (strings).

    If any of the fields has conflicting values, returns two semi-colon separated
    string instead (in the same order so you can match taxid to species, sorting
    on the genus-species string).
    """
    tax = consoliate_and_sort_taxonomy(
        {(t.genus, t.species, t.ncbi_taxid) for t in taxon_entries}
    )
    if not tax:
        return 0, "", "No DB match"
    return (
        unique_or_separated([t[2] for t in tax]),
        unique_or_separated([genus_species_name(t[0], t[1]) for t in tax]),
        "",  # Not very useful to report # of entries
    )


def perfect_match_in_db(session, marker_name, seq, debug=False):
    """Lookup sequence in DB, returns taxid, genus_species, note as tuple.

    If the 100% matches in the DB give multiple species, then taxid and
    genus_species will be semi-colon separated strings.
    """
    assert seq == seq.upper(), seq
    # Now, does this equal any of the marker sequences in our DB?
    return taxid_and_sp_lists(
        session.query(Taxonomy)
        .join(SeqSource)
        .join(MarkerDef, SeqSource.marker_definition)
        .filter(MarkerDef.name == marker_name)
        .join(MarkerSeq)
        .filter(MarkerSeq.sequence == seq)
        .distinct()
    )


def perfect_substr_in_db(session, marker_name, seq, debug=False):
    """Lookup sequence in DB, returns taxid, genus_species, note as tuple.

    If the matches containing the sequence as a substring give multiple species,
    then taxid and genus_species will be semi-colon separated strings.
    """
    assert seq == seq.upper(), seq
    return taxid_and_sp_lists(
        session.query(Taxonomy)
        .join(SeqSource)
        .join(MarkerDef, SeqSource.marker_definition)
        .filter(MarkerDef.name == marker_name)
        .join(MarkerSeq)
        .filter(MarkerSeq.sequence.like("%" + seq + "%"))
        .distinct()
    )


def apply_method_to_file(
    method_fn,
    fasta_file,
    session,
    marker_name,
    read_report,
    min_abundance=0,
    debug=False,
):
    """Call given method on each sequence in the FASTA file."""
    count = 0
    tax_counts = Counter()

    with open(fasta_file) as handle:
        for title, seq in SimpleFastaParser(handle):
            idn = title.split(None, 1)[0]
            abundance = abundance_from_read_name(idn)
            if min_abundance and abundance < min_abundance:
                continue
            count += abundance
            taxid, genus_species, note = method_fn(
                session, marker_name, seq.upper(), debug=debug
            )
            tax_counts[genus_species] += abundance
            if debug:
                read_report.write(f"{idn}\t{str(taxid)}\t{genus_species}\t{note}\n")
            else:
                read_report.write(f"{idn}\t{str(taxid)}\t{genus_species}\n")
    assert count == sum(tax_counts.values())
    return tax_counts


def method_identity(
    fasta_file,
    session,
    marker_name,
    read_report,
    tmp_dir,
    shared_tmp_dir,
    min_abundance=0,
    debug=False,
    cpu=0,
):
    """Classify using perfect identity.

    This is a deliberately simple approach, in part for testing
    purposes. It looks for a perfect identical entry in the database.
    """
    return apply_method_to_file(
        perfect_match_in_db,
        fasta_file,
        session,
        marker_name,
        read_report,
        min_abundance=min_abundance,
        debug=debug,
    )


def method_substr(
    fasta_file,
    session,
    marker_name,
    read_report,
    tmp_dir,
    shared_tmp_dir,
    min_abundance=0,
    debug=False,
    cpu=0,
):
    """Classify using perfect identity including as a sub-string.

    Like the 'identity' method, but allows for a database where the marker
    has not been trimmed, or has been imperfectly trimmed (e.g. primer
    mismatch).
    """
    return apply_method_to_file(
        perfect_substr_in_db,
        fasta_file,
        session,
        marker_name,
        read_report,
        min_abundance=min_abundance,
        debug=debug,
    )


def setup_onebp(session, marker_name, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a set of all the DB marker sequences as upper case strings.

    Also setup dict of any entries in the DB with a single ambiguous base.
    """
    global db_seqs
    global max_dist_genus
    max_dist_genus = 1

    check_rapidfuzz()

    db_seqs = {}
    for source in (
        session.query(SeqSource)
        .join(MarkerDef, SeqSource.marker_definition)
        .filter(MarkerDef.name == marker_name)
    ):
        seq = source.marker_seq.sequence.upper()
        genus = source.taxonomy.genus
        assert seq and genus, source
        try:
            db_seqs[seq].add(genus)
        except KeyError:
            db_seqs[seq] = {genus}

    # Update the genus to taxid cache as needed:
    for genus in db_seqs.values():
        for g in genus:
            if g not in genus_taxid:
                query = (
                    session.query(Taxonomy.ncbi_taxid)
                    .filter(Taxonomy.genus == g)
                    .filter(Taxonomy.species == "")
                    .one_or_none()
                )
                genus_taxid[g] = query.ncbi_taxid if query else 0


def setup_dist2(session, marker_name, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a set of all DB marker sequences; set dist to 2."""
    global max_dist_genus
    setup_onebp(session, marker_name, shared_tmp_dir, debug=False, cpu=0)
    max_dist_genus = 2


def setup_dist3(session, marker_name, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a set of all DB marker sequences; set dist to 3."""
    global max_dist_genus
    setup_onebp(session, marker_name, shared_tmp_dir, debug=False, cpu=0)
    max_dist_genus = 3


def setup_dist4(session, marker_name, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a set of all DB marker sequences; set dist to 4."""
    global max_dist_genus
    setup_onebp(session, marker_name, shared_tmp_dir, debug=False, cpu=0)
    max_dist_genus = 4


def setup_dist5(session, marker_name, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a set of all DB marker sequences; set dist to 5."""
    global max_dist_genus
    setup_onebp(session, marker_name, shared_tmp_dir, debug=False, cpu=0)
    max_dist_genus = 5


def method_dist(
    fasta_file,
    session,
    marker_name,
    read_report,
    tmp_dir,
    shared_tmp_dir,
    min_abundance=0,
    debug=False,
    cpu=0,
):
    """Classify using edit distance."""
    global db_seqs
    global genus_taxid
    global max_dist_genus
    assert max_dist_genus >= 1, max_dist_genus

    count = 0
    seqs = {}
    with open(fasta_file) as handle:
        for title, seq in SimpleFastaParser(handle):
            idn = title.split(None, 1)[0]
            abundance = abundance_from_read_name(idn)
            if min_abundance and abundance < min_abundance:
                continue
            count += abundance
            seqs[idn] = seq.upper()

    # Compute all the query vs DB distances in one call
    all_dists = cdist(
        seqs.values(),
        db_seqs.keys(),
        scorer=Levenshtein.distance,
        dtype=int8,
        score_cutoff=max_dist_genus,
    )

    results = {}
    for (idn, seq), dists in zip(seqs.items(), all_dists):
        min_dist = min(dists)
        results[idn] = 0, "", f"No matches up to distance {max_dist_genus}"
        if min_dist == 0:
            assert (
                seq in db_seqs
            ), f"Expected {idn} to be in DB as min dist zero:\n{seq}"
            # If seq in db_seqs might be genus only, and
            # we'd prefer a species level match 1bp away:
            taxid, genus_species, note = perfect_match_in_db(session, marker_name, seq)
            results[idn] = taxid, genus_species, note
            if any(species_level(_) for _ in genus_species.split(";")):
                # Found 100% identical match(es) in DB at species level, done :)
                continue
            # What's next if we ignore the perfect genus-only match?
            min_dist = min(d for d in dists if d != 0)
        if min_dist == 1:
            # No species level exact matches, so do we have 1bp off match(es)?
            # What if the 1bp genus differ from any 0bp genus? Take both!
            matches = {s for (s, d) in zip(db_seqs, dists) if d <= min_dist}
            assert matches
            taxid, genus_species, _ = taxid_and_sp_lists(
                session.query(Taxonomy)
                .join(SeqSource)
                .join(MarkerDef, SeqSource.marker_definition)
                .filter(MarkerDef.name == marker_name)
                .join(MarkerSeq)
                .filter(MarkerSeq.sequence.in_(matches))
                .distinct()
            )
            assert genus_species
            results[idn] = taxid, genus_species, ""
        elif min_dist <= max_dist_genus and not results[idn][1]:
            # Thus far code was shared for onebp and the 1sXg classifiers.
            # We now come to the genus-only fall back for 1sXg if relevant:
            results[idn] = 0, "?", "Skipped"
            # Nothing within 1bp, take genus only info from Xbp away:
            matches = {s for (s, d) in zip(db_seqs, dists) if d <= min_dist}
            assert matches
            note = f"{len(matches)} matches at distance {min_dist}"
            genus = {g for s in matches for g in db_seqs[s]}
            # genus = {g for ((s, g_set), d) in zip(db_seqs.items(), dists)
            #          if d==min_dist for g in g_set}
            assert genus
            results[idn] = (
                ";".join(str(genus_taxid.get(g, 0)) for g in sorted(genus)),
                ";".join(sorted(genus)),
                note,
            )
        assert results[idn]

    tax_counts = Counter()
    for idn in seqs:
        taxid, genus_species, note = results[idn]
        # TODO - avoid calling abundance_from_read_name twice
        tax_counts[genus_species] += abundance_from_read_name(idn)
        if debug:
            read_report.write(f"{idn}\t{str(taxid)}\t{genus_species}\t{note}\n")
        else:
            read_report.write(f"{idn}\t{str(taxid)}\t{genus_species}\n")
    assert count == sum(tax_counts.values())
    return tax_counts


def setup_blast(session, marker_name, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a BLAST DB from the marker sequence DB entries."""
    view = (
        session.query(MarkerSeq)
        .join(SeqSource)
        .join(MarkerDef, SeqSource.marker_definition)
        .filter(MarkerDef.name == marker_name)
    )
    db_fasta = os.path.join(shared_tmp_dir, "blast_db.fasta")
    blast_db = os.path.join(shared_tmp_dir, "blast_db")
    count = 0
    with open(db_fasta, "w") as handle:
        for marker in view:
            md5 = marker.md5
            marker_seq = marker.sequence
            handle.write(f">{md5}\n{marker_seq}\n")
            count += 1
    sys.stderr.write(
        f"Wrote {count} unique sequences from DB to FASTA file for BLAST database.\n"
    )
    cmd = ["makeblastdb", "-dbtype", "nucl", "-in", db_fasta, "-out", blast_db]
    run(cmd, debug)


def method_blast(
    fasta_file,
    session,
    marker_name,
    read_report,
    tmp_dir,
    shared_tmp_dir,
    min_abundance=0,
    debug=False,
    cpu=0,
):
    """Classify using BLAST.

    Another simplistic classifier, run the reads through blastn
    against a BLAST database of our marker sequence database entries.
    """
    assert os.path.isdir(tmp_dir)
    assert os.path.isdir(shared_tmp_dir)

    blast_out = os.path.join(shared_tmp_dir, "blast.tsv")
    blast_db = os.path.join(shared_tmp_dir, "blast_db")
    if not (
        os.path.isfile(blast_db + ".nhr")
        and os.path.isfile(blast_db + ".nin")
        and os.path.isfile(blast_db + ".nsq")
    ):
        sys.exit(f"ERROR: Missing generated BLAST database {blast_db}.n*\n")
    if min_abundance:
        old = fasta_file
        fasta_file = os.path.join(tmp_dir, os.path.basename(old))
        if debug:
            sys.stderr.write(
                f"DEBUG: Applying minimum abundance filter to {old}, "
                f"making {fasta_file}\n"
            )
        abundance_filter_fasta(old, fasta_file, min_abundance)
        del old
    cmd = [
        "blastn",
        "-db",
        blast_db,
        "-query",
        fasta_file,
        "-perc_identity",
        "95",
        "-outfmt",
        "6",
        "-out",
        blast_out,
    ]
    if cpu:
        cmd += ["-num_threads", str(cpu)]
    run(cmd, debug)

    if not os.path.isfile(blast_out):
        sys.exit(f"ERROR: BLAST did not produce expected output file {blast_out}\n")

    # We want to report on entries without a BLAST hit,
    # and they will be missing in the BLAST output.
    # Therefore must look at the FASTA input file too.

    query_length = OrderedDict()
    with open(fasta_file) as handle:
        for title, seq in SimpleFastaParser(handle):
            idn = title.split(None, 1)[0]
            query_length[idn] = len(seq)

    # Load the top-equal BLAST results into a dict, values are lists hit MD5,
    # and the associated score in a second dict
    blast_hits = {}
    blast_score = {}
    score = None
    with open(blast_out) as handle:
        for line in handle:
            # if debug:
            #     sys.stderr.write(line)
            parts = line.rstrip("\n").split("\t")
            idn = parts[0]
            if float(parts[3]) / query_length[idn] < MIN_BLAST_COVERAGE:
                # Too short
                continue
            if idn not in blast_hits:
                blast_hits[idn] = [parts[1]]
                score = float(parts[11])
                blast_score[idn] = parts[11]  # as string
            elif score == float(parts[11]):
                # Tied hit
                blast_hits[idn].append(parts[1])

    tax_counts = Counter()
    for idn in query_length:
        abundance = abundance_from_read_name(idn)
        if min_abundance:
            assert min_abundance <= abundance, idn
        if idn in blast_hits:
            db_md5s = blast_hits[idn]
            score = blast_score[idn]
            taxid, genus_species, note = taxid_and_sp_lists(
                session.query(Taxonomy)
                .join(SeqSource)
                .join(MarkerDef, SeqSource.marker_definition)
                .filter(MarkerDef.name == marker_name)
                .join(MarkerSeq)
                .filter(MarkerSeq.md5.in_(db_md5s))
                .distinct()
            )
            note = (f"{len(db_md5s)} BLAST hits (bit score {score}). {note}").strip()
        else:
            taxid = 0
            genus_species = ""
            note = "No DB match"
        if debug:
            read_report.write(f"{idn}\t{str(taxid)}\t{genus_species}\t{note}\n")
        else:
            read_report.write(f"{idn}\t{str(taxid)}\t{genus_species}\n")
        tax_counts[genus_species] += abundance
    return tax_counts


def method_cleanup():
    """Free any memory and/or delete any files on disk.

    Currently no need to generalise this for the different classifiers, but
    could if for example we also needed to delete any files on disk.
    """
    global db_seqs, max_dist_genus
    db_seqs = None  # global variable for onbep and 1s?g classifiers
    max_dist_genus = None  # global variable for 1s?g distance classifier


method_tool_check = {
    "blast": ["makeblastdb", "blastn"],
    "identity": [],
    "onebp": [],
    "1s2g": [],
    "1s3g": [],
    "1s4g": [],
    "1s5g": [],
    "substr": [],
}

method_classify_file = {
    "blast": method_blast,
    "identity": method_identity,
    "onebp": method_dist,
    "1s2g": method_dist,
    "1s3g": method_dist,
    "1s4g": method_dist,
    "1s5g": method_dist,
    "substr": method_substr,
}

# Not all methods define a setup function:
method_setup = {
    "blast": setup_blast,
    "onebp": setup_onebp,
    "1s2g": setup_dist2,
    "1s3g": setup_dist3,
    "1s4g": setup_dist4,
    "1s5g": setup_dist5,
}


def main(
    fasta,
    session,
    marker_name,
    method,
    out_dir,
    ignore_prefixes,
    tmp_dir,
    min_abundance=0,
    debug=False,
    cpu=0,
):
    """Implement the ``thapbi_pict classify`` command.

    For use in the pipeline command, returns a filename list of the TSV
    classifier output.

    The FASTA files should have been prepared with the same or a lower minimum
    abundance - this acts as an additional filter useful if exploring the best
    threshold.
    """
    global genus_taxid
    assert isinstance(fasta, list)

    if method not in method_classify_file:
        sys.exit(
            f"ERROR: Invalid method name {method!r},"
            f" should be one of: {', '.join(sorted(method_classify_file))}\n"
        )
    classify_file_fn = method_classify_file[method]
    try:
        setup_fn = method_setup[method]
    except KeyError:
        setup_fn = None
    try:
        req_tools = method_tool_check[method]
    except KeyError:
        req_tools = []
    check_tools(req_tools, debug)

    count = session.query(Taxonomy).distinct(Taxonomy.genus, Taxonomy.species).count()
    if debug:
        sys.stderr.write(f"Taxonomy table contains {count} distinct species.\n")
    if not count:
        sys.exit("ERROR: Taxonomy table empty, cannot classify anything.\n")
    genus_taxid = {}  # reset any values from a previous DB

    if not marker_name:
        view = session.query(MarkerDef)
        if view.count() > 1:
            sys.exit("ERROR: Need -k / --marker when DB has multiple amplicon markers")
        marker = view.one_or_none()
        if not marker:
            sys.exit("ERROR: Need DB has no amplicon markers defined")
        marker_name = marker.name
        if debug:
            sys.stderr.write(
                f"DEBUG: Assuming want {marker_name} as only marker in the DB\n"
            )
        del marker, view

    # Now want to get the number of species associated with marker DB entries,
    #
    # $ sqlite3 L5-2019-01-01.sqlite "SELECT DISTINCT taxonomy.genus, taxonomy.species
    # FROM taxonomy JOIN its1_source ON taxonomy.id == its1_source.taxonomy_id
    # ORDER BY taxonomy.genus, taxonomy.species;" | wc -l
    # 143
    #
    # $ sqlite3 L5-2019-01-01.sqlite "SELECT DISTINCT taxonomy.genus, taxonomy.species
    # FROM taxonomy ORDER BY taxonomy.genus, taxonomy.species;" | wc -l
    # 253
    #
    # Note the number with an NCBI taxid could be lower...
    view = (
        session.query(Taxonomy)
        .distinct(Taxonomy.genus, Taxonomy.species)
        .join(SeqSource, SeqSource.taxonomy_id == Taxonomy.id)
        .join(MarkerDef, SeqSource.marker_definition)
        .filter(MarkerDef.name == marker_name)
    )
    db_sp_list = sorted({genus_species_name(t.genus, t.species) for t in view})
    assert "" not in db_sp_list
    if debug:
        sys.stderr.write(
            f"{marker_name} DB entries linked to {len(db_sp_list)} distinct species.\n"
        )
    if not db_sp_list:
        sys.exit(
            f"ERROR: Have no {marker_name} sequences in DB with species information.\n"
        )

    count = (
        session.query(MarkerSeq.sequence)
        .join(SeqSource)
        .join(MarkerDef, SeqSource.marker_definition)
        .filter(MarkerDef.name == marker_name)
        .distinct()
        .count()
    )
    if debug:
        sys.stderr.write(f"DB contains {count} distinct sequences for {marker_name}.\n")
    if not count:
        sys.exit(f"ERROR: No {marker_name} sequences, cannot classify anything.\n")

    fasta_files = find_requested_files(fasta, ".fasta", ignore_prefixes, debug=debug)
    if debug:
        sys.stderr.write(f"Classifying {len(fasta_files)} input FASTA files\n")

    if out_dir and out_dir != "-" and not os.path.isdir(out_dir):
        sys.stderr.write(f"Making output directory {out_dir!r}\n")
        os.mkdir(out_dir)

    if tmp_dir:
        # Up to the user to remove the files
        tmp_obj = None
        shared_tmp = tmp_dir
    else:
        tmp_obj = tempfile.TemporaryDirectory()
        shared_tmp = tmp_obj.name

    if debug:
        sys.stderr.write(f"DEBUG: Shared temp folder {shared_tmp}\n")

    classifier_output = []  # return value

    seq_count = 0
    match_count = 0
    skipped_samples = set()
    for filename in fasta_files:
        sys.stdout.flush()
        sys.stderr.flush()

        folder, stem = os.path.split(filename)
        stem = os.path.splitext(stem)[0]
        if not out_dir:
            # Use input folder
            output_name = os.path.join(folder, f"{stem}.{method}.tsv")
        elif out_dir == "-":
            output_name = None
        else:
            output_name = os.path.join(out_dir, f"{stem}.{method}.tsv")

        classifier_output.append(output_name)

        if output_name is not None and os.path.isfile(output_name):
            skipped_samples.add(output_name)
            if debug:
                sys.stderr.write(f"Skipping {output_name} as already done\n")
            # Count the number of sequences and matches
            with open(output_name) as handle:
                for line in handle:
                    if line.startswith("#"):
                        continue
                    parts = line.rstrip("\n").split("\t")
                    # MD5_abundance, taxids, species, notes
                    a = abundance_from_read_name(parts[0])
                    if min_abundance and a < min_abundance:
                        del a
                        continue
                    seq_count += a
                    if parts[2].strip():
                        match_count += a
                    del a, parts
            continue

        if setup_fn:
            # There are some files still to process, do setup now (once only)
            setup_fn(session, marker_name, shared_tmp, debug, cpu)
            setup_fn = None

        sys.stderr.write(f"Running {method} classifier on {filename}\n")
        if debug:
            sys.stderr.write(f"DEBUG: Output {output_name}\n")

        tmp = os.path.join(shared_tmp, stem)
        if not os.path.isdir(tmp):
            # If using tempfile.TemporaryDirectory() for shared_tmp
            # this will be deleted automatically, otherwise user must:
            os.mkdir(tmp)

        if debug:
            sys.stderr.write(f"DEBUG: Temp folder of {stem} is {tmp}\n")
        # Using same file names, but in tmp folder:
        tmp_pred = os.path.join(tmp, f"{stem}.{method}.tsv")
        # Run the classifier and write the sequence report:
        if output_name is None:
            pred_handle = sys.stdout
        else:
            pred_handle = open(tmp_pred, "w")

        # Could write one column per db_sp_list entry, but would be very sparse.
        if debug:
            pred_handle.write(
                f"#{marker_name}/sequence-name\ttaxid\tgenus-species\tnote\n"
            )
        else:
            pred_handle.write(f"#{marker_name}/sequence-name\ttaxid\tgenus-species\n")
        # Will get empty dict for empty file
        headers = load_fasta_header(filename)
        if not headers or headers["abundance"]:
            # There are sequences to classify
            assert os.path.isdir(tmp), tmp
            tax_counts = classify_file_fn(
                filename,
                session,
                marker_name,
                pred_handle,
                tmp,
                shared_tmp,
                min_abundance=min_abundance,
                debug=debug,
                cpu=cpu,
            )
        else:
            sys.stderr.write(
                f"WARNING: Skipping {method} classifier on {filename}"
                " as zero sequences\n"
            )
            tax_counts = Counter()
            pred_handle.write("#(no sequences to classify)\n")

        # Record the taxonomy counts
        count = sum(tax_counts.values())
        seq_count += count
        match_count += count - tax_counts.get("", 0)

        if output_name is not None:
            pred_handle.close()
            # Move our temp file into position...
            shutil.move(tmp_pred, output_name)

    method_cleanup()

    if skipped_samples:
        sys.stderr.write(
            f"Skipped {len(skipped_samples)} previously classified samples\n"
        )

    if tmp_dir:
        sys.stderr.write(
            f"WARNING: Please remove temporary files written to {tmp_dir}\n"
        )
    else:
        tmp_obj.cleanup()

    sys.stderr.write(
        f"{method} classifier assigned species/genus to {match_count}"
        f" of {seq_count} sequences from {len(fasta_files)} files\n"
    )

    sys.stdout.flush()
    sys.stderr.flush()
    return classifier_output
