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
from rapidfuzz.string_metric import levenshtein

from .db_orm import MarkerDef
from .db_orm import MarkerSeq
from .db_orm import SeqSource
from .db_orm import Taxonomy
from .utils import abundance_filter_fasta
from .utils import abundance_from_read_name
from .utils import find_requested_files
from .utils import genus_species_name
from .utils import load_fasta_header
from .utils import md5seq_16b
from .utils import onebp_variants
from .utils import run
from .utils import species_level
from .versions import check_tools

MIN_BLAST_COVERAGE = 0.85  # percentage of query length

use_fuzzy_only = False  # global variable for onebp classifier
fuzzy_matches = None  # global variable for onebp classifier
db_seqs = None  # global variable for onebp and 1s?g distance classifiers
max_dist_genus = None  # global variable for 1s?g distance classifiers


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
assert (
    consoliate_and_sort_taxonomy(
        [
            ("Genie", "alpha", 0),
            ("Genie", "", 0),
            ("Genie", "", 100),
            ("Genie", "alpha", 101),
        ]
    )
    == [("Genie", "alpha", 101)]
)


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
    global use_fuzzy_only
    global fuzzy_matches

    db_seqs = {
        marker.sequence.upper()
        for marker in session.query(MarkerSeq.sequence)
        .join(SeqSource)
        .join(MarkerDef, SeqSource.marker_definition)
        .filter(MarkerDef.name == marker_name)
        .distinct()
    }

    fuzzy_matches = {}
    if len(db_seqs) <= 1000:
        # Build a cloud, likely to be faster until being run on a trivial
        # sample set. Will be much faster on an tiny ad-hoc DB.
        use_fuzzy_only = True
        for marker_seq in db_seqs:
            for variant in onebp_variants(marker_seq):
                md5_16b = md5seq_16b(variant)
                try:
                    # This variant is 1bp different to multiple DB entries...
                    fuzzy_matches[md5_16b].append(marker_seq)
                except KeyError:
                    # Thus far this variant is only next to one DB entry:
                    fuzzy_matches[md5_16b] = [marker_seq]
        sys.stderr.write(
            f"Expanded {len(db_seqs)} marker sequences from DB into cloud of"
            f" {len(fuzzy_matches)} 1bp different variants\n"
        )
    else:
        unambiguous = set("ACGT")
        count = 0
        for seq in db_seqs:
            bad = [_ for _ in seq if _ not in unambiguous]
            if len(bad) == 1:
                # This has a chance to be matched with the onebp classifier
                count += 1
                for letter in unambiguous:
                    md5_16b = md5seq_16b(seq.replace(bad[0], letter))
                    try:
                        # There must be another similar bad sequence
                        # e.g. ACGTwACGT and ACGTnACGT
                        fuzzy_matches[md5_16b].append(seq)
                    except KeyError:
                        fuzzy_matches[md5_16b] = [seq]
        if debug and fuzzy_matches:
            sys.stderr.write(
                f"DEBUG: Cloud of {len(fuzzy_matches)} from {count} DB entries"
                " with one ambiguity\n"
            )


def method_onebp(
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
    """Classify using identity or 1bp difference.

    This is a deliberately simple approach, based on the perfect
    identity classifier. It compares each sequence to a dictionary
    of all the database entries and their 1bp variants.
    """
    return apply_method_to_file(
        onebp_match_in_db,
        fasta_file,
        session,
        marker_name,
        read_report,
        min_abundance=min_abundance,
        debug=debug,
    )


def onebp_match_in_db(session, marker_name, seq, debug=False):
    """Look in database for a perfect match or with 1bp edit.

    Returns taxid (integer or string), genus-species (string), note (string).
    If there are multiple matches, semi-colon separated strings are returned.

    If there is a perfect genus only match, but a species level match one
    base pair away, takes that instead.
    """
    global fuzzy_matches
    global db_seqs
    taxid, genus_species, note = perfect_match_in_db(session, marker_name, seq)
    if any(species_level(_) for _ in genus_species.split(";")):
        # Found 100% identical match(es) in DB at species level, done :)
        return taxid, genus_species, note

    # No species level exact matches, so do we have 1bp off match(es)?
    assert db_seqs
    variants = set(fuzzy_matches.get(md5seq_16b(seq), [])).union([seq])
    if not use_fuzzy_only:
        # Checking variants of the query sequence, in addition to the query
        # sequence itself and any 1bp variants of DB entries with a single
        # ambiguous base:
        variants.update(onebp_variants(seq))
    taxid, genus_species, _ = taxid_and_sp_lists(
        session.query(Taxonomy)
        .join(SeqSource)
        .join(MarkerDef, SeqSource.marker_definition)
        .filter(MarkerDef.name == marker_name)
        .join(MarkerSeq)
        .filter(MarkerSeq.sequence.in_(variants))
        .distinct()
    )
    if not genus_species:
        taxid = 0
        note = "No DB matches, even with 1bp diff"
    return taxid, genus_species, note


def setup_dist2(session, marker_name, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a set of all DB marker sequences; set dist to 2."""
    global max_dist_genus
    max_dist_genus = 2
    setup_onebp(session, marker_name, shared_tmp_dir, debug=False, cpu=0)


def setup_dist3(session, marker_name, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a set of all DB marker sequences; set dist to 3."""
    global max_dist_genus
    max_dist_genus = 3
    setup_onebp(session, marker_name, shared_tmp_dir, debug=False, cpu=0)


def setup_dist4(session, marker_name, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a set of all DB marker sequences; set dist to 4."""
    global max_dist_genus
    max_dist_genus = 4
    setup_onebp(session, marker_name, shared_tmp_dir, debug=False, cpu=0)


def setup_dist5(session, marker_name, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a set of all DB marker sequences; set dist to 5."""
    global max_dist_genus
    max_dist_genus = 5
    setup_onebp(session, marker_name, shared_tmp_dir, debug=False, cpu=0)


def dist_in_db(session, marker_name, seq, debug=False):
    """Species up to 1bp, genus up to given distance away."""
    global db_seqs
    assert seq and db_seqs
    # If seq in db_seqs might be genus only, and
    # we'd prefer a species level match 1bp away:
    taxid, genus_species, note = onebp_match_in_db(session, marker_name, seq)
    if genus_species:
        return taxid, genus_species, note
    assert seq not in db_seqs

    # Else any matches from 2 bp up to X bp away, and will take genus only:
    min_dist = 0
    best = set()
    # Any matches are at least 2bp away, will take genus only.
    # Fall back on brute force! But only on a minority of cases
    for db_seq in db_seqs:
        dist = levenshtein(seq, db_seq)
        if dist > max_dist_genus:
            pass
        elif dist == min_dist:
            # Best equal
            best.add(db_seq)
        elif dist < min_dist or min_dist == 0:
            # New winner
            min_dist = dist
            best = {db_seq}
    if not best:
        assert min_dist == 0, min_dist
        return 0, "", f"No matches up to distance {max_dist_genus}"
    note = f"{len(best)} matches at distance {min_dist}"
    assert min_dist > 1  # Should have caught via onebp_match_in_db!

    # Nothing within 1bp, take genus only info from Xbp away:
    genus = {
        _.genus
        for _ in session.query(Taxonomy.genus)
        .join(SeqSource)
        .join(MarkerDef, SeqSource.marker_definition)
        .filter(MarkerDef.name == marker_name)
        .join(MarkerSeq)
        .filter(MarkerSeq.sequence.in_(best))
        .distinct()
    }
    assert genus
    # TODO - look up genus taxid
    return 0, ";".join(sorted(genus)), note


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
    return apply_method_to_file(
        dist_in_db,
        fasta_file,
        session,
        marker_name,
        read_report,
        min_abundance=min_abundance,
        debug=debug,
    )


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
    global fuzzy_matches, db_seqs, max_dist_genus
    fuzzy_matches = None  # global variable for onebp classifier
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
    "onebp": method_onebp,
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
