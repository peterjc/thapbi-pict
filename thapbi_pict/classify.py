# Copyright 2018-2021 by Peter Cock, The James Hutton Institute.
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

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Levenshtein import distance as levenshtein
from sqlalchemy.orm import aliased
from sqlalchemy.orm import contains_eager

from .db_orm import connect_to_db
from .db_orm import ITS1
from .db_orm import SequenceSource
from .db_orm import Taxonomy
from .utils import abundance_filter_fasta
from .utils import abundance_from_read_name
from .utils import cmd_as_string
from .utils import find_requested_files
from .utils import genus_species_name
from .utils import load_fasta_header
from .utils import md5seq
from .utils import md5seq_16b
from .utils import onebp_variants
from .utils import run
from .utils import species_level
from .versions import check_tools


MIN_BLAST_COVERAGE = 0.85  # percentage of query length

fuzzy_matches = None  # global variable for onebp classifier
db_seqs = None  # global variable for 1s3g classifier
max_dist_genus = None  # global variable for 1s?g distance classifier


def md5_to_taxon(md5_list, session):
    """Return all the taxon entries linked to given marker sequences.

    Each marker sequence (matched by MD5 checksum) should be used
    in one more more source table entries, each of which will have
    a current taxonomy entry.
    """
    # This is deliberately following the query style used in dump
    # (since at some point we'll want to define database subsets
    # consistently). Still, might be able to refactor this query...
    cur_tax = aliased(Taxonomy)
    its1_seq = aliased(ITS1)
    view = (
        session.query(SequenceSource)
        .join(its1_seq, SequenceSource.its1)
        .join(cur_tax, SequenceSource.current_taxonomy)
        .options(contains_eager(SequenceSource.its1, alias=its1_seq))
        .filter(its1_seq.md5.in_(md5_list))
        .options(contains_eager(SequenceSource.current_taxonomy, alias=cur_tax))
    )
    return list({_.current_taxonomy for _ in view})  # dedulicate


def unique_or_separated(values, sep=";"):
    """Return sole element, or a string joining all elements using the separator."""
    if len(set(values)) == 1:
        return values[0]
    else:
        return sep.join(str(_) for _ in values)


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
    if not taxon_entries:
        # Unexpected - this is perhaps worth an assert statement / debug msg?
        return 0, "", "No taxonomy entries"

    genus_with_species = {t.genus for t in taxon_entries if t.species}
    if genus_with_species:
        # Want either species level predictions, or genus level without a
        # child species level prediction.
        #
        # Note ignoring the TaxID here - would need to know the parent/child
        # relationship to confirm the genus we're removing does have species
        # level children in the prediction set.
        taxon_entries = [
            t for t in taxon_entries if t.species or t.genus not in genus_with_species
        ]

    if len(taxon_entries) == 1:
        t = taxon_entries[0]
        return (
            t.ncbi_taxid,
            genus_species_name(t.genus, t.species),
            "Unique taxonomy match",
        )

    # Discard clade, and now remove duplicates, sort on genus-sepcies
    tax = sorted({(t.genus, t.species, t.ncbi_taxid) for t in taxon_entries})

    return (
        unique_or_separated([t[2] for t in tax]),
        unique_or_separated([genus_species_name(t[0], t[1]) for t in tax]),
        "",  # Not useful to report # of entries as redundant info
    )


def perfect_match_in_db(session, seq, debug=False):
    """Lookup sequence in DB, returns taxid, genus_species, note as tuple.

    If the 100% matches in the DB give multiple species, then taxid and
    genus_species will be semi-colon separated strings.
    """
    assert seq == seq.upper(), seq
    # Now, does this equal any of the ITS1 seq in our DB?
    its1 = session.query(ITS1).filter(ITS1.sequence == seq).one_or_none()
    if its1 is None:
        return 0, "", "No DB match"
    assert its1.sequence == seq
    # its1 -> one or more SequenceSource
    # each SequenceSource -> one current taxonomy
    # TODO: Refactor the query to get the DB to apply disinct?
    t = list(
        {_.current_taxonomy for _ in session.query(SequenceSource).filter_by(its1=its1)}
    )
    if not t:
        sys.exit(
            f"ERROR: perfect_match_in_db, no taxonomy for id={its1.id} md5={its1.md5}"
            f" sequence={its1.sequence}\n"
        )
    return taxid_and_sp_lists(t)


def perfect_substr_in_db(session, seq, debug=False):
    """Lookup sequence in DB, returns taxid, genus_species, note as tuple.

    If the matches containing the sequence as a substring give multiple species,
    then taxid and genus_species will be semi-colon separated strings.
    """
    assert seq == seq.upper(), seq
    # Now, does this match any of the ITS1 seq in our DB (as a substring)
    # This didn't work:
    #
    # its1 = session.query(ITS1).filter(ITS1.sequence.match(seq)).one_or_none()
    #
    # Gave:
    #
    #     sqlalchemy.exc.OperationalError: (sqlite3.OperationalError) unable to
    #     use function MATCH in the requested context
    #     [SQL: SELECT its1_sequence.id AS its1_sequence_id,
    #     its1_sequence.md5 AS its1_sequence_md5,
    #     its1_sequence.sequence AS its1_sequence_sequence
    #     FROM its1_sequence
    #     WHERE its1_sequence.sequence MATCH ?]
    t = set()
    for its1 in session.query(ITS1).filter(ITS1.sequence.like("%" + seq + "%")):
        assert its1 is not None
        assert seq in its1.sequence
        t.update(
            _.current_taxonomy
            for _ in session.query(SequenceSource).filter_by(its1=its1)
        )
    if not t:
        return 0, "", "No DB match"
    else:
        return taxid_and_sp_lists(t)


def apply_method_to_file(
    method_fn, fasta_file, session, read_report, min_abundance=0, debug=False
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
            taxid, genus_species, note = method_fn(session, seq.upper(), debug=debug)
            tax_counts[genus_species] += abundance
            read_report.write(f"{idn}\t{str(taxid)}\t{genus_species}\t{note}\n")
    assert count == sum(tax_counts.values())
    return tax_counts


def method_identity(
    fasta_file,
    session,
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
        read_report,
        min_abundance=min_abundance,
        debug=debug,
    )


def seq_method_identity(seq, session, debug=False):
    """Look for a perfect match in the database.

    Returns taxid (integer or string), genus-species (string), note (string).
    If there are multiple matches, semi-colon separated strings are returned.
    """
    its1 = session.query(ITS1).filter(ITS1.sequence == seq).one_or_none()
    if its1 is None:
        return 0, "", "No DB match"
    else:
        return perfect_match_in_db(session, seq)


def method_substr(
    fasta_file,
    session,
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
        read_report,
        min_abundance=min_abundance,
        debug=debug,
    )


def setup_onebp(session, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a dictionary of DB variants from the DB marker sequences.

    Because MD5 checksums are shorter than the typical marker sequences,
    they are used as the dictionary keys to save memory. We assume there
    are no collisions.
    """
    global fuzzy_matches
    fuzzy_matches = {}

    view = session.query(ITS1)
    count = 0
    for its1 in view:
        count += 1
        if count == 4000:
            # Warn once, this is a somewhat arbitrary threshold
            sys.stderr.write(
                "WARNING: Database contains so many unique sequences "
                "that building a cloud of 1bp variants may run out of memory...\n"
            )
        its1_seq = its1.sequence
        for variant in onebp_variants(its1_seq):
            md5_16b = md5seq_16b(variant)
            try:
                # This variant is 1bp different to multiple DB entries...
                fuzzy_matches[md5_16b].append(its1_seq)
            except KeyError:
                # Thus far this variant is only next to one DB entry:
                fuzzy_matches[md5_16b] = [its1_seq]

    sys.stderr.write(
        f"Expanded {count} ITS1 sequences from DB into cloud of"
        f" {len(fuzzy_matches)} 1bp different variants\n"
    )


def method_onebp(
    fasta_file,
    session,
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
        read_report,
        min_abundance=min_abundance,
        debug=debug,
    )


def onebp_match_in_db(session, seq, debug=False):
    """Look in database for a perfect match or with 1bp edit.

    Returns taxid (integer or string), genus-species (string), note (string).
    If there are multiple matches, semi-colon separated strings are returned.

    If there is a perfect genus only match, but a species level match one
    base pair away, takes that instead.
    """
    global fuzzy_matches
    taxid, genus_species, note = perfect_match_in_db(session, seq)
    if any(species_level(_) for _ in genus_species.split(";")):
        # Found 100% identical match(es) in DB at species level, done :)
        return taxid, genus_species, note

    # No species level exact matches, so do we have 1bp off match(es)?
    md5_16b = md5seq_16b(seq)
    if md5_16b in fuzzy_matches:
        t = set()
        # TODO - Refactor this 2-query-per-loop into one lookup?
        # Including [seq] here in order to retain any perfect genus match.
        # If there are any *different* genus matches 1bp away, they'll be
        # reported too, but that would most likely be a DB problem...
        for db_seq in [seq] + fuzzy_matches[md5_16b]:
            its1 = session.query(ITS1).filter(ITS1.sequence == db_seq).one_or_none()
            assert db_seq, f"Could not find {db_seq} ({md5seq(db_seq)}) in DB?"
            t.update(
                _.current_taxonomy
                for _ in session.query(SequenceSource).filter_by(its1=its1)
            )
        note = (
            f"{len(fuzzy_matches[md5_16b])} ITS1 matches with up to 1bp diff,"
            f" {len(t)} taxonomy entries"
        )
        if not t:
            sys.exit(
                f"ERROR: onebp: {len(fuzzy_matches[md5_16b])} matches"
                f" but no taxonomy entries for {seq}\n"
            )
        taxid, genus_species, _ = taxid_and_sp_lists(list(t))
    elif not genus_species:
        note = "No DB matches, even with 1bp diff"
    return taxid, genus_species, note


def _setup_dist(session, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a set of all DB marker sequences."""
    global db_seqs
    db_seqs = set()
    view = session.query(ITS1)
    for its1 in view:
        db_seqs.add(its1.sequence.upper())
    # Now set fuzzy_matches too...
    setup_onebp(session, shared_tmp_dir, debug, cpu)


def setup_dist2(session, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a set of all DB marker sequences; set dist to 2."""
    global max_dist_genus
    max_dist_genus = 2
    _setup_dist(session, shared_tmp_dir, debug=False, cpu=0)


def setup_dist3(session, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a set of all DB marker sequences; set dist to 3."""
    global max_dist_genus
    max_dist_genus = 3
    _setup_dist(session, shared_tmp_dir, debug=False, cpu=0)


def setup_dist4(session, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a set of all DB marker sequences; set dist to 4."""
    global max_dist_genus
    max_dist_genus = 4
    _setup_dist(session, shared_tmp_dir, debug=False, cpu=0)


def setup_dist5(session, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a set of all DB marker sequences; set dist to 5."""
    global max_dist_genus
    max_dist_genus = 5
    _setup_dist(session, shared_tmp_dir, debug=False, cpu=0)


def dist_in_db(session, seq, debug=False):
    """Species up to 1bp, genus up to 3bp away."""
    global db_seqs
    global fuzzy_matches
    assert seq and db_seqs and fuzzy_matches
    # If seq in db_seqs might be genus only, and
    # we'd prefer a species level match 1bp away:
    if (seq in db_seqs) or (md5seq_16b(seq) in fuzzy_matches):
        return onebp_match_in_db(session, seq)
    min_dist = 0
    best = {}
    # Any matches are at least 2bp away, will take genus only.
    # Fall back on brute force! But only on a minority of cases
    for db_seq in db_seqs:
        dist = levenshtein(seq, db_seq)
        if dist > max_dist_genus:
            pass
        elif dist == min_dist:
            # Best equal
            assert dist > 0
            best.add(db_seq)
        elif dist < min_dist or min_dist == 0:
            # New winner
            min_dist = dist
            best = {
                db_seq,
            }
    if not min_dist:
        return 0, "", f"No matches up to distance {max_dist_genus}"
    note = f"{len(best)} matches at distance {min_dist}"
    assert min_dist > 1  # Should have caught via fuzzy list!

    genus = set()
    for db_seq in best:
        its1 = session.query(ITS1).filter(ITS1.sequence == db_seq).one_or_none()
        assert db_seq, f"Could not find {db_seq} ({md5seq(db_seq)}) in DB?"
        genus.update(
            _.current_taxonomy.genus
            for _ in session.query(SequenceSource).filter_by(its1=its1)
        )
    assert genus
    # TODO - look up genus taxid
    return 0, ";".join(sorted(genus)), note


def method_dist(
    fasta_file,
    session,
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
        read_report,
        min_abundance=min_abundance,
        debug=debug,
    )


def setup_blast(session, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a BLAST DB from the marker sequence DB entries."""
    view = session.query(ITS1)
    db_fasta = os.path.join(shared_tmp_dir, "blast_db.fasta")
    blast_db = os.path.join(shared_tmp_dir, "blast_db")
    count = 0
    with open(db_fasta, "w") as handle:
        for its1 in view:
            md5 = its1.md5
            its1_seq = its1.sequence
            handle.write(f">{md5}\n{its1_seq}\n")
            count += 1
    sys.stderr.write(
        f"Wrote {count} unique sequences from DB to FASTA file for BLAST database.\n"
    )
    cmd = ["makeblastdb", "-dbtype", "nucl", "-in", db_fasta, "-out", blast_db]
    run(cmd, debug)


def method_blast(
    fasta_file,
    session,
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
            t = md5_to_taxon(db_md5s, session)
            if not t:
                sys.exit(f"ERROR: No taxon entry for {idn}")
            taxid, genus_species, note = taxid_and_sp_lists(t)
            note = (f"{len(db_md5s)} BLAST hits (bit score {score}). {note}").strip()
        else:
            taxid = 0
            genus_species = ""
            note = "No BLAST hit"
        read_report.write(f"{idn}\t{str(taxid)}\t{genus_species}\t{note}\n")
        tax_counts[genus_species] += abundance
    return tax_counts


def make_swarm_db_fasta(db_fasta, session):
    """Prepare a SWARM style ITS1 FASTA sequence from our DB."""
    view = session.query(ITS1)
    count = 0
    skip = 0
    unambig = set("ACGT")
    with open(db_fasta, "w") as handle:
        for its1 in view:
            md5 = its1.md5
            its1_seq = its1.sequence
            if not unambig.issuperset(its1_seq):
                skip += 1
                continue
            # Using md5_ref_abundance to avoid duplicating any exact
            # sequence matching read which would be md5_abundance
            abundance = 1  # need to look at sequence_source join...
            handle.write(f">{md5}_db_{abundance}\n{its1_seq}\n")
            count += 1
    sys.stderr.write(
        f"Wrote {count} unique sequences from DB to FASTA file for SWARM.\n"
    )
    if skip:
        sys.stderr.write(
            f"WARNING: Skipped {skip} DB sequences due to ambiguous bases.\n"
        )
    return count


def run_swarm(input_fasta, output_clusters, diff=1, debug=False, cpu=0):
    """Run swarm on a prepared FASTA file."""
    # Seems swarm only looks at one positional argument,
    # silently ignores any second or third FASTA file.
    # However, it will take stdin, e.g. swarm -d 1 <(cat *.nu)
    # This may be useful as we want to feed in the combination
    # of a database dump (with idn_abundance naming), plus the
    # prepared sequencing reads.
    cmd = ["swarm"]
    if cpu:
        cmd += ["-t", str(cpu)]  # -t, --threads
    cmd += ["-d", str(diff)]  # -d, --differences
    cmd += ["-o", output_clusters]  # -o, --output-file
    if isinstance(input_fasta, str):
        cmd += [input_fasta]
    else:
        # Swarm defaults to stdin, so don't need this:
        # cmd += ["/dev/stdin"]
        # Using grep to remove the header lines from 'FASTA' file:
        cmd = 'cat "%s" | grep -v "^#" | %s' % (
            '" "'.join(input_fasta),
            cmd_as_string(cmd),
        )
    return run(cmd, debug=debug)


def setup_swarm(session, shared_tmp_dir, debug=False, cpu=0):
    """Prepare files we'll reuse when running SWARM.

    Dumps the database to a swarm-ready FASTA file, which will
    later be used as input to swarm together with the prepared
    read sequnces.
    """
    db_fasta = os.path.join(shared_tmp_dir, "swarm_db.fasta")
    make_swarm_db_fasta(db_fasta, session)


def method_swarm_core(
    fasta_file,
    session,
    read_report,
    tmp_dir,
    shared_tmp_dir,
    min_abundance=0,
    identity=False,
    debug=False,
    cpu=0,
):
    """Classify using SWARM.

    Uses the previously generated dump of the database to a
    swarm-ready FASTA file, and non-redundant input FASTA, as
    input to swarm.

    Uses the database sequences to assign species to clusters,
    and thus input sequences within a cluster to that species.

    If identity=True (i.e. the 'swarmid' classifier), it will
    override the species for individual sequences with that of
    any 100% identical matches in the DB. In this mode, acts
    like the 'identity' classifier with the 'swarm' classifier
    as a fallback.
    """
    db_fasta = os.path.join(shared_tmp_dir, "swarm_db.fasta")
    if not os.path.isfile(db_fasta):
        sys.exit(f"ERROR: Missing generated file {db_fasta}\n")

    if min_abundance:
        old = fasta_file
        fasta_file = os.path.join(tmp_dir, os.path.basename(old))
        if debug:
            sys.stderr.write(
                "DEBUG: Applying minimum abundance filter, making {fasta_file}\n"
            )
        abundance_filter_fasta(old, fasta_file, min_abundance)
        del old

    if identity:
        seq_dict = SeqIO.index(fasta_file, "fasta")

    swarm_clusters = os.path.join(tmp_dir, "swarm_clusters.txt")
    run_swarm([fasta_file, db_fasta], swarm_clusters, diff=1, debug=debug, cpu=cpu)

    if not os.path.isfile(swarm_clusters):
        sys.exit(
            f"ERROR: Swarm did not produce expected output file {swarm_clusters}\n"
        )
    cluster_count = 0
    count = 0
    tax_counts = Counter()
    with open(swarm_clusters) as handle:
        for line in handle:
            cluster_count += 1
            idns = line.strip().split()
            # This split is safe if the sequence came though our prepare-reads
            read_idns = [_ for _ in idns if "_db_" not in _]
            db_md5s = [_.split("_db_", 1)[0] for _ in idns if "_db_" in _]
            del idns
            abundance = sum(abundance_from_read_name(_) for _ in read_idns)
            count += abundance
            if not read_idns:
                # DB only cluster, ignore
                continue
            if db_md5s:
                t = md5_to_taxon(db_md5s, session)
                taxid, genus_species, note = taxid_and_sp_lists(t)
                note = (
                    f"Cluster #{cluster_count} - {len(read_idns)} seqs"
                    f" and {len(db_md5s)} DB entries. {note}"
                ).strip()
            else:
                # Cannot classify, report
                taxid = 0
                genus_species = ""
                note = (
                    f"Cluster #{cluster_count} - {len(read_idns)} seqs"
                    " but no DB entry"
                )
            for idn in read_idns:
                if identity:
                    # Does this match any of the ITS1 seq in our DB?
                    seq = str(seq_dict[idn].seq).upper()
                    taxid2, genus_species2, _ = perfect_match_in_db(session, seq)
                    if genus_species2:
                        read_report.write(
                            "%s\t%s\t%s\t%s\n"
                            % (
                                idn,
                                str(taxid2),
                                genus_species2,
                                f"Cluster #{cluster_count} - {len(read_idns)} seqs,"
                                " but this seq itself in DB",
                            )
                        )
                        continue
                read_report.write(f"{idn}\t{str(taxid)}\t{genus_species}\t{note}\n")
            tax_counts[genus_species] += abundance
    sys.stderr.write(f"Swarm generated {cluster_count} clusters\n")
    assert count == sum(tax_counts.values())
    return tax_counts


def method_swarm(
    fasta_file,
    session,
    read_report,
    tmp_dir,
    shared_tmp_dir,
    min_abundance=0,
    debug=False,
    cpu=0,
):
    """SWARM classifier.

    Uses the previously generated dump of the database to a
    swarm-ready FASTA file, and the prepared non-redundant input
    FASTA to input to swarm.

    Uses the database sequences to assign species to clusters,
    and thus input sequences within a cluster to that species.
    """
    return method_swarm_core(
        fasta_file,
        session,
        read_report,
        tmp_dir,
        shared_tmp_dir,
        identity=False,
        min_abundance=min_abundance,
        debug=debug,
        cpu=cpu,
    )


def method_swarmid(
    fasta_file,
    session,
    read_report,
    tmp_dir,
    shared_tmp_dir,
    min_abundance=0,
    debug=False,
    cpu=0,
):
    """SWARM classifier with 100% identity special case.

    Proceeds as per the simple SWARM classifier, with each cluster
    being assigned species. However, before using the cluster species,
    checks if the sequence has a 100% identical match in the DB and
    if so, uses that in preference.

    i.e. This is like method_identity falling back on method_swarm.
    """
    return method_swarm_core(
        fasta_file,
        session,
        read_report,
        tmp_dir,
        shared_tmp_dir,
        identity=True,
        min_abundance=min_abundance,
        debug=debug,
        cpu=cpu,
    )


method_tool_check = {
    "blast": ["makeblastdb", "blastn"],
    "identity": [],
    "onebp": [],
    "1s2g": [],
    "1s3g": [],
    "1s4g": [],
    "1s5g": [],
    "substr": [],
    "swarm": ["swarm"],
    "swarmid": ["swarm"],
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
    "swarm": method_swarm,
    "swarmid": method_swarmid,
}

# Not all methods define a setup function:
method_setup = {
    "blast": setup_blast,
    "onebp": setup_onebp,
    "1s2g": setup_dist2,
    "1s3g": setup_dist3,
    "1s4g": setup_dist4,
    "1s5g": setup_dist5,
    "swarm": setup_swarm,
    "swarmid": setup_swarm,  # can share the setup with swarm
}


def main(
    fasta,
    db_url,
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

    # Connect to the DB,
    Session = connect_to_db(db_url, echo=False)  # echo=debug is too distracting now
    session = Session()

    count = session.query(Taxonomy).distinct(Taxonomy.genus, Taxonomy.species).count()
    if debug:
        sys.stderr.write(f"Taxonomy table contains {count} distinct species.\n")
    if not count:
        sys.exit("ERROR: Taxonomy table empty, cannot classify anything.\n")

    # Now want to get the number of species associated with ITS1 DB entries,
    #
    # $ sqlite3 L5-2019-01-01.sqlite "SELECT DISTINCT taxonomy.genus, taxonomy.species
    # FROM taxonomy JOIN its1_source ON taxonomy.id == its1_source.current_taxonomy_id
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
        .join(SequenceSource, SequenceSource.current_taxonomy_id == Taxonomy.id)
    )
    db_sp_list = sorted({genus_species_name(t.genus, t.species) for t in view})
    assert "" not in db_sp_list
    if debug:
        sys.stderr.write(
            f"Marker entries in DB linked to {len(db_sp_list)} distrinct species.\n"
        )
    if not db_sp_list:
        sys.exit("ERROR: Have no marker sequences in DB with species information.\n")

    count = session.query(ITS1).count()
    if debug:
        sys.stderr.write(f"ITS1 table contains {count} distinct sequences.\n")
    if not count:
        sys.exit("ERROR: ITS1 table empty, cannot classify anything.\n")

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
            setup_fn(session, shared_tmp, debug, cpu)
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
        pred_handle.write("#sequence-name\ttaxid\tgenus-species\tnote\n")
        # Will get empty dict for empty file
        headers = load_fasta_header(filename)
        if not headers or headers["abundance"]:
            # There are sequences to classify
            assert os.path.isdir(tmp), tmp
            tax_counts = classify_file_fn(
                filename,
                session,
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
