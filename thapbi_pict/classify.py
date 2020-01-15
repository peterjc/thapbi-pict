# Copyright 2018-2019 by Peter Cock, The James Hutton Institute.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.

"""Classifying prepared ITS1 sequences using an ITS1 database.

This implements the ``thapbi_pict classify ...`` command.
"""

import os
import shutil
import sys
import tempfile

from collections import Counter, OrderedDict

from sqlalchemy.orm import aliased, contains_eager

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

from .db_orm import connect_to_db
from .db_orm import ITS1, SequenceSource, Taxonomy
from .utils import abundance_from_read_name
from .utils import cmd_as_string, run
from .utils import find_requested_files
from .utils import genus_species_name
from .utils import md5seq
from .utils import onebp_variants
from .utils import species_level
from .versions import check_tools


MIN_BLAST_COVERAGE = 0.85  # percentage of query length

fuzzy_matches = None  # global variable for onebp classifier


def md5_to_taxon(md5_list, session):
    """Return all the taxon entries linked to given ITS1 sequences.

    Each ITS1 sequence (matched by MD5 checksum) should be used
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
    # Now, does this match any of the ITS1 seq in our DB?
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


def apply_method_to_file(method_fn, fasta_file, session, read_report, debug=False):
    """Call given method on each sequence in the FASTA file."""
    count = 0
    tax_counts = Counter()

    with open(fasta_file) as handle:
        for title, seq in SimpleFastaParser(handle):
            idn = title.split(None, 1)[0]
            abundance = abundance_from_read_name(idn)
            count += abundance
            taxid, genus_species, note = method_fn(session, seq.upper(), debug=debug)
            tax_counts[genus_species] += abundance
            read_report.write(f"{idn}\t{str(taxid)}\t{genus_species}\t{note}\n")
    assert count == sum(tax_counts.values())
    return tax_counts


def method_identity(
    fasta_file, session, read_report, tmp_dir, shared_tmp_dir, debug=False, cpu=0
):
    """Classify using perfect identity.

    This is a deliberately simple approach, in part for testing
    purposes. It looks for a perfect identical entry in the database.
    """
    return apply_method_to_file(
        perfect_match_in_db, fasta_file, session, read_report, debug=debug
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


def setup_onebp(session, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a dictionary of DB variants from the DB marker sequences.

    Because MD5 checksums are shorter than the typical ITS1 markers,
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
            md5 = md5seq(variant)
            try:
                # This variant is 1bp different to multiple DB entries...
                fuzzy_matches[md5].append(its1_seq)
            except KeyError:
                # Thus far this variant is only next to one DB entry:
                fuzzy_matches[md5] = [its1_seq]

    sys.stderr.write(
        f"Expanded {count:d} ITS1 sequences from DB into cloud of"
        f" {len(fuzzy_matches):d} 1bp different variants\n"
    )


def method_onebp(
    fasta_file, session, read_report, tmp_dir, shared_tmp_dir, debug=False, cpu=0
):
    """Classify using identity or 1bp difference.

    This is a deliberately simple approach, based on the perfect
    identity classifier. It compares each sequence to a dictionary
    of all the database entries and their 1bp variants.
    """
    return apply_method_to_file(
        onebp_match_in_db, fasta_file, session, read_report, debug=debug
    )


def onebp_match_in_db(session, seq, debug=False):
    """Look in database for a perfect match or with 1bp edit.

    Returns taxid (integer or string), genus-species (string), note (string).
    If there are multiple matches, semi-colon separated strings are returned.
    """
    global fuzzy_matches
    taxid, genus_species, note = perfect_match_in_db(session, seq)
    if any(species_level(_) for _ in genus_species.split(";")):
        # Found 100% identical match(es) in DB at species level, done :)
        return taxid, genus_species, note

    # No species level exact matches, so do we have 1bp off match(es)?
    md5 = md5seq(seq)
    if md5 in fuzzy_matches:
        t = []
        # TODO - Refactor this 2-query-per-loop into one lookup?
        # Including [seq] here in order to retain any perfect genus match.
        # If there are any *different* genus matches 1bp away, they'll be
        # reported too, but that would most likely be a DB problem...
        for db_seq in [seq] + fuzzy_matches[md5]:
            its1 = session.query(ITS1).filter(ITS1.sequence == db_seq).one_or_none()
            assert db_seq, f"Could not find {db_seq} ({md5seq(db_seq)}) in DB?"
            t.extend(
                _.current_taxonomy
                for _ in session.query(SequenceSource).filter_by(its1=its1)
            )
        t = list(set(t))
        note = (
            f"{len(fuzzy_matches[md5]):d} ITS1 matches with up to 1bp diff,"
            f" {len(t):d} taxonomy entries"
        )
        if not t:
            sys.exit(
                f"ERROR: onebp: {len(fuzzy_matches[md5]):d} matches"
                f" but no taxonomy entries for {seq}\n"
            )
        taxid, genus_species, _ = taxid_and_sp_lists(t)
    elif not genus_species:
        note = "No DB matches, even with 1bp diff"
    return taxid, genus_species, note


def setup_blast(session, shared_tmp_dir, debug=False, cpu=0):
    """Prepare a BLAST DB from the ITS1 DB entries."""
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
        f"Wrote {count:d} unique sequences from DB to FASTA file for BLAST database.\n"
    )
    cmd = ["makeblastdb", "-dbtype", "nucl", "-in", db_fasta, "-out", blast_db]
    run(cmd, debug)


def method_blast(
    fasta_file, session, read_report, tmp_dir, shared_tmp_dir, debug=False, cpu=0
):
    """Classify using BLAST.

    Another simplistic classifier, run the ITS1 reads through blastn
    against a BLAST database of our ITS1 database entries.
    """
    blast_out = os.path.join(shared_tmp_dir, "blast.tsv")
    blast_db = os.path.join(shared_tmp_dir, "blast_db")
    if not (
        os.path.isfile(blast_db + ".nhr")
        and os.path.isfile(blast_db + ".nin")
        and os.path.isfile(blast_db + ".nsq")
    ):
        sys.exit(f"ERROR: Missing generated BLAST database {blast_db}.n*\n")
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
        if idn in blast_hits:
            db_md5s = blast_hits[idn]
            score = blast_score[idn]
            t = md5_to_taxon(db_md5s, session)
            if not t:
                sys.exit(f"ERROR: No taxon entry for {idn}")
            taxid, genus_species, note = taxid_and_sp_lists(t)
            note = (f"{len(db_md5s):d} BLAST hits (bit score {score}). {note}").strip()
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
            handle.write(f">{md5}_db_{abundance:d}\n{its1_seq}\n")
            count += 1
    sys.stderr.write(
        f"Wrote {count:d} unique sequences from DB to FASTA file for SWARM.\n"
    )
    if skip:
        sys.stderr.write(
            f"WARNING: Skipped {skip:d} DB sequences due to ambiguous bases.\n"
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
        cmd = 'cat "%s" | %s' % ('" "'.join(input_fasta), cmd_as_string(cmd))
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
                    f"Cluster #{cluster_count:d} - {len(read_idns):d} seqs"
                    f" and {len(db_md5s):d} DB entries. {note}"
                ).strip()
            else:
                # Cannot classify, report
                taxid = 0
                genus_species = ""
                note = (
                    f"Cluster #{cluster_count:d} - {len(read_idns):d} seqs"
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
    sys.stderr.write(f"Swarm generated {cluster_count:d} clusters\n")
    assert count == sum(tax_counts.values())
    return tax_counts


def method_swarm(
    fasta_file, session, read_report, tmp_dir, shared_tmp_dir, debug=False, cpu=0
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
        debug=debug,
        cpu=cpu,
    )


def method_swarmid(
    fasta_file, session, read_report, tmp_dir, shared_tmp_dir, debug=False, cpu=0
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
        debug=debug,
        cpu=cpu,
    )


method_tool_check = {
    "blast": ["makeblastdb", "blastn"],
    "identity": [],
    "onebp": [],
    "swarm": ["swarm"],
    "swarmid": ["swarm"],
}

method_classify_file = {
    "blast": method_blast,
    "identity": method_identity,
    "onebp": method_onebp,
    "swarm": method_swarm,
    "swarmid": method_swarmid,
}

method_setup = {
    "blast": setup_blast,
    # "identify": setup_identify,
    "onebp": setup_onebp,
    "swarm": setup_swarm,
    "swarmid": setup_swarm,  # can share the setup
}


def main(fasta, db_url, method, out_dir, ignore_prefixes, tmp_dir, debug=False, cpu=0):
    """Implement the ``thapbi_pict classify`` command.

    For use in the pipeline command, returns a filename list of the TSV
    classifier output.
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
        sys.stderr.write(f"Taxonomy table contains {count:d} distinct species.\n")
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
            f"ITS1 entries in DB linked to {len(db_sp_list)} distrinct species.\n"
        )
    if not db_sp_list:
        sys.exit("ERROR: Have no ITS1 entries in DB with species information.\n")

    count = session.query(ITS1).count()
    if debug:
        sys.stderr.write(f"ITS1 table contains {count:d} distinct sequences.\n")
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
    if setup_fn:
        setup_fn(session, shared_tmp, debug, cpu)

    classifier_output = []  # return value

    seq_count = 0
    match_count = 0
    for filename in fasta_files:
        sys.stderr.write(f"Running {method} classifier on {filename}\n")
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
            sys.stderr.write(f"WARNING: Skipping {output_name} as already exists\n")
            # TODO: Count the number of sequences and matches?
            continue

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
        pred_handle.write(
            "#sequence-name\ttaxid\tgenus-species:%s\tnote\n" % ";".join(db_sp_list)
        )
        if os.path.getsize(filename):
            # There are sequences to classify
            tax_counts = classify_file_fn(
                filename, session, pred_handle, tmp, shared_tmp, debug, cpu
            )
        else:
            # No sequences, no taxonomy assignments
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

    if tmp_dir:
        sys.stderr.write(
            f"WARNING: Please remove temporary files written to {tmp_dir}\n"
        )
    else:
        tmp_obj.cleanup()

    sys.stderr.write(
        f"{method} classifier assigned species/genus to {match_count:d}"
        f" of {seq_count:d} sequences from {len(fasta_files):d} files\n"
    )

    sys.stdout.flush()
    sys.stderr.flush()
    return classifier_output
