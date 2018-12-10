"""Classifying prepared ITS1 reads using an ITS1 database.

This implementes the ``thapbi_pict classify-reads ...`` command.
"""

import os
import sys

from .db_orm import connect_to_db
from .db_orm import ITS1, SequenceSource, Taxonomy
from .hmm import filter_for_ITS1


def taxonomy_consensus(taxon_entries):
    """Return LCA summary of the taxonomy objects from DB.

    Expects a de-duplicated list of Taxonomy table entries.

    Returns a tuple of strings, starting with genus, species and clade.
    """
    if not taxon_entries:
        return "", "", "", "No taxonomy entries"
    if len(taxon_entries) == 1:
        t = taxon_entries[0]
        return t.genus, t.species, t.clade, "Unique taxonomy match"

    tmp = list(set(_.genus for _ in taxon_entries))
    if "" in tmp:
        tmp.remove("")
    genus = tmp[0] if len(tmp) == 1 else ""

    if not genus:
        return "", "", "", "Conflicting genera"

    # e.g. Clades of "", "8a" --> "8a"
    # but any conflict -> ""
    tmp = list(set(_.clade for _ in taxon_entries))
    if "" in tmp:
        tmp.remove("")
    clade = tmp[0] if len(tmp) == 1 else ""

    tmp = list(set(_.species for _ in taxon_entries))
    if "" in tmp:
        tmp.remove("")
    species = tmp[0] if len(tmp) == 1 else ""

    return genus, species, clade, "Consensus"


def method_identity(fasta_file, session, read_report, debug=False):
    """Classify using perfect identity.

    This is a deliberately simple approach, in part for testing
    purposes. It uses HMMER3 to find any ITS1 match, and then
    looks for a perfect identical entry in the database.
    """
    # The FASTA file should have long sequences which might
    # contain a known ITS1 sequence as a substring. If the
    # search were inverted, the SQL LIKE command could likely
    # be used (e.g. via SQLalchemey's contains operator).
    #
    # Plan B is brute force - we can run hmmscan to find any
    # ITS1 matchs, and then look for 100% equality in the DB.
    count = 0
    matched = 0
    for title, seq, its1_seq in filter_for_ITS1(fasta_file):
        count += 1
        idn = title.split(None, 1)[0]
        genus = species = clade = note = ""
        if not its1_seq:
            # TODO - Handle multiple HMM matches here
            note = "No single ITS1 HMM match"
        else:
            assert its1_seq == its1_seq.upper()
            # Now, does this match any of the ITS1 seq in our DB?
            its1 = session.query(ITS1).filter_by(
                sequence=its1_seq).one_or_none()
            if its1 is None:
                note = "No ITS1 database match"
            else:
                # its1 -> one or more SequenceSource
                # each SequenceSource -> one current taxonomy
                # TODO: Refactor the query to get the DB to apply disinct?
                genus, species, clade, note = taxonomy_consensus(
                    list(set(_.current_taxonomy for _ in session.query(
                        SequenceSource).filter_by(its1=its1))))
                matched += 1
        read_report.write(
            "%s\t%s\t%s\t%s\t%s\n" % (idn, genus, species, clade, note))
    return count, matched


methods = {
    "identity": method_identity
}


def find_fasta_files(filenames_or_folders, ext=".fasta", debug=False):
    """Interpret a list of filenames and/or foldernames."""
    answer = []
    for x in filenames_or_folders:
        if os.path.isdir(x):
            if debug:
                sys.stderr.write("Walking directory %r\n" % x)
            for f in os.listdir(x):
                if f.endswith(ext):
                    # Check not a directory?
                    answer.append(os.path.join(x, f))
        elif os.path.isfile(x):
            answer.append(x)
        else:
            sys.exit("ERROR: %r is not a file or a directory\n" % x)
    # Warn if there were duplicates?
    return sorted(set(answer))


def main(fasta, db_url, method, read_report, debug=False):
    """Implement the thapbi_pict classify-reads command."""
    assert isinstance(fasta, list)

    if method not in methods:
        sys.exit(
            "ERROR: Invalid method name %r, should be one of: %s\n"
            % (method, ", ".join(sorted(methods))))
    method_fn = methods[method]

    # Connect to the DB,
    Session = connect_to_db(db_url, echo=debug)
    session = Session()

    count = session.query(Taxonomy).distinct(
        Taxonomy.genus, Taxonomy.species).count()
    if debug:
        sys.stderr.write(
            "Taxonomy table contains %i distinct species.\n" % count)
    if not count:
        sys.exit(
            "ERROR: Taxonomy table empty, cannot classify anything.\n")

    count = session.query(ITS1).count()
    if debug:
        sys.stderr.write(
            "ITS1 table contains %i distinct sequences.\n" % count)
    if not count:
        sys.exit(
            "ERROR: ITS1 table empty, cannot classify anything.\n")

    fasta_files = find_fasta_files(fasta, debug=debug)
    if debug:
        sys.stderr.write(
            "Classifying %i input FASTA files\n" % len(fasta_files))

    read_count = 0
    match_count = 0
    if read_report == "-":
        read_handle = sys.stdout
    else:
        read_handle = open(read_report, "w")
    for filename in fasta_files:
        if debug:
            sys.stderr.write("%s classifer on %s\n" % (method, filename))
        r_count, m_count = method_fn(filename, session, read_handle, debug)
        read_count += r_count
        match_count += m_count

    sys.stderr.write(
        "%s classifier assigned species to %i of %i reads from %i files\n"
        % (method, match_count, read_count, len(fasta_files)))
    if read_report != "-":
        read_handle.close()

    return 0
