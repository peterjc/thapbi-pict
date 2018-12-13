"""Classifying prepared ITS1 reads using an ITS1 database.

This implementes the ``thapbi_pict classify-reads ...`` command.
"""

import os
import shutil
import sys
import tempfile

from collections import Counter

from .db_orm import connect_to_db
from .db_orm import ITS1, SequenceSource, Taxonomy
from .hmm import filter_for_ITS1


def taxonomy_consensus(taxon_entries):
    """Return LCA summary of the taxonomy objects from DB.

    Expects a de-duplicated list of Taxonomy table entries.

    Returns a tuple of strings, starting with genus, species and clade.

    Currently also reports a text note for debugging.
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

    note = "Consensus from %i taxonomy entries" % len(taxon_entries)

    # e.g. Clades of "", "8a" --> "8a"
    # but any conflict -> ""
    c_list = list(set(_.clade for _ in taxon_entries))
    if "" in c_list:
        c_list.remove("")
    clade = c_list[0] if len(c_list) == 1 else ""
    if not clade:
        note += " (clades: %s)" % ",".join(sorted(c_list))

    s_list = list(set(_.species for _ in taxon_entries))
    if "" in s_list:
        s_list.remove("")
    species = s_list[0] if len(s_list) == 1 else ""
    if not species:
        note += " (species: %s)" % ", ".join(sorted(s_list))

    return (genus, species, clade, note)


def method_identity(fasta_file, session, read_report, tax_report,
                    tmp_dir, debug=False, cpu=0):
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
    tax_counts = Counter()
    for title, seq, its1_seqs in filter_for_ITS1(fasta_file):
        count += 1
        idn = title.split(None, 1)[0]
        genus = species = clade = note = ""
        if not its1_seqs:
            note = "No ITS1 HMM matches"
        else:
            for its1_seq in its1_seqs:
                assert its1_seq == its1_seq.upper()
            # Now, does this match any of the ITS1 seq in our DB?
            its1 = session.query(ITS1).filter(
                ITS1.sequence.in_(its1_seqs)).one_or_none()
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
                tax_counts[(genus, species, clade)] += 1
        read_report.write(
            "%s\t%s\t%s\t%s\t%s\n" % (idn, genus, species, clade, note))
    for (genus, species, clade), tax_count in sorted(tax_counts.items()):
        tax_report.write(
            "%s\t%s\t%s\t%i\n" % (genus, species, clade, tax_count))
    return count, matched


methods = {
    "identity": method_identity,
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


def main(fasta, db_url, method, out_dir, debug=False, cpu=0):
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
    for filename in fasta_files:
        sys.stderr.write(
            "Running %s classifer on %s\n" % (method, filename))

        folder, stem = os.path.split(filename)
        stem = os.path.splitext(stem)[0]
        if out_dir and out_dir != "-":
            folder = out_dir
        read_name = os.path.join(
            folder, "%s.%s-reads.tsv" % (stem, method))
        tax_name = os.path.join(
            folder, "%s.%s-tax.tsv" % (stem, method))

        if os.path.isfile(read_name) and os.path.isfile(tax_name):
            sys.stderr.write(
                "WARNING: Skipping %s and %s as already exist\n"
                % (read_name, tax_name))
            # TODO: Count the number of reads and matches?
            continue

        if debug:
            sys.stderr.write(
                "DEBUG: Outputs %s and %s\n" % (read_name, tax_name))

        # Context manager should remove the temp dir:
        with tempfile.TemporaryDirectory() as tmp:
            if debug:
                sys.stderr.write(
                    "DEBUG: Temp folder of %s is %s\n" % (stem, tmp))
            # Using same file names, but in tmp folder:
            tmp_reads = os.path.join(
                tmp, "%s.%s-reads.tsv" % (stem, method))
            tmp_tax = os.path.join(
                tmp, "%s.%s-tax.tsv" % (stem, method))
            # Run the classifier...
            with open(tmp_reads, "w") as reads_handle:
                with open(tmp_tax, "w") as tax_handle:
                    r_count, m_count = method_fn(
                        filename, session,
                        reads_handle, tax_handle,
                        tmp, debug)
            read_count += r_count
            match_count += m_count
            # Move our temp files into position...
            shutil.move(tmp_reads, read_name)
            shutil.move(tmp_tax, tax_name)

    sys.stderr.write(
        "%s classifier assigned species to %i of %i reads from %i files\n"
        % (method, match_count, read_count, len(fasta_files)))

    return 0
