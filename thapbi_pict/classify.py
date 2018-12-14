"""Classifying prepared ITS1 reads using an ITS1 database.

This implementes the ``thapbi_pict classify-reads ...`` command.
"""

import os
import shutil
import sys
import tempfile

from collections import Counter

from sqlalchemy.orm import aliased, contains_eager

from .db_orm import connect_to_db
from .db_orm import ITS1, SequenceSource, Taxonomy
from .hmm import filter_for_ITS1
from .utils import cmd_as_string, run


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
    view = session.query(SequenceSource).join(
        its1_seq, SequenceSource.its1).join(
        cur_tax, SequenceSource.current_taxonomy).options(
        contains_eager(SequenceSource.its1, alias=its1_seq)).filter(
        its1_seq.md5.in_(md5_list)).options(
        contains_eager(SequenceSource.current_taxonomy, alias=cur_tax))
    return list(set(_.current_taxonomy for _ in view))  # depulicate


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


def method_identity(fasta_file, session, read_report,
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
        tax_counts[(genus, species, clade)] += 1
        read_report.write(
            "%s\t%s\t%s\t%s\t%s\n" % (idn, genus, species, clade, note))
    assert count == sum(tax_counts.values())
    return tax_counts


def make_swarm_db_fasta(db_fasta, session):
    """Prepare a SWARM style ITS1 FASTA sequence from our DB."""
    view = session.query(ITS1)
    count = 0
    skip = 0
    unambig = set('ACGT')
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
            handle.write(">%s_db_%i\n%s\n" % (md5, abundance, its1_seq))
            count += 1
    sys.stderr.write(
        "Wrote %i unique sequences from DB to FASTA file for SWARM.\n"
        % count)
    if skip:
        sys.stderr.write(
            "WARNING: Skipped %i DB sequences due to ambiguous bases.\n"
            % skip)
    return count


def run_swarm(input_fasta, output_clusters,
              diff=1,
              debug=False, cpu=0):
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


def method_swarm(fasta_file, session, read_report,
                 tmp_dir, debug=False, cpu=0):
    """Classify using SWARM.

    Dumps the database to a swarm-ready FASTA file, and gives
    this and the prepared non-redundant input FASTA to swarm.
    Uses the database sequences to assign species to clusters,
    and thus reads within a cluster to that species.
    """
    db_fasta = os.path.join(tmp_dir, "swarm_db.fasta")
    swarm_clusters = os.path.join(tmp_dir, "swarm_clusters.txt")
    make_swarm_db_fasta(db_fasta, session)
    run_swarm([fasta_file, db_fasta], swarm_clusters, diff=1,
              debug=debug, cpu=cpu)

    if not os.path.isfile(swarm_clusters):
        sys.exit(
            "Swarm did not produce expected output file %s\n"
            % swarm_clusters)
    cluster_count = 0
    count = 0
    tax_counts = Counter()
    with open(swarm_clusters) as handle:
        for line in handle:
            cluster_count += 1
            idns = line.strip().split()
            # This split is safe if the reads came though our prepare-reads
            read_idns = [_ for _ in idns if "_db_" not in _]
            db_md5s = [_.split("_db_", 1)[0] for _ in idns if "_db_" in _]
            del idns
            read_count = len(read_idns)
            count += read_count
            if not read_idns:
                # DB only cluster, ignore
                continue
            if db_md5s:
                genus, species, clade, note = taxonomy_consensus(
                    md5_to_taxon(db_md5s, session))
                note = ("Cluster of %i reads and %i DB entries. %s"
                        % (len(read_idns), len(db_md5s), note))
            else:
                # Cannot classify, report
                genus = species = clade = ""
                note = "Cluster of %i reads but no DB entry" % len(read_idns)
            for idn in read_idns:
                read_report.write(
                    "%s\t%s\t%s\t%s\t%s\n"
                    % (idn, genus, species, clade, note))
            tax_counts[(genus, species, clade)] += read_count
    sys.stderr.write("Swarm generated %i clusters\n" % cluster_count)
    assert count == sum(tax_counts.values())
    return tax_counts


methods = {
    "identity": method_identity,
    "swarm": method_swarm,
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
            # Run the classifier and write the read report:
            with open(tmp_reads, "w") as reads_handle:
                tax_counts = method_fn(
                    filename, session,
                    reads_handle,
                    tmp, debug)
            # Record the taxonomy counts
            count = sum(tax_counts.values())
            read_count += count
            match_count += count - tax_counts.get(("", "", ""))
            with open(tmp_tax, "w") as tax_handle:
                for (genus, species, clade), tax_count in sorted(
                        tax_counts.items()):
                    tax_handle.write(
                        "%s\t%s\t%s\t%i\n"
                        % (genus, species, clade, tax_count))

            # Move our temp files into position...
            shutil.move(tmp_reads, read_name)
            shutil.move(tmp_tax, tax_name)

    sys.stderr.write(
        "%s classifier assigned species to %i of %i reads from %i files\n"
        % (method, match_count, read_count, len(fasta_files)))

    return 0
