"""Extract HMM matches from sequences (using hmmer3).

This is part of the THAPBI *Phytophthora* ITS1 Classifier Tool (PICT).
We have a Hidden Markov Model (HMM) based on an alignment of manually
curated *Phytophthora* ITS1 sequences, and we will use this to find
ITS1 regions within sequences of potential interest.

The first use case is to trim down ITS1 sequences during import to our
database (for example, sequences submitted to GenBank will often be
longer than the region we are focused on). Second, we can apply this
during our classification pipeline - both to sort sequencing reads into
ITS1 matching or not, but also to trim to a standard region for simpler
comparison.

The code uses Biopython's SearchIO module to parse the output from the
command line tool hmmscan from the hmmer3 suite.
"""

import os
import shutil
import sys
import tempfile
import time
import warnings

from ..utils import run

from Bio import BiopythonExperimentalWarning
from Bio import SeqIO

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    # This experimenal code (beta) warning was dropped in Biopython 1.73
    from Bio import SearchIO


def run_and_parse_hmmscan(
    hmm_file,
    fasta_input_file,
    hmmscan="hmmscan",
    bitscore_threshold=None,
    debug=False,
    cpu=0,
):
    """Run hmmscan and parse the output with Bio.SearchIO.

    This will run the command line tool ``hmmscan`` (assumed to be on the
    ``$PATH`` but an explicit path can be given), writing the output to a
    temporary file which will be deleted after parsing it.

    This will yield 2-tuples of Biopython objects via the FASTA parser in
    Bio.SeqIO, and the HMMER3 text parser in Bio.SearchIO.
    """
    if not os.path.isfile(fasta_input_file):
        raise ValueError("Missing FASTA input %r" % fasta_input_file)
    if not os.path.isfile(hmm_file):
        raise ValueError("Missing HMM input %r" % hmm_file)

    # TODO - Once this is working, refactor to parse stdout
    tmp_dir = tempfile.mkdtemp()
    hmm_out = os.path.join(tmp_dir, "hmmscan.txt")
    cmd = [hmmscan, "--noali"]
    if cpu:
        cmd += ["--cpu", str(cpu)]
    if bitscore_threshold is not None:
        cmd += ["-T", bitscore_threshold, "--domT", bitscore_threshold]
    cmd += ["-o", hmm_out, hmm_file, fasta_input_file]
    # cmd = "'%s' --noali -o '%s' '%s' '%s'" % (hmmscan, hmm_out, hmm_file,
    #                                           fasta_input_file)
    run(cmd, debug=debug, attempts=3)

    if debug:
        sys.stderr.write("DEBUG: hmmscan finished, about to parse the output\n")

    for record, result in zip(
        SeqIO.parse(fasta_input_file, "fasta"), SearchIO.parse(hmm_out, "hmmer3-text")
    ):
        if record.id != result.id:
            raise ValueError("FASTA %s vs HMMER %s" % (record.id, result.id))
        yield record, result

    os.remove(hmm_out)
    shutil.rmtree(tmp_dir)


def hmm_cache(hmm_file, cache_dir, debug=False):
    """Cache HMM files for calling hmmscan."""
    old_dir, stem = os.path.split(hmm_file)

    new = os.path.join(cache_dir, stem)
    if os.path.isfile(new):
        # Already cached
        return new

    if debug:
        sys.stderr.write("DEBUG: Caching HMM at %s\n" % new)

    # Copy the extra files first:
    for f in os.listdir(old_dir):
        if f.startswith(stem + ".h"):
            shutil.copyfile(os.path.join(old_dir, f), os.path.join(cache_dir, f))

    # Do this last in case interupted (as later only check this exists):
    shutil.copyfile(hmm_file, new)

    return new


def load_result_cache(cached_results_file):
    """Load simple TSV file of cached ITS1 trim results as a dict.

    Will silently ignore duplicate entries, taking the later copy.
    e.g. Rival processes appended the same entry.
    """
    answer = dict()
    claim_lock(cached_results_file)
    for line in open(cached_results_file):
        seq, start, end = line.rstrip("\n").split("\t")
        answer[seq] = (int(start), int(end))
    release_lock(cached_results_file)
    return answer


def claim_lock(cached_results_file, timeout=240, debug=False):
    """Wait until can creat log file, or abort."""
    wait = 0
    pause = 1
    while wait < timeout and os.path.isfile(cached_results_file + ".lock"):
        if debug:
            sys.stderr.write("Waiting for lock to write to %s\n" % cached_results_file)
        time.sleep(pause)
        wait += pause
        pause = max(60, pause * 2)
    if os.path.isfile(cached_results_file + ".lock"):
        sys.exit("Could not acquire lock to update %s" % cached_results_file)
    with open(cached_results_file + ".lock", "w") as handle:
        handle.write("In use\n")


def release_lock(cached_results_file):
    """Delete lock file."""
    try:
        os.remove(cached_results_file + ".lock")
    except OSError:
        # Probably a race condition... cache could be incomplete now
        sys.stderr.write("ERROR: Problem releasing %s.lock\n" % cached_results_file)


def save_result_cache(cached_results_file, value_dict, timeout=60, debug=False):
    """Save ITS1 trim results as a dict."""
    claim_lock(cached_results_file, timeout=60, debug=debug)
    if os.path.isfile(cached_results_file):
        # Reload first in case altered in meantime
        for line in open(cached_results_file):
            seq, start, end = line.rstrip("\n").split("\t")
            value_dict[seq] = (int(start), int(end))
    if debug:
        sys.stderr.write(
            "DEBUG: Updating %s to have %i entries\n"
            % (cached_results_file, len(value_dict))
        )
    with open(cached_results_file, "w") as handle:
        for seq, (start, end) in value_dict.items():
            handle.write("%s\t%i\t%i\n" % (seq, start, end))
    release_lock(cached_results_file)


def append_result_cache(cached_results_file, new_values, timeout=60, debug=False):
    """Append new ITS1 trim results to the cache file.

    This could lead to rival processess appending the same value,
    however we discard duplicates on loading the cache.
    """
    claim_lock(cached_results_file, timeout=60, debug=debug)
    if debug:
        sys.stderr.write(
            "DEBUG: Updating %s with %i new entries\n"
            % (cached_results_file, len(new_values))
        )
    with open(cached_results_file, "a") as handle:
        for seq, start, end in new_values:
            handle.write("%s\t%i\t%i\n" % (seq, start, end))
    release_lock(cached_results_file)


def update_result_cache(
    cached_results_file,
    value_dict,
    missing_seqs,
    tmp_dir,
    model_cache_dir,
    debug=False,
    cpu=0,
):
    """Call hmmscan, update dict in memory, and write it to disk.

    Has simple file locking.
    """
    if not missing_seqs:
        return
    tmp_fasta = os.path.join(tmp_dir, "input_to_hmmscan.fasta")
    with open(tmp_fasta, "w") as handle:
        for title, seq in missing_seqs:
            handle.write(">%s\n%s\n" % (title, seq))
    if debug:
        sys.stderr.write(
            "DEBUG: Created temp file %s for calling hmmscan\n" % tmp_fasta
        )
    new_values = []
    for _title, seq, trimmed in filter_for_ITS1(
        tmp_fasta, model_cache_dir, bitscore_threshold="6", debug=debug, cpu=cpu
    ):
        start = seq.index(trimmed)
        end = start + len(trimmed)
        value_dict[seq] = (start, end)
        new_values.append((seq, start, end))
    append_result_cache(cached_results_file, new_values, debug=debug)


def cached_filter_for_ITS1(
    title_and_seqs, cached_results_file, model_cache_dir, tmp_dir, debug=False, cpu=0
):
    """Return tuples of (input title, input sequence, ITS1 sequence) in arbitrary order.

    Takes advantage of the specified cache (TSV file mapping full sequences to the
    ITS1 region), and for novel sequences runs hmmscan and updates the cache.

    If the ITS1 region is not found, then an empty string is used.
    """
    if os.path.isfile(cached_results_file):
        result_cache = load_result_cache(cached_results_file)
    else:
        result_cache = dict()
        if debug:
            sys.stderr.write("DEBUG: No cache at %s\n" % cached_results_file)

    cached_seqs = set()  # only for debug msg
    missing_seqs = set()
    for title, seq in title_and_seqs:
        try:
            start, end = result_cache[seq]
            yield title, seq, seq[start:end]
            cached_seqs.add(seq)
        except KeyError:
            missing_seqs.add((title, seq))

    if debug:
        sys.stderr.write(
            "DEBUG: Using hmmscan to trim %i unique sequences (%i cached, %i new)\n"
            % (
                len(cached_seqs) + len(missing_seqs),
                len(cached_seqs),
                len(missing_seqs),
            )
        )

    if missing_seqs:
        update_result_cache(
            cached_results_file,
            result_cache,
            missing_seqs,
            tmp_dir,
            model_cache_dir,
            debug=debug,
            cpu=cpu,
        )
        if debug:
            sys.stderr.write(
                "DEBUG: Now have %i cached ITS1 trim results from hmmscan\n"
                % len(result_cache)
            )
        for title, seq in missing_seqs:
            start, end = result_cache[seq]
            yield title, seq, seq[start:end]


def filter_for_ITS1(input_fasta, cache_dir, bitscore_threshold="6", debug=False, cpu=0):
    """Search for the expected single ITS1 sequence within FASTA entries.

    The arbitrary low bitscore_threshold default is based on ensuring
    all the sequences in our legacy FASTA files pass without false
    positive additional weak hits complicating things.
    """
    hmm = os.path.join(os.path.split(__file__)[0], "phytophthora_its1.hmm")
    if cache_dir:
        hmm = hmm_cache(hmm, cache_dir, debug=debug)

    for record, result in run_and_parse_hmmscan(
        hmm, input_fasta, bitscore_threshold=bitscore_threshold, debug=debug, cpu=cpu
    ):
        title = record.description
        seq = str(record.seq).upper()
        if not result:
            yield title, seq, ""
            continue

        assert len(result) == 1
        hit = result[0]
        assert hit.id == "phytophthora_its1"

        if len(hit) == 0:
            yield title, seq, ""
        elif len(hit) == 1:
            hsp = hit[0]
            yield title, seq, seq[hsp.query_start : hsp.query_end]
        else:
            # Merge them - does not seem useful to insist non-overlapping
            start = min(_.query_start for _ in hit)
            end = max(_.query_end for _ in hit)
            if len(hit) > 2:
                sys.stderr.write(
                    "WARNING: Taking span %i-%i from %i HMM hits for:\n"
                    ">%s\n%s\n" % (start, end, len(hit), title, seq)
                )
            yield title, seq, seq[start:end]
