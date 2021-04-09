#!/usr/bin/env python
"""Run console snippets in RST documentation, check their output.

Example usage::

    $ cd examples/woody_hosts/
    $ ../scripts/rst_doc_test.py ../../docs/examples/woody_hosts/*.rst

Takes one or more RST filenames. These are scanned for ``.. code:: console``
entries which are then parsed to pull out commands (starting ``$ `` with slash
line continuation supported) and expected terminal output. The special output
line ``...`` can be used to indicate omitted output. The commands are run in
the current directory.

STRONG WARNING: Do not use this on untrusted input files!
"""
# TODO - Configurable skipped commands prefix/suffix list?
import os
import subprocess
import sys
import tempfile

if len(sys.argv) == 1 or not os.path.isfile(sys.argv[1]):
    sys.exit("ERROR: Requires one or more RST filesnames as input.")


def scan_rst(filename):
    """Look for console blocks in an RST file."""
    with open(filename) as handle:
        lines = list(handle)
        while lines:
            line = lines.pop(0)
            if line == ".. code:: console\n":
                line = lines.pop(0)
                assert not line.strip()
                block = []
                while lines:
                    line = lines.pop(0)
                    if not line.strip("\n").strip():
                        break
                    if not block and line.startswith("     "):
                        # Over indented
                        sys.exit(
                            "ERROR: Console entry not four space indented"
                            f" in {filename}:\n{line}"
                        )
                        break
                    elif not line.startswith("    "):
                        # Under indented
                        sys.exit(
                            "ERROR: Console entry not four space indented"
                            f" in {filename}:\n{line}"
                        )
                        break
                    elif line.strip() == "<SEE TABLE BELOW>":
                        # Magic happens, expect a blank line ending console
                        # section, a paragraph of text, and then an RST Table
                        line = lines.pop(0)
                        if line != "\n":
                            sys.exit(
                                "ERROR: Expected blank line after <SEE TABLE BELOW>"
                                f" in {filename}:\n{line!r}"
                            )
                        line = lines.pop(0)
                        if line[0] == " ":
                            sys.exit(
                                "ERROR: Expected paragraph after <SEE TABLE BELOW>"
                                f" in {filename}:\n{line!r}"
                            )
                        while line.strip():
                            line = lines.pop(0)
                        assert not line.strip()
                        line = lines.pop(0)
                        # Now should be a table!
                        if line.startswith("    ="):
                            sys.exit(
                                "ERROR: Don't indent the table after <SEE TABLE BELOW>"
                                f" in  {filename}"
                            )
                        if not line.startswith("="):
                            sys.exit(
                                "ERROR: Expected table after <SEE TABLE BELOW>"
                                f" in {filename}:\n{line!r}"
                            )
                        while line.strip():
                            block.append(line)
                            line = lines.pop(0)
                        lines.insert(0, line)  # Put the blank line back
                    else:
                        block.append(line[4:])
                if not block[0].startswith("$ "):
                    sys.exit(
                        "ERROR: Console entry does not start four space dollar space"
                        f" in {filename}:\n{block[0]}"
                    )
                assert "\n" not in block, block
                yield block
            elif ".. code:: console" in line:
                sys.exit(f"ERROR: {filename} has this line:\n{repr(line)}")


def parse_block(block):
    """Look for command and output pairs in an RST console block."""
    block = list(block)  # copy
    assert block[0].startswith("$ ")
    while block:
        out = ""
        line = block.pop(0)
        assert line.startswith("$ "), line
        cmd = line[2:].lstrip().rstrip()
        while cmd.endswith("\\"):
            line = block.pop(0)
            assert not line.startswith("$ "), line
            cmd = cmd[:-1].strip() + " " + line.strip()
        while block and not block[0].startswith("$ "):
            out += block.pop(0)
        yield cmd, out


def run_cmd(cmd):
    """Run a shell command, return stdout and sterr as strings."""
    child = subprocess.run(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    return child.stdout, child.stderr


def tsv_align(text, min_pad=2):
    """Align columns of TSV output."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as handle:
        handle.write(text)
        filename = handle.name
    output = subprocess.getoutput(f"xsv table -p {min_pad} -d '\t' {filename}")
    os.remove(filename)
    # Remove any trailing whitespace for RST readiness
    output = "\n".join(line.rstrip() for line in output.split("\n"))
    return output.rstrip("\n") + "\n"


def tsv_to_rst(text):
    """Turn TSV output into a simple RST table."""
    first_line, rest = text.split("\n", 1)
    headers = first_line.split("\t")
    text = tsv_align(text, min_pad=1)
    lines = text.split("\n")
    width = max(len(_) for _ in lines)
    first_line, rest = text.split("\n", 1)

    cuts = []
    index = len(headers[0])
    for col_name in headers[1:]:
        index = first_line.index(col_name, index)
        assert first_line[index] != " " and index > 0 and first_line[index - 1] == " "
        cuts.append(index - 1)
    del index

    header = "".join(" " if i in cuts else "=" for i in range(width)) + "\n"
    return header + first_line + "\n" + header + rest + header


def fasta_wrap(text):
    """Line wrap FASTA output at 80 characters."""
    new = []
    lines = text.split("\n")
    while lines:
        line = lines.pop(0)
        if line.startswith(">"):
            while len(line) > 80 and " " in line:
                cut = line
                while len(cut) > 80 and " " in cut:
                    cut, rest = cut.rsplit(" ", 1)
                new.append(cut)
                line = line[len(cut) + 1 :].lstrip()
                del cut
            new.append(line)
            line = lines.pop(0)
            while len(line) > 80:
                new.append(line[:80])
                line = line[80:]
        new.append(line)
    return "\n".join(new)


assert (
    fasta_wrap(">Silly\n" + "ATCG" * 25) == ">Silly\n" + "ATCG" * 20 + "\n" + "ATCG" * 5
), fasta_wrap(">Silly\n" + "ATCG" * 100)

assert (
    fasta_wrap(">Silly" + " blah" * 15 + "\n" + "ATCG")
    == ">Silly" + " blah" * 14 + "\nblah\nATCG"
), fasta_wrap(">Silly" + " blah" * 15 + "\n" + "ATCG")
assert (
    fasta_wrap(">Silly!" + " blah" * 15 + "\n" + "ATCG")
    == ">Silly!" + " blah" * 14 + "\nblah\nATCG"
), fasta_wrap(">Silly!" + " blah" * 15 + "\n" + "ATCG")
assert (
    fasta_wrap(">Silly!!" + " blah" * 15 + "\n" + "ATCG")
    == ">Silly!!" + " blah" * 14 + "\nblah\nATCG"
), fasta_wrap(">Silly!!" + " blah" * 15 + "\n" + "ATCG")
assert (
    fasta_wrap(">Silly!!?" + " blah" * 15 + "\n" + "ATCG")
    == ">Silly!!?" + " blah" * 14 + "\nblah\nATCG"
), fasta_wrap(">Silly!!?" + " blah" * 15 + "\n" + "ATCG")


cur_dir = os.path.abspath(os.curdir)
tmp_dir = "/tmp/"
errors = 0
for filename in sys.argv[1:]:
    if not filename.endswith(".rst"):
        sys.stderr.write(f"WARNING: Ignoring non-RST file {filename}\n")
        continue
    os.chdir(cur_dir)

    print(filename)

    for block in scan_rst(filename):
        for cmd, old_out in parse_block(block):
            print("$ " + cmd)
            if cmd.startswith(
                (
                    "thapbi_pict ... ",
                    "md5sum -c ",
                    "curl ",
                    "pip ",
                    "conda ",
                    "pre-commit ",
                    "git ",
                    "sudo ",
                    "apt-get ",
                )
            ):
                continue
            if cmd.endswith("# Are you sure?"):
                continue
            if cmd.startswith("cd "):
                os.chdir(cmd[3:])
                continue

            new_out, err_out = run_cmd(cmd)
            if new_out and not new_out.endswith("\n"):
                new_out += "\n"

            # Transform the output
            if new_out.startswith(">") or "\n>" in new_out:
                new_out = fasta_wrap(new_out)
            elif "\t" in new_out and old_out.startswith("="):
                new_out = tsv_to_rst(new_out)
            elif "\t" in new_out:
                new_out = tsv_align(new_out)

            # Compare the output
            if old_out == "...\n":
                pass
            elif old_out == new_out:
                if err_out:
                    # Warning?
                    print(err_out)
                pass
            elif old_out == new_out + err_out or old_out == err_out + new_out:
                pass
            elif old_out.startswith("...\n") and err_out.endswith(old_out[3:]):
                pass
            elif (
                old_out.startswith("...\n")
                and old_out.endswith("...\n")
                and (old_out[3:-4] in new_out or old_out[3:-4] in err_out)
            ):
                pass
            elif cmd.startswith("ls ") and set(old_out.split("\n")).issubset(
                new_out.split("\n")
            ):
                # Ignore extra files (e.g. from running other commands)
                pass
            else:
                print("---- Expected:")
                print(old_out)
                print("---- New stdout:")
                print(new_out)
                print("---- New stderr:")
                print(err_out)
                print("---- End")
                errors += 1
    print("File done")
    print()
if errors:
    sys.stderr.write(f"ERROR: {errors} examples failed\n")
sys.exit(errors)
