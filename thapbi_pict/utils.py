"""Helper functions for THAPB-PICT code."""

import subprocess
import sys


def cmd_as_string(cmd):
    """Express a list command as a suitably quoted string.

    Intended for using in debugging or error messages.
    """
    if isinstance(cmd, list):
        # Quote any entries with spaces
        return ' '.join('"%s"' % _ if ' ' in _ else _ for _ in cmd)
    else:
        return cmd


def run(cmd, debug=False):
    """Run a command via subprocess, abort if fails."""
    if debug:
        sys.stderr.write("Calling command: %s\n" % cmd_as_string(cmd))
    try:
        # On Python 3.7 onwards, could use capture_output=True
        # rather than stdout=PIPE and stderr=PIPE
        if isinstance(cmd, list):
            return subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                check=True)
        else:
            return subprocess.run(
                cmd, shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                check=True)
    except subprocess.CalledProcessError as e:
        if debug:
            # Used universal_newlines=True above so that this just works
            # (equivalent to text=True in Python 3.7 onwards):
            sys.stdout.write(e.stdout)
            sys.stderr.write(e.stderr)
        sys.exit("This command failed with return code %i:\n%s\n"
                 % (e.returncode, cmd_as_string(cmd)))
