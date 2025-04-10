# Copyright 2018-2024 by Peter Cock, The James Hutton Institute.
# Revisions copyright 2024 by Peter Cock, University of Strathclyde.
# All rights reserved.
# This file is part of the THAPBI Phytophthora ITS1 Classifier Tool (PICT),
# and is released under the "MIT License Agreement". Please see the LICENSE
# file that should have been included as part of this package.
"""setuptools based setup script for THAPBI PICT.

This uses setuptools which is now the standard python mechanism for installing
packages (and is used internally by the tool pip). If you have downloaded and
uncompressed the THAPBI PICT source code, or fetched it from git, for the
simplest installation try one of these commands::

    pip install .

Or::

    python setup.py install

Typically however, since we have released this software on PyPI, you don't
need to download it first. Instead just run::

    pip install thapbi_pict

Once installed, you should be able to run the tool using:

   thapbi_pict

Or start it as a module::

   python -m thapbi_pict

Either should find the installed copy of the Python code.
"""

import sys

try:
    from setuptools import find_packages
    from setuptools import setup
except ImportError:
    sys.exit(
        "We need the Python library setuptools to be installed. "
        "Try running: python -m ensurepip"
    )


# Make sure we have the right Python version.
if sys.version_info[:2] < (3, 10):  # noqa: UP036
    sys.exit(
        "THAPBI PICT requires Python 3.10 or later. "
        "Python %d.%d detected.\n" % sys.version_info[:2]
    )

# We define the version number in thapbi_pict/__init__.py
# Here we can't use "import thapbi_pict" then "thapbi_pict.__version__"
# as that would tell us the version already installed (if any).
__version__ = "Undefined"
with open("thapbi_pict/__init__.py") as handle:
    for line in handle:
        if line.startswith("__version__ = "):
            exec(line.strip())
            break

# Load our rsStructuredText file README.rst as the long description.
#
# Without declaring an encoding, decoding a problematic character in the file
# may fail on Python 3 depending on the user's locale. By explicitly checking
# it is ASCII (could use latin1 or UTF8 if needed later), if any invalid
# character does appear in our README, this will fail for everyone.
with open("README.rst", "rb") as handle:
    readme_rst = handle.read().decode("ascii")

setup(
    name="thapbi_pict",
    version=__version__,
    author="Peter Cock",
    author_email="peter.cock@hutton.ac.uk",
    url="https://github.com/peterjc/thapbi-pict",  # For now at least
    download_url="https://github.com/peterjc/thapbi-pict",
    description="THAPBI Phytophthora ITS1 Classifier Tool (PICT).",
    project_urls={
        "Documentation": "https://thapbi-pict.readthedocs.io/",
        "Source": "https://github.com/peterjc/thapbi-pict/",
        "Tracker": "https://github.com/peterjc/thapbi-pict/issues",
    },
    long_description=readme_rst,
    long_description_content_type="text/x-rst",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Software Development :: Libraries :: Python Modules",
    ],
    license="MIT License",
    license_files=["LICENSE.rst"],
    python_requires=">=3.10",
    entry_points={"console_scripts": ["thapbi_pict = thapbi_pict.__main__:main"]},
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "biopython >=1.84",
        "cutadapt >=4.9",
        "xopen >=2.0.2",
        "matplotlib >=3.9.2",
        "networkx >=3.4.2",
        "pydot >=3.0.1",
        "rapidfuzz >=3.10.1",
        "sqlalchemy >=2.0.36",
        "xlsxwriter >=3.2.0",
    ],
)
