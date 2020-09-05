"""Configuration file for the Sphinx documentation builder."""
import os
import sys

import thapbi_pict

# This file only contains a selection of the most common options. For a full
# list see the documentation:
# http://www.sphinx-doc.org/en/master/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
#
# Rather than messing with the path, just requiring thapbi_pict be installed.


# -- Project information -----------------------------------------------------

project = "THAPBI PICT"
copyright = "2020, Peter Cock, James Hutton Institute;"
author = "Peter Cock"

# The full version, including alpha/beta/rc tags
release = thapbi_pict.__version__
# Should we shorten this to major.minor only?:
version = thapbi_pict.__version__

# Sphinx vs ReadTheDocs conflict on default contents.rst vs index.rst
master_doc = "index"

# -- General configuration ---------------------------------------------------

# autodoc options require Sphinx 1.8.0b1 or later:
needs_sphinx = "1.8"

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
#
# We are using SVG images which work fine in the HTML output, but need
# 'sphinx.ext.imgconverter' for it to work in the PDF output on RTD.
extensions = ["sphinx.ext.imgconverter", "sphinx.ext.autodoc", "sphinx.ext.autosummary"]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for autodoc --------------------------------------------------

# This requires Sphinx 1.8.0b1 or later:
autodoc_default_values = {
    "members": None,
    "undoc-members": None,
    "special-members": None,
    "show-inheritance": None,
    "member-order": "bysource",
    "exclude-members": "__dict__,__weakref__,__module__",
}


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# Sphinx default was html_theme = "alabaster"
html_theme = "sphinx_rtd_theme"

# Sphinx Read The Docs theme settings, see
# https://sphinx-rtd-theme.readthedocs.io/en/latest/configuring.html
html_theme_options = {"prev_next_buttons_location": "both"}

# Based on:
# https://github.com/readthedocs/sphinx_rtd_theme/issues/231#issuecomment-126447493
html_context = {
    "display_github": True,  # Add 'Edit on Github' link instead of 'View page source'
    "github_user": "peterjc",
    "github_repo": "thapbi-pict",
    "github_version": "master",
    "conf_py_path": "/docs/",
    # "source_suffix": source_suffix,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# -- Options for numpydoc -------------------------------------------------

numpydoc_class_members_toctree = False

# -- Magic to run sphinx-apidoc automatically -----------------------------

# See https://github.com/rtfd/readthedocs.org/issues/1139
# on which this is based.


def insert_github_link(filename):
    """Insert file specific :github_url: metadata for theme breadcrumbs."""
    assert "/" not in filename and filename.endswith(".rst")
    with open("api/" + filename) as handle:
        text = handle.read()
    if ":github_url:" in text:
        return

    python = filename[:-4].replace(".", "/") + "/__init__.py"
    if not os.path.isfile(os.path.join("../", python)):
        python = filename[:-4].replace(".", "/") + ".py"
    if not os.path.isfile(os.path.join("../", python)):
        sys.stderr.write(
            "WARNING: Could not map %s to a Python file, e.g. %s\n" % (filename, python)
        )
        return

    text = ":github_url: https://github.com/%s/%s/blob/%s/%s\n\n%s" % (
        html_context["github_user"],
        html_context["github_repo"],
        html_context["github_version"],
        python,
        text,
    )
    with open("api/" + filename, "w") as handle:
        handle.write(text)


def run_apidoc(_):
    """Call sphinx-apidoc on Bio and BioSQL modules."""
    from sphinx.ext.apidoc import main as apidoc_main

    apidoc_main(["-M", "-e", "-F", "-o", "api/", "../thapbi_pict"])
    os.remove("api/thapbi_pict.rst")  # replaced with index.rst

    for f in os.listdir("api/"):
        if f.endswith(".rst"):
            insert_github_link(f)


def setup(app):
    """Over-ride Sphinx setup to trigger sphinx-apidoc."""
    app.connect("builder-inited", run_apidoc)
