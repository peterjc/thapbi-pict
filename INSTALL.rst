Installation
============

We recommend installing this tool on Linux or macOS using the `Conda
<https://conda.io/>`__ packaging system, via the `BioConda
<https://bioconda.github.io/>`__ channel, which will handle
*all* the dependencies:

.. code:: console

   $ conda install thapbi_pict

The availability of all the command line dependencies likely rules out
using Windows.

Alternatively, since the software is on the `Python Package Index (PyPI)
<https://pypi.python.org/>`__, the following command will install it along
with its Python dependencies:

.. code:: console

   $ pip install thapbi_pict

However, in this case you will still need to install various external command
line tools like ``hmmer``, and others which are only used for some classifiers
(like ``blast`` and ``swarm``). If you have BioConda setup, use the following:

.. code:: console

   $ conda install --file requirements-ext.txt

If you are not using Conda,  then on a typical Linux system most of the tools
required will be available via the default distribution packages, although not
always under the same package name.

On Debian (with the efforts of DebianMed), or Ubuntu Linux, try:

.. code:: console

   $ sudo apt-get install ncbi-blast+ cutadapt hmmer swarm trimmomatic

If you want to install the very latest unreleased code, you must download the
source code from the `repository on GitHub
<https://github.com/peterjc/thapbi-pict>`_ - see the ``CONTRIBUTING.rst`` file
for more details.

Once installed, you should be able to run the tool using:

.. code:: console

   $ thapbi_pict

This should automatically find the installed copy of the Python code. Use
``thapbi_pict -v`` to report the version, or ``thapbi_pict -h`` for help.
