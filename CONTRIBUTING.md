# THAPBI *Phytophthora* ITS1 Classifier Tool (PICT) - Development Notes

## Dependencies

We will use the [Conda](https://conda.io/) packaging system, including the
[BioConda](https://bioconda.github.io/) channel for bioinformatics packages,
and [PyPI](https://pypi.python.org/) packages to define a reproducible
environment for running this software.

## Python style conventions

The Python code follows [PEP8](https://www.python.org/dev/peps/pep-0008/)
and [PEP257 docstring](https://www.python.org/dev/peps/pep-0257/) style,
guided by the [Zen of Python](https://www.python.org/dev/peps/pep-0020/).

Practically, coding style is enforced with the command line tools
[black](https://github.com/ambv/black) (which can automatically edit
your code) via [flake8](http://flake8.pycqa.org/) (which in addition to
its own style checking, has a range of plugins - including a plugin to
call black from flake8).

You can install these tools using:

```console
pip install black flake8 flake8-black flake8-blind-except flake8-docstrings flake8-rst-docstrings restructuredtext-lint flake8-bugbear
```

You can run the checks using:

```console
$ flake8 .
```

You can ask black to edit your files with:

```console
$ black .
```

We recommend enabling the flake8 git pre-commit hook as follows:

```console
$ flake8 --install-hook git
$ git config --bool flake8.strict true
```

If your editor can be configured to run flake8 and/or black automatically,
even better.
