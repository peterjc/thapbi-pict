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

This is enforced using command line tool [flake8](http://flake8.pycqa.org/)
for which we recommend enabling the git pre-commit hook as follows:

```console
$ flake8 --install-hook git
$ git config --bool flake8.strict true
```
