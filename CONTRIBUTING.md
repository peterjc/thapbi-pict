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

## Release process

For a release, start from a clean git checkout (to reduce the chance of
bundling any stray local files despite a cautious ``MANIFEST.in``).

```bash
rm -rf thapbi_pict/ITS1_DB.sqlite
sqlite3 thapbi_pict/ITS1_DB.sqlite < database/ITS1_DB.sql
chmod a-w thapbi_pict/ITS1_DB.sqlite
python setup.py sdist --formats=gztar
python setup.py bdist_wheel
twine upload dist/thapbi_pict-X.Y.Z.tar.gz dist/thapbi_pict-X.Y.Z-py3-none-any.whl
git tag vX.Y.Z
git push origin master --tags
```
