# ------------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# ------------------------------------------------------------------------------

SHELL=/bin/bash -o pipefail

test:
	py.test -v --cov=aeneas --doctest-modules tests/*.py

install:
	pip install .

devenv:
	pip install pytest pytest-cov pep8

style:
	pep8 aeneas/*.py tests/*.py

loc:
	cloc --exclude-list-file=<(echo aeneas/_version.py) aeneas/*.py
	cloc tests/test_*.py

clean:
	find . -type d -name __pycache__ -maxdepth 2 -exec rm -r {} \;
