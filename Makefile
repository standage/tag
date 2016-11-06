# ------------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# ------------------------------------------------------------------------------

test:
	py.test -v --cov=aeneas --doctest-modules aeneas/*.py

install:
	pip install .

devenv:
	pip install pytest pytest-cov pep8

style:
	pep8 aeneas/*.py
