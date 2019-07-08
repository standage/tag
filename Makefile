# ------------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# ------------------------------------------------------------------------------

SHELL=/bin/bash -o pipefail

test:
	pytest --cov=tag --doctest-modules $(shell ls tag/*.py | grep -v __init__) tag/tests/test_*.py

doc:
	cd docs && make html

install:
	pip install .

devenv:
	pip install 'pytest>=3.6,<5.0' pytest-cov pycodestyle sphinx

style:
	pycodestyle tag/*.py tag/tests/*.py tag/cli/*.py

loc:
	cloc --exclude-list-file=<(echo tag/_version.py) tag/*.py
	cloc tests/test_*.py

clean:
	find . -type d -name '*__pycache__*' -maxdepth 2 -exec rm -rf {} \;
	rm -rf tag.egg-info/ build/ dist/
