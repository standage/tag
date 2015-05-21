#!/usr/bin/env bash
set -eo pipefail

pep8 aeneas/*.py
nosetests -v --with-coverage aeneas/*.py
