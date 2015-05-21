#!/usr/bin/env bash
pep8 aeneas/*.py
nosetests -v --with-coverage aeneas/*.py
