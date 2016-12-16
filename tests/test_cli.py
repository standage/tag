#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import pytest
import tag
import sys


def test_cli_args():
    oldstdout = sys.stdout

    parser = tag.cli.parser()
    sys.stdout = StringIO()
    with pytest.raises(SystemExit) as se:
        args = parser.parse_args(['-v'])
    assert tag.__version__ in sys.stdout.getvalue()

    parser = tag.cli.parser()
    sys.stdout = StringIO()
    with pytest.raises(SystemExit) as se:
        args = parser.parse_args(['gff3', '-h'])
    assert 'input file in GFF3 format' in sys.stdout.getvalue()

    sys.stdout = oldstdout


def test_gff3_strict():
    args = type('', (), {})()
    args.out = StringIO()
    args.gff3 = 'tests/testdata/mito-trna.gff3'
    args.strict = True

    with pytest.raises(AssertionError) as ae:
        tag.cli.gff3.main(args)
    assert 'not contained within its span' in str(ae)

    args.out == StringIO()
    args.strict = False
    tag.cli.gff3.main(args)
    testout = tag.pkgdata('mito-trna-out.gff3').read()
    assert args.out.getvalue() == testout


def test_occ():
    oldstdout = sys.stdout

    sys.stdout = StringIO()
    args = type('', (), {})()
    args.gff3 = 'tests/testdata/oluc-20kb.gff3'
    args.type = 'CDS'
    args.strict = True
    tag.cli.occ.main(args)
    assert sys.stdout.getvalue() == '14100\n'

    sys.stdout = oldstdout
