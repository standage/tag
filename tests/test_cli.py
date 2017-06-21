#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import pytest
import tag
import sys


def test_cli_args():
    oldstdout = sys.stdout
    oldstderr = sys.stderr

    parser = tag.cli.parser()
    sys.stdout = StringIO()
    sys.stderr = StringIO()
    with pytest.raises(SystemExit) as se:
        args = parser.parse_args(['-v'])
    assert tag.__version__ in sys.stdout.getvalue() or \
        tag.__version__ in sys.stderr.getvalue()

    parser = tag.cli.parser()
    sys.stdout = StringIO()
    sys.stderr = StringIO()
    with pytest.raises(SystemExit) as se:
        args = parser.parse_args(['gff3', '-h'])
    assert 'input file in GFF3 format' in sys.stdout.getvalue()

    sys.stdout = oldstdout
    sys.stderr = oldstderr

    parser = tag.cli.parser()
    args = parser.parse_args(['gff3', '-r', 'tests/testdata/mito-trna.gff3'])
    args.out = StringIO()
    tag.__main__.main(args)
    testout = tag.pkgdata('mito-trna-out.gff3').read()
    assert args.out.getvalue() == testout

    with pytest.raises(SystemExit):
        tag.__main__.main()


def test_gff3_strict():
    args = type('', (), {})()
    args.out = StringIO()
    args.gff3 = 'tests/testdata/mito-trna.gff3'
    args.strict = True
    args.sorted = False

    with pytest.raises(AssertionError) as ae:
        tag.cli.gff3.main(args)
    assert 'not contained within its span' in str(ae)

    args.out == StringIO()
    args.strict = False
    tag.cli.gff3.main(args)
    testout = tag.pkgdata('mito-trna-out.gff3').read()
    assert args.out.getvalue() == testout


@pytest.mark.parametrize('gff3,ftype,expected_output', [
    ('oluc-20kb.gff3', 'CDS', '14100\n'),
    ('bogus-aligns.gff3', 'cDNA_match', '7006\n'),
    ('bogus-genes.gff3', 'gene', '18000\n'),
    ('bogus-genes.gff3', 'mRNA', '18000\n'),
    ('bogus-genes.gff3', 'exon', '11000\n'),
])
def test_occ(gff3, ftype, expected_output):
    oldstdout = sys.stdout

    sys.stdout = StringIO()
    args = type('', (), {})()
    args.gff3 = 'tests/testdata/' + gff3
    args.type = ftype
    args.strict = True
    tag.cli.occ.main(args)
    assert sys.stdout.getvalue() == expected_output

    sys.stdout = oldstdout


def test_pmrna(capsys):
    args = type('', (), {})()
    args.gff3 = 'tests/testdata/nanosplice.gff3'
    args.strict = True
    tag.cli.pmrna.main(args)

    out, err = capsys.readouterr()
    exp_out = tag.pkgdata('nanosplice-primary.gff3').read()
    assert out.strip() == exp_out.strip()


def test_sum(capsys):
    infile = 'tests/testdata/GCF_001639295.1_ASM163929v1_genomic.gff.gz'
    args = tag.cli.parser().parse_args(['sum', infile])
    tag.cli.sum.main(args)

    out, err = capsys.readouterr()
    assert out.strip() == tag.pkgdata('sum-test-out.txt').read().strip()
