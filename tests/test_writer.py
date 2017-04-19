#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import pytest
import tag
from tag import GFF3Reader
from tag import GFF3Writer


def test_write_grape(capsys):
    reader = GFF3Reader(infilename='tests/testdata/grape-cpgat.gff3')
    writer = GFF3Writer(reader)
    writer.write()

    out, err = capsys.readouterr()
    testout = tag.pkgdata('grape-cpgat-writer.gff3').read()
    assert out.strip() == testout.strip()


def test_write_stdout():
    reader = GFF3Reader(infilename='tests/testdata/grape-cpgat.gff3')
    writer = GFF3Writer(reader)
    writer.write()


def test_write_file():
    reader = GFF3Reader(infilename='tests/testdata/grape-cpgat.gff3')
    output = tag.pkgdata('grape-cpgat-writer-out.gff3', 'w')
    writer = GFF3Writer(reader, output)
    writer.write()
    del writer

    obs_out = tag.pkgdata('grape-cpgat-writer-out.gff3').read()
    exp_out = tag.pkgdata('grape-cpgat-writer.gff3').read()
    assert obs_out.strip() == exp_out.strip()


@pytest.mark.parametrize('gff3', [
    'minimus.gff3',
    'prokka.gff3',
])
def test_write_in_out(gff3, capsys):
    reader = GFF3Reader(tag.pkgdata(gff3))
    writer = GFF3Writer(reader)
    writer.write()

    out, err = capsys.readouterr()
    assert out.strip() == tag.pkgdata(gff3).read().strip()


def test_sort_multifeat(capsys):
    reader = GFF3Reader(tag.pkgdata('psyllid-cdnamatch.gff3'))
    writer = GFF3Writer(reader)
    writer.write()

    out, err = capsys.readouterr()
    testout = tag.pkgdata('psyllid-cdnamatch-sorted.gff3').read()
    assert out.strip() == testout.strip()


def test_sort_multifeat_reverse(capsys):
    reader = GFF3Reader(tag.pkgdata('psyllid-cdnamatch-reverse.gff3'))
    writer = GFF3Writer(reader)
    writer.write()

    out, err = capsys.readouterr()
    testout = tag.pkgdata('psyllid-cdnamatch-reverse-sorted.gff3').read()
    assert out.strip() == testout.strip()
