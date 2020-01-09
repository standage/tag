#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import tag
from tag import GFF3Reader
from tag import GFF3Writer
from tag.tests import data_file, data_stream
from tempfile import NamedTemporaryFile


def test_write_grape(capsys):
    reader = GFF3Reader(infilename=data_file('grape-cpgat.gff3'))
    writer = GFF3Writer(reader)
    writer.write()
    terminal = capsys.readouterr()
    testout = data_stream('grape-cpgat-writer.gff3').read()
    assert terminal.out.strip() == testout.strip()


def test_write_grape_retainids():
    with NamedTemporaryFile(suffix='.gff3', mode='w+t') as outfile:
        reader = GFF3Reader(infilename=data_file('grape-cpgat-unsorted.gff3'))
        writer = GFF3Writer(reader, outfile=outfile)
        writer.retainids = True
        writer.write()
        outfile.seek(0)
        output = outfile.read().strip()
    testfile = data_file('grape-cpgat-unsorted-retainids-out.gff3')
    with open(testfile, 'r') as fh:
        testout = fh.read().strip()
    assert output == testout


def test_write_stdout(capsys):
    reader = GFF3Reader(infilename=data_file('grape-cpgat.gff3'))
    writer = GFF3Writer(reader)
    writer.write()
    terminal = capsys.readouterr()
    assert 'CpGAT\tCDS\t4916\t5081' in terminal.out


def test_write_file():
    reader = GFF3Reader(infilename=data_file('grape-cpgat.gff3'))
    with NamedTemporaryFile(suffix='.gff3', mode='w+t') as outfile:
        writer = GFF3Writer(reader, outfile)
        writer.write()
        outfile.seek(0)
        obs_out = outfile.read()
    exp_out = data_stream('grape-cpgat-writer.gff3').read()
    assert obs_out.strip() == exp_out.strip()


def test_no_complex_separators():
    reader = GFF3Reader(infilename=data_file('grape-cpgat.gff3'))
    with NamedTemporaryFile(suffix='.gff3', mode='w+t') as outfile:
        writer = GFF3Writer(reader, outfile)
        writer.complex_separators = False
        writer.write()
        outfile.seek(0)
        obs_out = outfile.read()
    exp_out = data_stream('grape-cpgat-no-sep.gff3').read()
    assert obs_out.strip() == exp_out.strip()


@pytest.mark.parametrize('gff3', [
    'minimus.gff3',
    'prokka.gff3',
])
def test_write_in_out(gff3, capsys):
    reader = GFF3Reader(data_stream(gff3))
    writer = GFF3Writer(reader)
    writer.write(blockitvl=10)
    terminal = capsys.readouterr()
    assert terminal.out.strip() == data_stream(gff3).read().strip()


def test_sort_multifeat(capsys):
    reader = GFF3Reader(data_stream('psyllid-cdnamatch.gff3'))
    writer = GFF3Writer(reader)
    writer.write()
    terminal = capsys.readouterr()
    testout = data_stream('psyllid-cdnamatch-sorted.gff3').read()
    assert terminal.out.strip() == testout.strip()


def test_sort_multifeat_reverse(capsys):
    reader = GFF3Reader(data_stream('psyllid-cdnamatch-reverse.gff3'))
    writer = GFF3Writer(reader)
    writer.write()
    terminal = capsys.readouterr()
    testout = data_stream('psyllid-cdnamatch-reverse-sorted.gff3').read()
    assert terminal.out.strip() == testout.strip()
