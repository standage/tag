#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
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
from tag import GFF3Reader
from tag import GFF3Writer


def test_write_grape():
    reader = GFF3Reader(infilename='tests/testdata/grape-cpgat.gff3')
    output = StringIO()
    writer = GFF3Writer(reader, output)
    writer.write()
    with open('tests/testdata/grape-cpgat-writer.gff3', 'r') as testout:
        testoutput = testout.read()
    assert output.getvalue() == testoutput, (output, testoutput)


def test_write_stdout():
    reader = GFF3Reader(infilename='tests/testdata/grape-cpgat.gff3')
    writer = GFF3Writer(reader)
    writer.write()


def test_write_file():
    reader = GFF3Reader(infilename='tests/testdata/grape-cpgat.gff3')
    output = 'tests/testdata/grape-cpgat-writer-out.gff3'
    writer = GFF3Writer(reader, output)
    writer.write()
    del writer

    with open('tests/testdata/grape-cpgat-writer.gff3', 'r') as testout:
        testoutput1 = testout.read()
    with open(output, 'r') as testout:
        testoutput2 = testout.read()
    assert testoutput1 == testoutput2, (testoutput1, testoutput2)


def test_write_minimus():
    reader = GFF3Reader(tag.pkgdata('minimus.gff3'))
    output = StringIO()
    writer = GFF3Writer(reader, output)
    writer.write()

    assert output.getvalue() == tag.pkgdata('minimus.gff3').read()
