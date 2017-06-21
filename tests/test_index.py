#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import tag


@pytest.fixture
def pcan():
    return index


def test_rice():
    reader = tag.reader.GFF3Reader(tag.pkgdata('osat-twoscaf.gff3.gz'))
    index = tag.index.Index()
    index.consume(reader)
    assert sorted(list(index.keys())) == ['NW_015379189.1', 'NW_015379208.1']
    assert len(index['NW_015379189.1']) == 5
    assert len(index['NW_015379208.1']) == 2


def test_volvox():
    data = '##sequence-region chr1 1 {}'
    entries = [tag.directive.Directive(data.format(x)) for x in range(2)]
    index = tag.index.Index()
    with pytest.raises(ValueError) as ve:
        index.consume(entries)
    assert 'duplicate sequence region' in str(ve)


@pytest.mark.parametrize('infer,s1,s2,s3,e1,e2,e3', [
    (True, 5583, 6021, 18994, 452177, 444332, 429634),
    (False, 1, 1, 1, 463498, 452451, 449039)
])
def test_declared(infer, s1, s2, s3, e1, e2, e3):
    reader = tag.reader.GFF3Reader(tag.pkgdata('pcan-123.gff3.gz'))
    index = tag.index.Index()
    index.consume(reader)

    index.yield_inferred = infer
    sr = [e for e in index if isinstance(e, tag.directive.Directive)]
    test = list()
    for seqid, start, end in zip([123, 124, 125], [s1, s2, s3], [e1, e2, e3]):
        data = '##sequence-region scaffold_{} {} {}'.format(seqid, start, end)
        test.append(data)
    assert [str(s) for s in sorted(sr)] == test


def test_consume():
    reader = tag.reader.GFF3Reader(tag.pkgdata('pdom-withseq.gff3'))
    entries = [e for e in reader]
    assert len(entries) == 6
    index = tag.index.Index()

    with pytest.raises(ValueError) as ve:
        index.consume_seqreg(entries[0])
    assert 'expected ##sequence-region directive' in str(ve)

    with pytest.raises(ValueError) as ve:
        index.consume_seqreg(entries[4])
    assert 'expected ##sequence-region directive' in str(ve)

    with pytest.raises(ValueError) as ve:
        index.consume_feature(entries[0])
    assert 'expected Feature object' in str(ve)

    index.consume_seqreg(entries[1])
    index.consume_feature(entries[3])
    index.consume_file('tests/testdata/pcan-123.gff3.gz')
    assert len(index) == 4


def test_extent():
    index = tag.index.Index()
    index.consume_file('tests/testdata/pcan-123.gff3.gz')
    assert index.extent('scaffold_124') == (6020, 444332)
    index.yield_inferred = False
    assert index.extent('scaffold_125') == (0, 449039)
