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
from tag.tests import data_file, data_stream


def test_rice():
    reader = tag.reader.GFF3Reader(data_stream('osat-twoscaf.gff3.gz'))
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
    reader = tag.reader.GFF3Reader(data_stream('pcan-123.gff3.gz'))
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
    reader = tag.reader.GFF3Reader(data_stream('pdom-withseq.gff3'))
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
    index.consume_file(data_file('pcan-123.gff3.gz'))
    assert len(index) == 4


def test_extent():
    index = tag.index.Index()
    index.consume_file(data_file('pcan-123.gff3.gz'))
    assert index.extent('scaffold_124') == (6020, 444332)
    index.yield_inferred = False
    assert index.extent('scaffold_125') == (0, 449039)


def test_query():
    index = tag.index.Index()
    index.consume_file(data_file('osat-twoscaf.gff3.gz'))
    assert len(index.query('NW_015379189.1', 5000)) == 2
    assert len(index.query('NW_015379189.1', 5000, 15000)) == 1
    assert len(index.query('NW_015379189.1', 5000, 15000, strict=True)) == 1
    assert len(index.query('NW_015379189.1', 5000, 15000, strict=False)) == 4


def test_named_index():
    index = tag.index.NamedIndex()
    index.consume_file(data_file('pdom-withseq.gff3'))
    assert index['mRNA2'].slug == 'mRNA@PdomSCFr1.2-0483[3830, 6206]'
    with pytest.raises(IndexError):
        feature = index['mRNA3']


def test_named_index_name():
    index = tag.index.NamedIndex()
    index.consume_file(data_file('oluc-20kb.gff3'), attribute='Name')
    assert index['OSTLU_40358'].slug == 'gene@NC_009355.1[18192, 18755]'
    print(list(index.names))
    assert list(index.names) == [
        'OSTLU_23817', 'OSTLU_23818', 'OSTLU_23820', 'OSTLU_23821',
        'OSTLU_28610', 'OSTLU_28615', 'OSTLU_28616', 'OSTLU_40330',
        'OSTLU_40358', 'OSTLU_85990', 'XM_001415321.1', 'XM_001415322.1',
        'XM_001415323.1', 'XM_001415324.1', 'XM_001415325.1', 'XM_001415326.1',
        'XM_001415671.1', 'XM_001415672.1', 'XM_001415673.1', 'XM_001415674.1',
        'XP_001415358.1', 'XP_001415359.1', 'XP_001415360.1', 'XP_001415361.1',
        'XP_001415362.1', 'XP_001415363.1', 'XP_001415708.1', 'XP_001415709.1',
        'XP_001415710.1', 'XP_001415711.1'
    ]
