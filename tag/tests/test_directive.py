#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
from tag.comment import Comment
from tag.directive import Directive
from tag.feature import Feature
from tag.range import Range


def test_basic():
    """Test basic object construction."""
    d = Directive('##gff-version 3')
    assert d.type == 'gff-version' and d.version == '3'
    d = Directive('##gff-version    3')  # 3 spaces
    assert d.type == 'gff-version' and d.version == '3'
    d = Directive('##gff-version	3')  # tab
    assert d.type == 'gff-version' and d.version == '3'
    assert '%r' % d == '##gff-version	3'

    with pytest.raises(AssertionError):
        d = Directive('')  # No data

    with pytest.raises(AssertionError) as ae:
        d = Directive('##gff-version   2.2')
    assert 'Only GFF version 3 is supported' in str(ae)

    with pytest.raises(AssertionError):
        d = Directive('not a directive')

    with pytest.raises(AssertionError):
        d = Directive('# Still not a directive')


def test_custom_directive():
    """Test custom directive type."""
    d1 = Directive('##bogus-directive')
    d2 = Directive('##bonus-directive   abc 1 2 3')
    d3 = Directive('##Type DNA NC_005213.1')

    assert d1.type is None and d1.dirtype == 'bogus-directive'
    assert d1.data is None
    assert d2.type is None and d2.dirtype == 'bonus-directive'
    assert d2.data == 'abc 1 2 3'
    assert d3.type is None and d3.dirtype == 'Type'
    assert d3.data == 'DNA NC_005213.1'
    assert d3.slug is None


def test_sequence_region():
    """Test sequence-region directive type."""
    r1 = Directive('##sequence-region ctg123 1 1497228')
    r2 = Directive('##sequence-region   ctg123 1 1497228')  # 3 spaces
    r3 = Directive('##sequence-region	ctg123 1 1497228')  # tab
    r4 = Directive('##sequence-region 1 1 1000')

    assert r1.type == 'sequence-region' and r1.seqid == 'ctg123' and \
        r1.range == Range(0, 1497228)
    assert r2.type == 'sequence-region' and r2.seqid == 'ctg123' and \
        r2.range == Range(0, 1497228)
    assert r3.type == 'sequence-region' and r3.seqid == 'ctg123' and \
        r3.range == Range(0, 1497228)
    assert r4.type == 'sequence-region' and r4.seqid == '1' and \
        r4.range == Range(0, 1000)
    assert r4.slug == 'sequence 1[1, 1000]'

    with pytest.raises(AssertionError) as ae:
        r5 = Directive('##sequence-region   BoGuScHr 123456 4321')
    assert '[123455, 4321] invalid, start must be <= end' in str(ae)


def test_ontology_directives():
    """Test ontology directives."""
    so_uri = ('http://song.cvs.sourceforge.net/viewvc/song/ontology/'
              'so.obo?revision=1.263')
    attr_uri = 'http://www.bogus.edu/attr-o.obo'
    src_uri = 'http://www.bogus.edu/src-o.obo'

    o1 = Directive('##feature-ontology ' + so_uri)
    o2 = Directive('##attribute-ontology   ' + attr_uri)
    o3 = Directive('##source-ontology	' + src_uri)

    assert o1.type == 'feature-ontology' and o1.uri == so_uri
    assert o2.type == 'attribute-ontology' and o2.uri == attr_uri
    assert o3.type == 'source-ontology' and o3.uri == src_uri


def test_species_directive():
    """Test species directive."""
    amel = 'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7460'
    pdom = 'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=743375'
    sinv = 'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=13686'

    s1 = Directive('##species ' + amel)
    s2 = Directive('##species   ' + pdom)
    s3 = Directive('##species	' + sinv)

    assert s1.type == 'species' and s1.uri == amel
    assert s2.type == 'species' and s2.uri == pdom
    assert s3.type == 'species' and s3.uri == sinv


def test_genome_build_directive():
    """Test genome-build directive."""
    b1 = Directive('##genome-build NCBI B36')
    b2 = Directive('##genome-build   WormBase ws110')
    b3 = Directive('##genome-build	FlyBase  r4.1')

    assert b1.type == 'genome-build' and b1.source == 'NCBI' and \
        b1.build_name == 'B36'
    assert b2.type == 'genome-build' and b2.source == 'WormBase' and \
        b2.build_name == 'ws110'
    assert b3.type == 'genome-build' and b3.source == 'FlyBase' and \
        b3.build_name == 'r4.1'


def test_sorting():
    """Test sorting and comparison"""
    gv = Directive('##gff-version   3')
    sr1 = Directive('##sequence-region chr1 500 2000')
    sr2 = Directive('##sequence-region chr1 3000 4000')
    sr3 = Directive('##sequence-region chr2 500 2000')
    sr4 = Directive('##sequence-region chr10 500 2000')
    d1 = Directive('##bonus-directive   abc 1 2 3')
    s1 = Directive('##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/'
                   'wwwtax.cgi?id=7460')
    f1 = Feature('chr', 'mRNA', 1000, 1420, strand='+')
    c1 = Comment('# The quick brown fox jumps over the lazy dog.')

    for record in [sr1, d1, s1, f1, c1]:
        assert gv < record
        assert gv <= record
        assert not gv > record
        assert not gv >= record

    assert sr1 > gv
    assert not sr1 < gv
    for record in [d1, s1, f1, c1]:
        assert sr1 < record
        assert sr1 <= record
        assert not sr1 > record
        assert not sr1 >= record
    assert sorted([sr4, sr1, sr3, sr2]) == [sr1, sr2, sr4, sr3], '%r %r' % \
        (sorted([sr4, sr1, sr3, sr2]), [sr1, sr2, sr4, sr3])
    assert sr1 > gv
    assert sr1 >= gv
    assert sr1 <= sr2
    assert sr1 <= sr3
    assert not sr1 >= sr2
    assert not sr1 >= sr3
    assert d1 < s1
    assert d1 <= s1
    assert d1 < f1
    assert d1 <= f1
    assert s1 < c1
    assert s1 <= c1

    for record in [s1, f1, c1]:
        assert d1 < record
