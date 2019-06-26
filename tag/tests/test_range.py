#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
from tag import Range
from tag.reader import GFF3Reader
from tag.tests import data_file, data_stream


def test_repr():
    """Test string representation of ranges."""
    r1 = Range(0, 10)
    r2 = Range(1234, 5678)
    r3 = Range(528, 901)
    r4 = Range(42, 42)

    assert str(r1) == '[0, 10)'
    assert str(r2) == '[1234, 5678)'
    assert str(r3) == '[528, 901)'
    assert str(r4) == '42'
    assert repr(r1) == '[0, 10)'
    assert repr(r2) == '[1234, 5678)'
    assert repr(r3) == '[528, 901)'
    assert repr(r4) == '42'


def test_len():
    """Test length."""
    assert len(Range(0, 10)) == 10
    assert len(Range(1234, 5678)) == 4444
    assert len(Range(528, 901)) == 373
    assert len(Range(42, 43)) == 1
    assert len(Range(42, 42)) == 0


def test_cmp_sort():
    """Test comparison and sorting."""
    r1 = Range(0, 10)
    r2 = Range(1234, 5678)
    r3 = Range(528, 901)
    r4 = Range(42, 42)

    assert r3.__eq__(r3)
    assert r3.__lt__(r2)
    assert r3.__le__(r2)
    assert r3.__gt__(r4)
    assert r3.__ge__(r4)

    assert hash(r3) == hash(Range(528, 901))
    assert hash(r3) != hash(r4)

    assert sorted([r1, r2, r3, r4]) == [r1, r4, r3, r2]


def test_getters():
    """Test getters and setters."""
    r1 = Range(0, 1)

    r1.end = 42
    assert r1.start == 0
    assert r1.end == 42

    r1.start = 42
    assert r1.start == 42
    assert r1.end == 42

    with pytest.raises(AssertionError) as ae:
        r1.start = 43
    assert 'new start coordinate 43 invalid' in str(ae)
    assert r1 == Range(42, 42)


def test_merge():
    """Test range merge."""
    assert Range(0, 10).merge(Range(11, 20)) == Range(0, 20)
    assert Range(100, 1000).merge(Range(2000, 2500)) == Range(100, 2500)
    assert Range(0, 10).merge(Range(9, 15)) == Range(0, 15)


def test_intersect():
    """Test range intersect."""
    assert Range(0, 10).intersect(Range(11, 20)) is None
    assert Range(0, 10).intersect(Range(10, 20)) is None
    assert Range(0, 10).intersect(Range(9, 20)) == Range(9, 10)
    assert Range(100, 1000).intersect(Range(2000, 2500)) is None
    assert Range(279886, 283581).intersect(Range(280065, 297216)) == \
        Range(280065, 283581)


def test_overlap():
    """Test range overlap."""
    assert Range(0, 10).overlap(Range(11, 20)) is False
    assert Range(0, 10).overlap(Range(10, 15)) is False
    assert Range(0, 10).overlap(Range(9, 15)) is True
    assert Range(100, 1000).overlap(Range(2000, 2500)) is False
    assert Range(279886, 283581).overlap(Range(280065, 297216)) is True


def test_overlap_extent():
    assert Range(0, 100).overlap_extent(Range(200, 300)) == 0
    assert Range(0, 100).overlap_extent(Range(90, 200)) == 10
    assert Range(20, 40).overlap_extent(Range(0, 100)) == 20


def test_overlap_atleast():
    with pytest.raises(AssertionError) as ae:
        Range(0, 100).overlap_atleast(Range(200, 300), minbp=0)
    assert 'must require at least 1bp overlap' in str(ae)
    assert Range(0, 100).overlap_atleast(Range(200, 300), minbp=1) is False
    assert Range(0, 100).overlap_atleast(Range(90, 200), minbp=1) is True
    assert Range(0, 100).overlap_atleast(Range(90, 200), minbp=10) is True
    assert Range(0, 100).overlap_atleast(Range(90, 200), minbp=11) is False
    assert Range(0, 100).overlap_atleast(Range(20, 120), minperc=0.5) is True
    assert Range(0, 100).overlap_atleast(Range(20, 120), minperc=0.95) is False
    assert Range(850, 1000).overlap_atleast(
        Range(900, 2000), minperc=0.25
    ) is False


def test_contains_within():
    """Containment tests."""
    assert Range(0, 10).contains(Range(0, 10)) is True
    assert Range(0, 10).within(Range(0, 10)) is True
    assert Range(0, 10).contains(Range(2, 5)) is True
    assert Range(2, 5).contains(Range(0, 10)) is False
    assert Range(279886, 283581).contains(Range(280065, 297216)) is False
    assert Range(279886, 283581).contains(Range(279986, 283481)) is True
    assert Range(279986, 283481).contains(Range(279886, 283581)) is False
    assert Range(279986, 283481).within(Range(279886, 283581)) is True
    assert Range(9876, 9999).contains_point(10) is False
    assert Range(9876, 9999).contains_point(9900) is True


def test_transform():
    """Transformation tests."""
    r1 = Range(0, 10)
    r2 = Range(2000, 3000)
    r3 = Range(100, 200)

    r1.transform(100)
    r2.transform(-1000)
    with pytest.raises(AssertionError) as ae:
        r3.transform(-300)
    assert 'offset -300 invalid' in str(ae)

    assert r1 == Range(100, 110)
    assert r2 == Range(1000, 2000)


def test_amel():
    reader = GFF3Reader(infilename=data_file('amel-cdna-multi.gff3'))
    records = [r for r in reader]
    assert len(records) == 4
    assert records[2]._range == Range(263429, 325784)
    assert records[2].type == 'gene'
    assert records[3]._range == Range(263429, 325784)
    assert records[3].type == 'cDNA_match'


def test_merge_overlapping():
    assert list(Range.merge_overlapping([])) == []
    assert list(Range.merge_overlapping([Range(0, 10)])) == [Range(0, 10)]

    r1 = [Range(1000, 2000), Range(3000, 4000)]
    assert list(Range.merge_overlapping(r1)) == r1
    r2 = [Range(1000, 2000), Range(1900, 2500)]
    assert list(Range.merge_overlapping(r2)) == [Range(1000, 2500)]
    r3 = [Range(10000, 20000), Range(15000, 16000)]
    assert list(Range.merge_overlapping(r3)) == [Range(10000, 20000)]

    r4 = [Range(20, 40), Range(50, 70), Range(80, 95)]
    assert sorted(Range.merge_overlapping(r4)) == r4
    r5 = [Range(5, 50), Range(40, 60), Range(80, 95)]
    assert sorted(Range.merge_overlapping(r5)) == [Range(5, 60), Range(80, 95)]
    r6 = [
        Range(16000, 18000), Range(17050, 19000), Range(18000, 18500),
        Range(20100, 21000), Range(22000, 24000), Range(23750, 29000),
    ]
    r6result = [Range(16000, 19000), Range(20100, 21000), Range(22000, 29000)]
    assert sorted(Range.merge_overlapping(r6)) == r6result
