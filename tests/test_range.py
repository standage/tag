#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import pytest
from aeneas.range import Range


def test_repr():
    """Test string representation of ranges."""
    r1 = Range(0, 10)
    r2 = Range(1234, 5678)
    r3 = Range(528, 901)
    r4 = Range(42, 42)

    assert str(r1) == '0-10'
    assert str(r2) == '1234-5678'
    assert str(r3) == '528-901'
    assert str(r4) == '42'
    assert repr(r1) == '0-10'
    assert repr(r2) == '1234-5678'
    assert repr(r3) == '528-901'
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
