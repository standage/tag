#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# -----------------------------------------------------------------------------
import pytest


class Range(object):
    """
    Represents a region of a biological sequence.

    Start and end coordinates correspond to the 0-based half-open interval:
    start inclusive, end exclusive (start=0 and end=10 corresonds to the first
    10 nucleotides).
    """

    def __init__(self, start, end):
        # Sanity checks

        assert start >= 0, ('start coordinate {} invalid, must be an '
                            'integer > 0'.format(start))
        assert end >= 0, ('end coordinate {} invalid,  must be an '
                          'integer > 0'.format(start))
        assert start <= end, ('coordinates [{}, {}] invalid, start must be '
                              '<= end'.format(start, end))
        self._start = start
        self._end = end

    def __str__(self):
        """String representation of a range"""
        if self._start == self._end:
            return str(self._start)
        return '{}-{}'.format(self._start, self._end)

    def __repr__(self):
        return str(self)

    def __len__(self):
        """Length of a range."""
        return self._end - self._start

    def __eq__(self, other):
        """
        Rich comparison operator for Python 3 support.

        Any object that defines a custom `__eq__()` opterator without also
        defining a custom `__hash__()` operator is unhashable. At least for
        now, I don't see this being a problem.
        """
        return self._start == other.start and self._end == other.end

    def __lt__(self, other):
        """Rich comparison operator for Python 3 support."""
        return self._start < other.start or (self._start == other.start and
                                             self._end < other.end)

    def __le__(self, other):
        """Rich comparison operator for Python 3 support."""
        return self.__eq__(other) or self.__lt__(other)

    def __gt__(self, other):
        """Rich comparison operator for Python 3 support."""
        return self._start > other.start or (self._start == other.start and
                                             self._end > other.end)

    def __ge__(self, other):
        """Rich comparison operator for Python 3 support."""
        return self.__eq__(other) or self.__gt__(other)

    @property
    def start(self):
        return self._start

    @start.setter
    def start(self, newstart):
        assert newstart > 0, ('new start coordinate %d invalid, must be an '
                              'integer > 0' % newstart)
        assert newstart <= self._end, ('new start coordinate %d invalid, must '
                                       'be <= end %d' % (newstart, self._end))
        self._start = newstart

    @property
    def end(self):
        return self._end

    @end.setter
    def end(self, newend):
        assert newend > 0, ('new end coordinate %d invalid, must be an '
                            'integer > 0' % newend)
        assert self._start <= newend, ('new end coordinate %d is invalid, '
                                       'must be >= start %d' %
                                       (newend, self._start))
        self._end = newend

    def merge(self, other):
        """
        Merge two ranges.

        Take the union of two ranges and create a new range.
        """
        newstart = min(self._start, other.start)
        newend = max(self._end, other.end)
        return Range(newstart, newend)

    def intersect(self, other):
        """
        Intersect two ranges.

        Calculate the overlap, if any, between two ranges and create a new
        range.
        """
        if not self.overlap(other):
            return None

        newstart = max(self._start, other.start)
        newend = min(self._end, other.end)
        return Range(newstart, newend)

    def overlap(self, other):
        """Determine whether two ranges overlap."""
        if self._start < other.end and self._end > other.start:
            return True
        return False

    def contains(self, other):
        """Determine whether the range contains another range."""
        if self._start <= other.start and self._end >= other.end:
            return True
        return False

    def within(self, other):
        """Determine whether the range falls completely within another."""
        return other.contains(self)

    def transform(self, offset):
        """Transform the range's coordinates using the given offset."""
        assert self._start + offset > 0, \
            ('offset %d invalid; resulting range [%d, %d] is undefined' %
             (offset, self._start + offset, self._end + offset))
        self._start += offset
        self._end += offset


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------

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
    assert '%r' % r1 == '0-10'
    assert '%r' % r2 == '1234-5678'
    assert '%r' % r3 == '528-901'
    assert '%r' % r4 == '42'
    assert '{}'.format(r1) == '0-10'
    assert '{}'.format(r2) == '1234-5678'
    assert '{}'.format(r3) == '528-901'
    assert '{}'.format(r4) == '42'


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
