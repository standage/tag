#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the ISC license: see LICENSE.txt.
# -----------------------------------------------------------------------------
"""Representation of a region of a biological sequence."""


class Region(object):
    """Start and end coordinates are 1-based, as in GFF3."""

    def __init__(self, start, end):
        # Sanity checks
        assert start > 0, ('start coordinate %d invalid, must be an '
                           'integer > 0' % start)
        assert end > 0, ('end coordinate %d invalid,  must be an '
                         'integer > 0' % start)
        assert start <= end, ('coordinates [%d, %d] invalid, start must be '
                              '<= end' % (start, end))
        self._start = start
        self._end = end

    def __str__(self):
        """String representation of a region"""
        if self._start == self._end:
            return '%d' % self._start
        return '%d-%d' % (self._start, self._end)

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        """Length of a region."""
        return self._end - self._start + 1

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
        Merge two regions.

        Take the union of two regions and create a new region.
        """
        newstart = min(self._start, other.start)
        newend = max(self._end, other.end)
        return Region(newstart, newend)

    def intersect(self, other):
        """
        Intersect two regions.

        Calculate the overlap, if any, between two regions and create a new
        region.
        """
        if not self.overlap(other):
            return None

        newstart = max(self._start, other.start)
        newend = min(self._end, other.end)
        return Region(newstart, newend)

    def overlap(self, other):
        """Determine whether two regions overlap."""
        if self._start <= other.end and self._end >= other.start:
            return True
        return False

    def contains(self, other):
        """Determine whether the region contains another region."""
        if self._start <= other.start and self._end >= other.end:
            return True
        return False

    def within(self, other):
        """Determine whether the region falls completely within another."""
        return other.contains(self)

    def transform(self, offset):
        """Transform the region's coordinates using the given offset."""
        assert self._start + offset > 0, \
            ('offset %d invalid; resulting region [%d, %d] is undefined' %
             (offset, self._start + offset, self._end + offset))
        self._start += offset
        self._end += offset


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------

def test_repr():
    """[aeneas::Region] Test string representation of regions."""
    r1 = Region(1, 10)
    r2 = Region(1234, 5678)
    r3 = Region(528, 901)
    r4 = Region(42, 42)

    assert '%s' % r1 == '1-10'
    assert '%r' % r1 == '1-10'
    assert '%s' % r2 == '1234-5678'
    assert '%r' % r2 == '1234-5678'
    assert '%s' % r3 == '528-901'
    assert '%r' % r3 == '528-901'
    assert '%s' % r4 == '42'
    assert '%r' % r4 == '42'


def test_len():
    """[aeneas::Region] Test length."""
    r1 = Region(1, 10)
    r2 = Region(1234, 5678)
    r3 = Region(528, 901)
    r4 = Region(42, 42)

    assert len(r1) == 10
    assert len(r2) == 4445
    assert len(r3) == 374
    assert len(r4) == 1


def test_cmp_sort():
    """[aeneas::Region] Test comparison and sorting."""
    r1 = Region(1, 10)
    r2 = Region(1234, 5678)
    r3 = Region(528, 901)
    r4 = Region(42, 42)

    assert r3.__eq__(r3)
    assert r3.__lt__(r2)
    assert r3.__le__(r2)
    assert r3.__gt__(r4)
    assert r3.__ge__(r4)

    assert sorted([r1, r2, r3, r4]) == [r1, r4, r3, r2]


def test_getters():
    """[aeneas::Region] Test getters and setters."""
    r1 = Region(1, 1)

    r1.end = 42
    assert r1.start == 1
    assert r1.end == 42

    r1.start = 42
    assert r1.start == 42
    assert r1.end == 42

    try:
        r1.start = 43
    except AssertionError as ae:
        pass
    assert r1 == Region(42, 42)


def test_merge():
    """[aeneas::Region] Test region merge."""
    assert Region(1, 10).merge(Region(11, 20)) == Region(1, 20)
    assert Region(100, 1000).merge(Region(2000, 2500)) == Region(100, 2500)
    assert Region(1, 10).merge(Region(9, 15)) == Region(1, 15)


def test_intersect():
    """[aeneas::Region] Test region intersect."""
    assert Region(1, 10).intersect(Region(11, 20)) is None
    assert Region(100, 1000).intersect(Region(2000, 2500)) is None
    assert Region(1, 10).intersect(Region(9, 15)) == Region(9, 10)
    assert Region(1, 10).intersect(Region(10, 15)) == Region(10, 10)
    assert Region(279886, 283581).intersect(Region(280065, 297216)) == \
        Region(280065, 283581)


def test_overlap():
    """[aeneas::Region] Test region overlap."""
    assert Region(1, 10).overlap(Region(11, 20)) == False
    assert Region(100, 1000).overlap(Region(2000, 2500)) == False
    assert Region(1, 10).overlap(Region(9, 15)) == True
    assert Region(1, 10).overlap(Region(10, 15)) == True
    assert Region(279886, 283581).overlap(Region(280065, 297216)) == True


def test_contains_within():
    """[aeneas::Region] Containment tests."""
    assert Region(1, 10).contains(Region(1, 10)) == True
    assert Region(1, 10).within(Region(1, 10)) == True
    assert Region(1, 10).contains(Region(2, 5)) == True
    assert Region(2, 5).contains(Region(1, 10)) == False
    assert Region(279886, 283581).contains(Region(280065, 297216)) == False
    assert Region(279886, 283581).contains(Region(279986, 283481)) == True
    assert Region(279986, 283481).contains(Region(279886, 283581)) == False
    assert Region(279986, 283481).within(Region(279886, 283581)) == True


def test_transform():
    """[aeneas::Region] Transformation tests."""
    r1 = Region(1, 10)
    r2 = Region(2000, 3000)
    r3 = Region(100, 200)

    r1.transform(100)
    assert r1 == Region(101, 110)
    r2.transform(-1000)
    assert r2 == Region(1000, 2000)
    try:
        r3.transform(-300)
    except AssertionError as ae:
        pass
    assert r3 == Region(100, 200)
