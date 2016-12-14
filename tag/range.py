#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------


class Range(object):
    """
    Represents a region of a biological sequence.

    Start and end coordinates correspond to the 0-based half-open interval:
    start inclusive, end exclusive (start=0 and end=10 corresonds to the first
    10 nucleotides).

    >>> rng = Range(0, 1000)
    >>> rng
    [0, 1000)
    >>> len(rng)
    1000
    >>> rng.start = 500
    >>> rng.start
    500
    >>> rng.end
    1000
    """

    def __init__(self, start, end):
        assert start >= 0, ('start coordinate {} invalid, must be an '
                            'integer >= 0'.format(start))
        assert end >= 0, ('end coordinate {} invalid,  must be an '
                          'integer >= 0'.format(end))
        assert start <= end, ('coordinates [{}, {}] invalid, start must be '
                              '<= end'.format(start, end))
        self._start = start
        self._end = end

    def __str__(self):
        if self._start == self._end:
            return str(self._start)
        return '[{}, {})'.format(self._start, self._end)

    def __repr__(self):
        return str(self)

    def __len__(self):
        return self._end - self._start

    # Begin rich comparison operators for Python 3 support
    def __eq__(self, other):
        """
        Rich comparison operator for Python 3 support.

        Any object that defines a custom `__eq__()` operator without also
        defining a custom `__hash__()` operator is unhashable. At least for
        now, I don't see this being a problem.
        """
        return self._start == other.start and self._end == other.end

    def __lt__(self, other):
        return self._start < other.start or (self._start == other.start and
                                             self._end < other.end)

    def __le__(self, other):
        return self == other or self < other

    def __gt__(self, other):
        return self._start > other.start or (self._start == other.start and
                                             self._end > other.end)

    def __ge__(self, other):
        return self == other or self > other
    # End rich comparison operators for Python 3 support

    @property
    def start(self):
        return self._start

    @start.setter
    def start(self, newstart):
        assert newstart >= 0, \
            ('new start coordinate {} invalid, must be an '
             'integer >= 0'.format(newstart))
        assert newstart <= self._end, \
            ('new start coordinate {} invalid, must '
             'be <= end {}'.format(newstart, self._end))
        self._start = newstart

    @property
    def end(self):
        return self._end

    @end.setter
    def end(self, newend):
        assert newend >= 0, \
            ('new end coordinate {} invalid, must be an '
             'integer >= 0'.format(newend))
        assert self._start <= newend, \
            ('new end coordinate {} is invalid, '
             'must be >= start {}'.format(newend, self._start))
        self._end = newend

    def merge(self, other):
        """
        Merge this range object with another (ranges need not overlap or abut).

        :returns: a new Range object representing the interval containing both
                  ranges.
        """
        newstart = min(self._start, other.start)
        newend = max(self._end, other.end)
        return Range(newstart, newend)

    def intersect(self, other):
        """
        Determine the interval of overlap between this range and another.

        :returns: a new Range object representing the overlapping interval,
                  or `None` if the ranges do not overlap.
        """
        if not self.overlap(other):
            return None

        newstart = max(self._start, other.start)
        newend = min(self._end, other.end)
        return Range(newstart, newend)

    def overlap(self, other):
        """Determine whether this range overlaps with another."""
        if self._start < other.end and self._end > other.start:
            return True
        return False

    def contains(self, other):
        """Determine whether this range contains another."""
        if self._start <= other.start and self._end >= other.end:
            return True
        return False

    def within(self, other):
        """Determine whether this range is contained within another."""
        return other.contains(self)

    def transform(self, offset):
        """
        Shift this range by the specified offset.

        Note: the resulting range must be a valid interval.
        """
        assert self._start + offset > 0, \
            ('offset {} invalid; resulting range [{}, {}) is '
             'undefined'.format(offset, self._start+offset, self._end+offset))
        self._start += offset
        self._end += offset
