#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import division
from itertools import combinations
from networkx import Graph, connected_components


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
    >>> rng.contains(Range(10, 20))
    True
    >>> rng.contains_point(2345)
    False
    >>> rng.start = 500
    >>> rng.start
    500
    >>> rng.end
    1000
    """

    @staticmethod
    def merge_overlapping(ranges):
        if len(ranges) < 2:
            for rng in ranges:
                yield rng
        elif len(ranges) == 2 and ranges[0].overlap(ranges[1]):
            newstart = min([r.start for r in ranges])
            newend = max([r.end for r in ranges])
            yield Range(newstart, newend)
        elif len(ranges) == 2 and not ranges[0].overlap(ranges[1]):
            for rng in ranges:
                yield rng
        else:
            graph = Graph()
            for rng in ranges:
                graph.add_node(rng)
            for r1, r2 in combinations(ranges, 2):
                if r1.overlap_extent(r2) > 0:
                    graph.add_edge(r1, r2)
            for cc in connected_components(graph):
                newstart = min([r.start for r in cc])
                newend = max([r.end for r in cc])
                yield Range(newstart, newend)

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
        return self._start == other._start and self._end == other._end

    def __ne__(self, other):
        return self._start != other._start or self._end != other._end

    def __hash__(self):
        return hash((self.start, self.end))

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
        """Determine the interval of overlap between this range and another.

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

    def overlap_extent(self, other):
        """Compute number of nucleotides of overlap between two ranges."""
        if not self.overlap(other):
            return 0
        coords = sorted((self._start, other._start, self._end, other._end))
        return coords[2] - coords[1]

    def overlap_atleast(self, other, minbp=1, minperc=0.0):
        """Determine whether two ranges overlap by at least some threshold.

        The thresholds can be specified using a minimum number of bases, or a
        minimum percentage of the length of the ranges, or both.
        """
        assert minbp >= 1, 'must require at least 1bp overlap'
        bpoverlap = self.overlap_extent(other)
        if minbp and bpoverlap < minbp:
            return False
        if minperc and bpoverlap / len(self) < minperc:
            return False
        if minperc and bpoverlap / len(other) < minperc:
            return False
        return True

    def contains(self, other):
        """Determine whether this range contains another."""
        return self._start <= other.start and self._end >= other.end

    def contains_point(self, point):
        """Determine whether this range contains the specified point."""
        return self._start <= point and self._end >= point

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

    def copy(self):
        return Range(self._start, self._end)
