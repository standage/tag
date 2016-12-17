#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
from intervaltree import IntervalTree
import tag
from tag.directive import Directive


class Index(defaultdict):

    def __init__(self):
        defaultdict.__init__(self, IntervalTree)
        self.declared_regions = dict()
        self.inferred_regions = dict()
        self.yield_inferred = True

    def consume(self, entrystream):
        for entry in entrystream:
            if isinstance(entry, tag.directive.Directive) and \
                    entry.type == 'sequence-region':
                if entry.seqid in self.declared_regions:
                    msg = 'duplicate sequence region "{}"'.format(entry.seqid)
                    raise ValueError(msg)
                self.declared_regions[entry.seqid] = entry.range
            elif isinstance(entry, tag.feature.Feature):
                self[entry.seqid][entry.start:entry.end] = entry
                if entry.seqid not in self.inferred_regions:
                    self.inferred_regions[entry.seqid] = tag.range.Range(
                        entry.start, entry.end
                    )
                newstart = min(entry.start,
                               self.inferred_regions[entry.seqid].start)
                newend = max(entry.end, self.inferred_regions[entry.seqid].end)
                self.inferred_regions[entry.seqid].start = newstart
                self.inferred_regions[entry.seqid].end = newend

    @property
    def outstream(self):
        regions = self.inferred_regions
        if not self.yield_inferred:
            regions = self.declared_regions

        for seqid in sorted(list(self.keys())):
            sr = regions[seqid]
            data = '##sequence-region {} {} {}'.format(seqid, sr.start, sr.end)
            yield Directive(data)

        for seqid in sorted(list(self.keys())):
            for interval in self[seqid]:
                yield interval.data
