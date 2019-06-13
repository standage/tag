#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2019 Battelle National Biodefense Institute.
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import tag


class LocusBuffer(object):
    def __init__(self, firstfeature):
        self.buffer = [firstfeature]
        self.seqid = firstfeature.seqid
        self.range = firstfeature.range

    def test(self, feature):
        if feature.seqid != self.seqid:
            return False
        if not feature.range.overlap(self.range):
            return False
        return True

    def add(self, feature):
        if not self.test(feature):
            return
        self.buffer.append(feature)
        self.range = self.range.merge(feature.range)


def loci(*sorted_streams, featuretype=None):
    locus = None
    merger = tag.select.merge(*sorted_streams)
    features = tag.select.features(merger, type=featuretype)
    for feature in features:
        if locus is None:
            locus = LocusBuffer(feature)
        elif locus.test(feature) is True:
            locus.add(feature)
        else:
            yield tuple(locus.buffer)
            locus = LocusBuffer(feature)
    yield tuple(locus.buffer)
