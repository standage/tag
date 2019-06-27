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

    def test(self, feature, minbp=25, minperc=0.25):
        if feature.seqid != self.seqid:
            return False
        insufficient_overlap = not feature.range.overlap_atleast(
            self.range, minbp=minbp, minperc=minperc
        )
        if insufficient_overlap:
            return False
        return True

    def add(self, feature):
        self.buffer.append(feature)
        self.range = self.range.merge(feature.range)


def loci(*sorted_streams, **kwargs):
    """Determine feature loci from two or more sorted annotation streams.

    Rather than simply relying on gene coordinates from a reference annotation,
    this function will determine coordinates for features of interest by
    considering an arbitrary number of annotation sources, whether they be
    references or predictions. This enables efficient and accurate assessment
    of both sensitivity and specificity at the level of individual nucleotides
    and entire features.
    """
    featuretype = kwargs['featuretype'] if 'featuretype' in kwargs else None
    minbp = kwargs['minbp'] if 'minbp' in kwargs else 25
    minperc = kwargs['minperc'] if 'minperc' in kwargs else 0.25
    locus = None
    merger = tag.select.merge(*sorted_streams)
    features = tag.select.features(merger, type=featuretype)
    for feature in features:
        if locus is None:
            locus = LocusBuffer(feature)
        elif locus.test(feature, minbp=minbp, minperc=minperc) is True:
            locus.add(feature)
        else:
            yield locus.seqid, locus.range, tuple(locus.buffer)
            locus = LocusBuffer(feature)
    yield locus.seqid, locus.range, tuple(locus.buffer)


def pocus(*sorted_streams, **kwargs):
    """Feature stream for locus parsing.

    From two or more sorted annotation streams, create a new stream that yields
    the features from the original streams, additional `locus` features, and
    separator directives (`###`) between loci.
    """
    delta = kwargs['delta'] if 'delta' in kwargs else 0
    for seqid, rng, features in loci(*sorted_streams, **kwargs):
        locus = tag.Feature(
            seqid, 'experimental_feature', rng.start - delta, rng.end + delta,
            source='tag::locus::pocus',
        )
        locus.add_attribute('Note', 'Interval locus')
        yield locus
        for feature in features:
            yield feature
        yield tag.Directive('###')
