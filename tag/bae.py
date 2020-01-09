#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2019 Battelle National Biodefense Institute.
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import division
from collections import defaultdict
import tag


def encodes_cds(feature):
    for subfeature in feature:
        if subfeature.type == 'CDS':
            return True
    return False


def collapse_locus(features):
    """Collapse redundant annotation information.

    Given a locus, collapse any gene predictions with identical coordinates
    into a single prediction, noting the annotation sources that support the
    collapsed gene model.
    """
    features_to_keep = list()
    cds_by_coord = defaultdict(list)
    for cds in tag.select.features(features, type='CDS', traverse=True):
        cds_by_coord[cds.range].append(cds)

    for coord in sorted(cds_by_coord):
        cdss = cds_by_coord[coord]
        support = ','.join([c.source for c in cdss])
        newcds = tag.Feature(
            cdss[0].seqid, 'CDS', coord.start, coord.end, source='tag::bae',
            score=len(cdss), strand=cdss[0].strand,
            attrstr='support_from=' + support
        )
        features_to_keep.append(newcds)

    for feature in features:
        if not encodes_cds(feature):
            features_to_keep.append(feature)

    for feature in sorted(features_to_keep):
        yield feature


def collapse_stream(locusstream):
    """Feature stream for collapsing bacterial annotation data."""
    for seqid, interval, locus in locusstream:
        for feature in collapse_locus(locus):
            yield feature
        yield tag.Directive('###')


def eval_locus(features, weights=[0.1, 0.5, 0.05, 0.05, 0.3]):
    """Evaluate congruence between gene predictions from different sources.

    Gene predictions are assumed to be microbial protein-coding genes, marked
    as type `CDS`. Any reference protein alignments can be included as features
    of type `translated_nucleotide_match`. Compute the agreement between
    different sources of start/stop codon position, entire ORF position, and
    coverage from reference proteins.
    """
    starts = defaultdict(int)
    ends = defaultdict(int)
    rpstarts = defaultdict(int)
    rpends = defaultdict(int)
    coverage = defaultdict(int)

    for feature in tag.select.features(features, type='CDS', traverse=True):
        starts[feature.start] += 1
        ends[feature.end] += 1

    aligns = list(
        tag.select.features(features, type='translated_nucleotide_match')
    )
    for alignment in aligns:
        rpstarts[alignment.start] += 1
        rpends[alignment.end] += 1

    for feat in tag.select.features(features, type='CDS', traverse=True):
        bestcov = 0.0
        for alignment in aligns:
            bpoverlap = feat.range.overlap_extent(alignment.range)
            if bpoverlap == 0:
                continue
            bpmerged = len(feat.range.merge(alignment.range))
            recip_coverage = bpoverlap / bpmerged
            if bestcov is None or recip_coverage > bestcov:
                bestcov = recip_coverage
        coverage[feat] = bestcov

    for feature in tag.select.features(features, type='CDS', traverse=True):
        start_confirmed = starts[feature.start] - 1
        end_confirmed = ends[feature.end] - 1
        rpstart_confirmed = rpstarts[feature.start]
        rpend_confirmed = rpends[feature.end]
        prot_coverage = coverage[feature]

        feature.add_attribute('start_confirmed', start_confirmed)
        feature.add_attribute('end_confirmed', end_confirmed)
        feature.add_attribute('refrprot_start_confirmed', rpstart_confirmed)
        feature.add_attribute('refrprot_end_confirmed', rpend_confirmed)
        feature.add_attribute('protein_recip_coverage', prot_coverage)

        def coord_score(support):
            if support <= 0:
                return 0.0
            elif support == 1:
                return 0.8
            else:
                return 1.0

        def coord_score_refrprot(support):
            if support > 0:
                return 1.0
            else:
                return 0.0

        score_components = (
            coord_score(start_confirmed), coord_score(end_confirmed),
            coord_score_refrprot(rpstart_confirmed),
            coord_score_refrprot(rpend_confirmed), prot_coverage
        )
        score = sum([c * w for c, w in zip(score_components, weights)])
        feature.score = score


def eval_stream(locusstream, weights=[0.1, 0.5, 0.05, 0.05, 0.3]):
    """Feature stream for bacterial annotation evaluation."""
    for seqid, interval, locus in locusstream:
        eval_locus(locus, weights=weights)
        for feature in locus:
            yield feature
        yield tag.Directive('###')
