#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import argparse
from collections import defaultdict
from intervaltree import IntervalTree
import tag


def interval_set_span(intset):
    begin = min([x for x, y, z in intset])
    end = max([y for x, y, z in intset])
    feats = set([z for x, y, z in intset])
    return begin, end, feats


def subparser(subparsers):
    subparser = subparsers.add_parser('occ')
    subparser.add_argument('-r', '--relax', action='store_false', default=True,
                           dest='strict', help='relax parsing stringency')
    subparser.add_argument('gff3', help='input file')
    subparser.add_argument('type', help='feature type')


def main(args):
    features = defaultdict(IntervalTree)
    reader = tag.reader.GFF3Reader(infilename=args.gff3, strict=args.strict)
    for feature in tag.select.features(reader, type=args.type, traverse=True):
        if feature.is_pseudo:
            for sub in feature:
                features[sub.seqid].addi(sub.start, sub.end, sub)
        else:
            features[feature.seqid].addi(feature.start, feature.end, feature)

    total_occ = 0
    ints_acct_for = defaultdict(IntervalTree)
    for seqid in features:
        for interval in features[seqid]:
            begin, end, feat = interval
            if ints_acct_for[seqid][begin:end] != set():
                continue

            feats = set([feat])
            overlapping = features[seqid][begin:end]
            testbegin, testend, testfeats = interval_set_span(overlapping)
            while set(feats) < testfeats:
                begin, end, feats = testbegin, testend, testfeats
                overlapping = features[seqid][begin:end]
                testbegin, testend, testfeats = interval_set_span(overlapping)
            total_occ += end - begin
            ints_acct_for[seqid].addi(begin, end, feats)
    print(total_occ)
