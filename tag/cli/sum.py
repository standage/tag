#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2017 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import argparse
from collections import defaultdict
import tag


def subparser(subparsers):
    desc = 'Briefly summarize a GFF3 file'
    subparser = subparsers.add_parser('sum', description=desc)
    subparser.add_argument('gff3', help='input file')


def main(args):
    annot = tag.index.Index()
    annot.consume_file(args.gff3)
    annot.yield_inferred = False

    seqs = list(annot.seqids)
    size = 0
    for seqid in seqs:
        start, end = annot.extent(seqid)
        size += end - start

    counts = defaultdict(int)
    maxlens = defaultdict(int)
    for feature in tag.select.features(annot):
        for subfeat in feature:
            counts[subfeat.type] += 1
            maxlens[subfeat.type] = max(len(subfeat), maxlens[subfeat.type])

    sumstr = 'Summary for file "{}":\n'.format(args.gff3)
    sumstr += '    - {} annotated sequences'.format(len(seqs))
    sumstr += ' for a total length of {} bp\n'.format(size)
    sumstr += '    - {} annotated features'.format(sum(counts.values()))
    sumstr += ' (or feature entries)\n'
    for ft in sorted(counts):
        sumstr += '        - {} entries of type {}'.format(counts[ft], ft)
        sumstr += ', maximum length: {} bp\n'.format(maxlens[ft])

    print(sumstr.strip())
