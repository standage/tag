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
import intervaltree
import tag


def interval_set_span(intset):
    begin = min([x for x, y, z in intset])
    end = max([y for x, y, z in intset])
    cds = set([z for x, y, z in intset])
    return begin, end, cds

parser = argparse.ArgumentParser()
parser.add_argument('gff3', type=argparse.FileType('r'))
args = parser.parse_args()

coding_exons = dict()
reader = tag.reader.GFF3Reader(instream=args.gff3)
for record in reader:
    if not isinstance(record, tag.feature.Feature):
        continue

    for feature in record:
        if feature.type != 'mRNA':
            continue

        for subfeat in feature:
            if subfeat.type == 'CDS':
                if subfeat.seqid not in coding_exons:
                    coding_exons[subfeat.seqid] = intervaltree.IntervalTree()
                coding_exons[subfeat.seqid].addi(subfeat.start - 1,
                                                 subfeat.end,
                                                 subfeat)

total_occ = 0
ints_acct_for = dict()
for seqid in coding_exons:
    ints_acct_for[seqid] = intervaltree.IntervalTree()
    for cds_interval in coding_exons[seqid]:
        begin, end, cds = cds_interval
        if ints_acct_for[seqid][begin:end] != set():
            continue

        cds = set(cds)
        overlapping = coding_exons[seqid][begin:end]
        testbegin, testend, testcds = interval_set_span(overlapping)
        while set(cds) < testcds:
            begin, end, cds = testbegin, testend, testcds
            overlapping = coding_exons[seqid][begin:end]
            testbegin, testend, testcds = interval_set_span(overlapping)
        total_occ += end - begin
        ints_acct_for[seqid].addi(begin, end, cds)
print(total_occ)
