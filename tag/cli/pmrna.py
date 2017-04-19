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


def subparser(subparsers):
    subparser = subparsers.add_parser('pmrna')
    subparser.add_argument('-r', '--relax', action='store_false', default=True,
                           dest='strict', help='relax parsing stringency')
    subparser.add_argument('gff3', help='input file')


def main(args):
    reader = tag.GFF3Reader(infilename=args.gff3, strict=args.strict)
    writer = tag.GFF3Writer(tag.transcript.primary_mrna(reader))
    writer.write()
