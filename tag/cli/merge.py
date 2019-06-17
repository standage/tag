#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import argparse
import tag
from tag import GFF3Reader, GFF3Writer


def subparser(subparsers):
    subparser = subparsers.add_parser('merge')
    subparser.add_argument(
        '-o', '--out', metavar='FILE', help='write output in GFF3 to FILE; '
        'default is terminal (stdout)'
    )
    subparser.add_argument(
        '-r', '--relax', action='store_false', default=True, dest='strict',
        help='relax parsing stringency'
    )
    subparser.add_argument(
        'gff3', nargs='+', help='input files in GFF3 format'
    )


def main(args):
    instreams = [
        GFF3Reader(infilename=fn, strict=args.strict, assumesorted=True)
        for fn in args.gff3
    ]
    merger = tag.select.merge(*instreams)
    writer = tag.writer.GFF3Writer(merger, args.out)
    writer.write()
