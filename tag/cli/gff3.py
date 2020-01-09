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


def subparser(subparsers):
    subparser = subparsers.add_parser('gff3')
    subparser.add_argument(
        '-n', '--no-sort', action='store_true', help='do not enforce sorting '
        'of output data'
    )
    subparser.add_argument(
        '-o', '--out', metavar='FILE', help='write output in GFF3 format to '
        'FILE; default is terminal (stdout)'
    )
    subparser.add_argument(
        '-r', '--relax', action='store_false', default=True, dest='strict',
        help='relax parsing stringency'
    )
    subparser.add_argument(
        '-i', '--retain-ids', action='store_true',
        help='retain feature IDs from input'
    )
    subparser.add_argument(
        '-s', '--sorted', action='store_true', help='assume the input data is '
        'sorted'
    )
    subparser.add_argument('gff3', help='input file in GFF3 format')


def main(args):
    reader = tag.reader.GFF3Reader(
        infilename=args.gff3, strict=args.strict, assumesorted=args.sorted,
        checkorder=not args.no_sort
    )
    writer = tag.writer.GFF3Writer(reader, args.out)
    writer.retainids = args.retain_ids
    writer.write()
