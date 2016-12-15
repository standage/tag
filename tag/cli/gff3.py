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
    subparser.add_argument('-o', '--out', metavar='FILE',
                           default='/dev/stdout', help='write output to the '
                           'specified file; default is terminal (stdout)')
    subparser.add_argument('-r', '--relax', action='store_false', default=True,
                           dest='strict', help='relax parsing stringency')
    subparser.add_argument('gff3', help='input file in GFF3 format')


def main(args):
    reader = tag.reader.GFF3Reader(infilename=args.gff3, strict=args.strict)
    writer = tag.writer.GFF3Writer(reader, args.out)
    writer.write()
