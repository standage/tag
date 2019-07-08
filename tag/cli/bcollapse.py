#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# -----------------------------------------------------------------------------
# Copyright (C) 2019 Battelle National Biodefense Institute.
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import argparse
import tag
from tag import GFF3Reader, GFF3Writer


def subparser(subparsers):
    subparser = subparsers.add_parser('bcollapse')
    subparser.add_argument(
        '-o', '--out', metavar='FILE', help='write output in GFF3 to FILE; '
        'default is terminal (stdout)'
    )
    subparser.add_argument(
        '-d', '--delta', metavar='Δ', type=int, default=0, help='extend locus '
        'interval by Δ bp in each direction; by default, Δ=0'
    )
    subparser.add_argument(
        '-n', '--min-bp', metavar='N', type=int, default=25, help='only group '
        'features together in the same locus if they overlap by at least N '
        'bp; by default N=25'
    )
    subparser.add_argument(
        '-p', '--min-perc', metavar='P', type=float, default=0.25,
        help='only group features together in the same locus if they overlap '
        'by a fraction of at least P of their overall length; by default '
        'P=0.25 (25%%)'
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
    locusstream = tag.locus.loci(
        *instreams, minbp=args.min_bp, minperc=args.min_perc
    )
    collapsestream = tag.bae.collapse_stream(locusstream)
    writer = tag.writer.GFF3Writer(collapsestream, args.out)
    writer.complex_separators = False
    writer.write()
