#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2020 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import argparse
import sys
import tag


def pep2nuc(genomestream, protstream, attr='ID', keepattr=None):
    index = tag.index.NamedIndex()
    index.consume(tag.GFF3Reader(genomestream), attribute=attr)
    reader = tag.GFF3Reader(protstream)
    for feature in tag.select.features(reader):
        if feature.seqid not in index:
            print(
                '[tag::pep2nuc] WARNING:',
                'protein identifier "{}" not defined'.format(feature.seqid),
                file=sys.stderr
            )
            continue
        gfeat = index[feature.seqid]
        newstart = gfeat.start + (feature.start * 3)
        newend = gfeat.start + (feature.end * 3)
        protid = feature.seqid
        feature.seqid = gfeat.seqid
        feature.set_coord(newstart, newend)
        if keepattr:
            feature.add_attribute(keepattr, protid)
        yield feature


def subparser(subparsers):
    subparser = subparsers.add_parser('pep2nuc')
    subparser.add_argument(
        '-o', '--out', metavar='FILE', default='-', help='file to which '
        'output will be written; by default, output is written to the '
        'terminal (stdout)'
    )
    subparser.add_argument(
        '-a', '--attr', metavar='ATTR', default='ID', help='CDS/protein '
        'attribute in the genome GFF3 file that corresponds to the protein '
        'identifier (column 1) of the protein GFF3 file; "ID" by default'
    )
    subparser.add_argument(
        '-k', '--keep-prot', metavar='ATTR', help='keep the original protein '
        'ID and write it to the specified attribute in the output'
    )
    subparser.add_argument(
        'genome', help='GFF3 file with CDS or protein features defined on a '
        'genomic (contig, scaffold, or chromosome) coordinate system'
    )
    subparser.add_argument(
        'protein', help='GFF3 file with features defined on a protein '
        'coordinate system'
    )


def main(args):
    with tag.open(args.genome, 'r') as genomestream:
        with tag.open(args.protein, 'r') as protstream:
            transformer = pep2nuc(
                genomestream, protstream, attr=args.attr,
                keepattr=args.keep_prot
            )
            writer = tag.GFF3Writer(transformer, outfile=args.out)
            writer.retainids = True
            writer.write()
