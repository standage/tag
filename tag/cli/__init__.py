#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import argparse
import sys
import tag
from . import gff3
from . import occ
from . import pmrna
from . import sum


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version',
                        version='tag v{}'.format(tag.__version__))
    parser.add_argument('-l', '--logfile', metavar='FILE', default=sys.stderr,
                        type=argparse.FileType('w'))
    subparsers = parser.add_subparsers(dest='cmd', metavar='cmd',
                                       help='gff3 | occ | pmrna | sum')
    tag.cli.gff3.subparser(subparsers)
    tag.cli.occ.subparser(subparsers)
    tag.cli.pmrna.subparser(subparsers)
    tag.cli.sum.subparser(subparsers)

    return parser
