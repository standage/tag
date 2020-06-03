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
from . import bae
from . import bcollapse
from . import gff3
from . import locuspocus
from . import merge
from . import occ
from . import pep2nuc
from . import pmrna
from . import sum

subparser_funcs = {
    'bae': bae.subparser,
    'bcollapse': bcollapse.subparser,
    'gff3': gff3.subparser,
    'locuspocus': locuspocus.subparser,
    'merge': merge.subparser,
    'occ': occ.subparser,
    'pep2nuc': pep2nuc.subparser,
    'pmrna': pmrna.subparser,
    'sum': sum.subparser,
}

mains = {
    'bae': bae.main,
    'bcollapse': bcollapse.main,
    'gff3': gff3.main,
    'locuspocus': locuspocus.main,
    'merge': merge.main,
    'occ': occ.main,
    'pmrna': pmrna.main,
    'pep2nuc': pep2nuc.main,
    'sum': sum.main,
}


def parser():
    subcmdstr = ', '.join(mains)
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--version', action='version',
                        version='tag v{}'.format(tag.__version__))
    subparsers = parser.add_subparsers(dest='cmd', metavar='cmd',
                                       help=subcmdstr)
    for subcmd, parserfunc in subparser_funcs.items():
        parserfunc(subparsers)
    return parser
