#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import tag


mains = {
    'gff3': tag.cli.gff3.main,
    'occ': tag.cli.occ.main,
}


def main(args=None):
    """
    Entry point for the tag CLI.

    Isolated as a method so that the CLI can be called by other Python code
    (e.g. for testing), in which case the arguments are passed to the function.
    If no arguments are passed to the function, parse them from the command
    line.
    """
    if args is None:
        args = tag.cli.parser().parse_args()

    assert args.cmd in mains
    mainmethod = mains[args.cmd]
    mainmethod(args)
