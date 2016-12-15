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


def main(args):
    assert args.cmd in mains
    mainmethod = mains[args.cmd]
    mainmethod(args)


if __name__ == '__main__':
    main(tag.cli.parser().parse_args())
