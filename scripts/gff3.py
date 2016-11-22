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
import tag

parser = argparse.ArgumentParser()
parser.add_argument('gff3', type=argparse.FileType('r'))
args = parser.parse_args()

reader = tag.reader.GFF3Reader(args.gff3)
for feature in reader:
    print(repr(feature))
