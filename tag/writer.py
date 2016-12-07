#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from collections import defaultdict
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import sys
import tag
from tag.feature import Feature


class GFF3Writer():
    """Writes sequence features and other GFF3 entries to a file."""

    def __init__(self, instream, outfile='-'):
        self._instream = instream
        self.outfilename = outfile
        self.outfile = None
        if outfile == '-':
            self.outfile == sys.stdout
        elif isinstance(outfile, StringIO):
            self.outfile = outfile
        else:
            self.outfile = tag.open(outfile, 'w')
        self.retainids = False
        self.feature_counts = defaultdict(int)

    def __del__(self):
        if self.outfilename != '-':
            self.outfile.close()

    def write(self, relax=False):
        for entry in self._instream:
            if isinstance(entry, Feature):
                for feature in entry:
                    if feature.num_children > 0:
                        self.feature_counts[feature.type] += 1
                        fid = '{}{}'.format(feature.type,
                                            self.feature_counts[feature.type])
                        feature.add_attribute('ID', fid)
                    else:
                        feature.drop_attribute('ID')
            print(repr(entry), file=self.outfile)
