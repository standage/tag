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
except ImportError:  # pragma: no cover
    from io import StringIO
import sys
import tag
from tag import Directive, Feature, Sequence, GFF3Reader


class GFF3Writer():
    """
    Writes sequence features and other GFF3 entries to a file.

    The :code:`instream` is expected to be an iterable of sequence features and
    other related objects. Set :code:`outfile` to :code:`-` to write output
    to stdout.

    >>> # Sort and tidy GFF3 file in 3 lines!
    >>> reader = GFF3Reader(infilename=tag.tests.data_file('grape-cpgat.gff3'))
    >>> writer = GFF3Writer(instream=reader, outfile='/dev/null')
    >>> writer.retainids = True
    >>> writer.write()
    """

    def __init__(self, instream, outfile='-'):
        self._instream = instream
        self.outfilename = outfile
        self.outfile = None
        if isinstance(outfile, str):
            self.outfile = tag.open(outfile, 'w')
        else:
            self.outfile = outfile
        self.retainids = False
        self.complex_separators = True
        self.feature_counts = defaultdict(int)
        self._seq_written = False
        self._block_count = 0

    def _write_separator(self, blockitvl):
        if not blockitvl:
            return
        if self._block_count < blockitvl:
            return
        print('###', file=self.outfile)
        self._block_count = 0

    def __del__(self):
        if self.outfilename != '-' and not isinstance(self.outfile, StringIO):
            self.outfile.close()

    def write(self, blockitvl=0):
        """Pull features from the instream and write them to the output.

        By default, separator tags are added at the end of complex features. To
        intersperse separators throughout blocks of simple features, specify a
        desired block size with `blockitvl`.
        """
        print(repr(Directive('##gff-version 3')), file=self.outfile)
        for entry in self._instream:
            if isinstance(entry, Directive):
                if entry.type == 'gff-version':
                    pass
                else:
                    print(repr(entry), file=self.outfile)
                continue
            if isinstance(entry, Feature):
                for feature in entry:
                    if self.retainids:
                        continue
                    if feature.num_children > 0 or feature.is_multi:
                        if feature.is_multi and feature != feature.multi_rep:
                            continue
                        self.feature_counts[feature.type] += 1
                        fid = '{}{}'.format(feature.type,
                                            self.feature_counts[feature.type])
                        feature.add_attribute('ID', fid)
                    else:
                        feature.drop_attribute('ID')
            if isinstance(entry, Sequence) and not self._seq_written:
                print('##FASTA', file=self.outfile)
                self._seq_written = True
            print(repr(entry), file=self.outfile)
            if isinstance(entry, Feature):
                if entry.is_complex:
                    self._block_count = 0
                    if self.complex_separators:
                        print('###', file=self.outfile)
                else:
                    self._block_count += 1
            self._write_separator(blockitvl)
