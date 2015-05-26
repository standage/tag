#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the ISC license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import sys
from .comment import Comment
from .directive import Directive
from .feature import Feature


class GFF3Reader():
    """Loads sequence features and other GFF3 entries into memory."""

    def __init__(self, instream=sys.stdin, infilename=None,
                 assumesorted=False):
        assert not instream or not infilename, ('provide either an instream '
                                                'or an infile name, not both')
        self.instream = instream
        if infilename:
            self.infilename = infilename
            self.instream = open(infilename, 'r')
        self.assumesorted = assumesorted

    def __iter__(self):
        """Generator function returns GFF3 entries."""
        self._reset()

        for line in self.instream:
            line = line.rstrip()
            if line == '':
                continue
            elif line == '###':
                if self.assumesorted:
                    for obj in self._resolve_features():
                        yield obj
            elif line.startswith('#'):
                if line.startswith('##') and line[2] != '#':
                    record = Directive(line)
                else:
                    record = Comment(line)
                self.records.append(record)
            else:
                feature = Feature(line)

                parentid = feature.get_attribute('Parent')
                if parentid is None:
                    self.records.append(feature)
                else:
                    if parentid not in self.featsbyparent:
                        self.featsbyparent[parentid] = list()
                    self.featsbyparent[parentid].append(feature)

                featureid = feature.get_attribute('ID')
                if featureid is not None:
                    if featureid in self.featsbyid:
                        # Validate multi-features
                        other = self.featsbyid[featureid]
                        assert feature.type == other.type, \
                            'feature type disagreement for ID="%s": %s vs %s' \
                            % (featureid, feature.type, other.type)
                        assert feature.seqid == other.seqid, \
                            'feature seq disagreement for ID="%s": %s vs %s' \
                            % (featureid, feature.seqid, other.seqid)
                        other.add_sibling(feature)
                    else:
                        self.featsbyid[featureid] = feature

        for obj in self._resolve_features():
            yield obj

    def _resolve_features(self):
        """Resolve Parent/ID relationships and yield all top-level features."""

        for parentid in self.featsbyparent:
            parent = self.featsbyid[parentid]
            for child in self.featsbyparent[parentid]:
                parent.add_child(child)
                if child.is_multi:
                    for sibling in child.multi_rep.siblings:
                        if sibling != child:
                            parent.add_child(sibling)

        for record in sorted(self.records):
            yield record
        self._reset()

    def _reset(self):
        """Clear internal data structure."""
        self.records = list()
        self.featsbyid = dict()
        self.featsbyparent = dict()
        self.countsbytype = dict()


def test_grape():
    """[aeneas::GFF3Reader] Sanity check. More coming soon."""
    infile = open('testdata/grape-cpgat.gff3')
    reader = GFF3Reader(instream=infile)
    output = ''
    for record in reader:
        output += '%r\n' % record
    assert output == open('testdata/grape-cpgat-sorted.gff3').read(), output
