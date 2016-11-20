#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import sys
from .comment import Comment
from .directive import Directive
from .feature import Feature
from .sequence import Sequence


def parse_fasta(data):  # pragma: no cover
    """
    Load sequences in Fasta format.

    This generator function yields a Sequence object for each sequence record
    in a GFF3 file. Implementation stolen shamelessly from
    http://stackoverflow.com/a/7655072/459780.
    """
    name, seq = None, []
    for line in data:
        line = line.rstrip()
        if line.startswith('>'):
            if name:
                yield Sequence(name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield Sequence(name, ''.join(seq))


class GFF3Reader():
    """Loads sequence features and other GFF3 entries into memory."""

    def __init__(self, instream=None, infilename=None,
                 assumesorted=False, strict=True):
        assert (not instream) != (not infilename), ('provide either an '
                                                    'instream or an infile '
                                                    'name, not both')
        self.instream = instream
        if infilename:
            self.infilename = infilename
            self.instream = open(infilename, 'r')
        self.assumesorted = assumesorted
        self.strict = strict

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
                if line == '##FASTA':
                    for sequence in parse_fasta(self.instream):
                        self.records.append(sequence)
                    break
                elif line.startswith('##') and line[2] != '#':
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
                parent.add_child(child, rangecheck=self.strict)

        for record in sorted(self.records):
            yield record
        self._reset()

    def _reset(self):
        """Clear internal data structure."""
        self.records = list()
        self.featsbyid = dict()
        self.featsbyparent = dict()
        self.countsbytype = dict()
