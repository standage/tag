#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import sys
import tag
from tag.range import Range
from tag.comment import Comment
from tag.directive import Directive
from tag.feature import Feature
from tag.sequence import Sequence


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
    """
    Loads sequence features and other GFF3 entries into memory.

    Features and other objects are obtained by iterating over the GFF3Reader
    object. See the :code:`tag.select` module for available filtering
    functions.

    >>> reader = GFF3Reader(infilename='tests/testdata/pbar-withseq.gff3')
    >>> for entry in reader:
    ...     print(type(entry))
    <class 'tag.directive.Directive'>
    <class 'tag.directive.Directive'>
    <class 'tag.directive.Directive'>
    <class 'tag.comment.Comment'>
    <class 'tag.feature.Feature'>
    <class 'tag.feature.Feature'>
    <class 'tag.feature.Feature'>
    <class 'tag.sequence.Sequence'>
    <class 'tag.sequence.Sequence'>
    >>> reader = GFF3Reader(infilename='tests/testdata/pbar-withseq.gff3')
    >>> for feature in tag.select.features(reader, type='gene'):
    ...     print(feature.slug)
    gene@NW_011929623.1[4557, 5749]
    gene@NW_011929624.1[3725, 4229]

    By default, the GFF3Reader assumes the input is not in sorted order and
    loads the entire annotation into memory to ensure sorting and resolution
    of ID/Parent relationships. If the input is sorted according to genomic
    position and independent features are separated by :code:`###` directives,
    memory consumption will be drastically reduced by setting the
    :code:`assumesorted` attribute to True.

    The :code:`strict` attribute enforces some additional sanity checks, which
    in some exceptional cases may need to be relaxed.
    """

    def __init__(self, instream=None, infilename=None,
                 assumesorted=False, strict=True):
        assert (not instream) != (not infilename), ('provide either an '
                                                    'instream or an infile '
                                                    'name, not both')
        self.instream = instream
        if infilename:
            self.infilename = infilename
            self.instream = tag.open(infilename, 'r')
        self.assumesorted = assumesorted
        self.strict = strict
        self.declared_regions = dict()
        self.inferred_regions = dict()
        self._counter = 0
        self._prevrecord = None

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
                        if self._prevrecord and self._prevrecord > obj:
                            msg = 'sorting error: '
                            msg += '{} > {}'.format(self._prevrecord, obj)
                            raise ValueError(msg)
                        self._prevrecord = obj
                        self._counter += 1
                        yield obj
            elif line.startswith('#'):
                if line == '##FASTA':
                    for sequence in parse_fasta(self.instream):
                        self.records.append(sequence)
                    break
                elif line.startswith('##') and line[2] != '#':
                    record = Directive(line)
                    if record.type == 'sequence-region':
                        if record.seqid in self.declared_regions:
                            m = ('sequence {} declared in multiple ##sequence'
                                 '-region entries'.format(record.seqid))
                            raise ValueError(m)
                        self.declared_regions[record.seqid] = record
                else:
                    record = Comment(line)
                self.records.append(record)
            else:
                feature = Feature(line)

                if feature.seqid in self.declared_regions:
                    seqregion = self.declared_regions[feature.seqid]
                    if not feature._range.within(seqregion.range):
                        msg = 'feature {} out-of-bounds'.format(feature.slug)
                        raise ValueError(msg)

                if feature.seqid not in self.inferred_regions:
                    self.inferred_regions[feature.seqid] = Range(feature.start,
                                                                 feature.end)
                rstr = self.inferred_regions[feature.seqid].start
                rend = self.inferred_regions[feature.seqid].end
                fstr = feature.start
                fend = feature.end
                self.inferred_regions[feature.seqid].start = min(rstr, fstr)
                self.inferred_regions[feature.seqid].end = max(rend, fend)

                parentid = feature.get_attribute('Parent')
                if parentid is None:
                    # All components of a multi-feature are added here...
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
            if isinstance(obj, Feature):
                if obj.is_multi and obj.multi_rep != obj:
                    # ...so we must filter multi-features here and only yield
                    # the multi-feature representative
                    continue

            if self._counter == 0:
                isv = isinstance(obj, Directive) and obj.type == 'gff-version'
                if not isv:
                    self._prevrecord = Directive('##gff-version 3')
                    self._counter += 1
                    yield self._prevrecord

            self._prevrecord = obj
            self._counter += 1
            yield obj

    def _resolve_features(self):
        """Resolve Parent/ID relationships and yield all top-level features."""

        for parentid in self.featsbyparent:
            parent = self.featsbyid[parentid]
            for child in self.featsbyparent[parentid]:
                parent.add_child(child, rangecheck=self.strict)

        if not self.assumesorted:
            for seqid in self.inferred_regions:
                if seqid not in self.declared_regions:
                    seqrange = self.inferred_regions[seqid]
                    srstring = '##sequence-region {:s} {:d} {:d}'.format(
                        seqid, seqrange.start + 1, seqrange.end
                    )
                    seqregion = Directive(srstring)
                    self.records.append(seqregion)

        for record in sorted(self.records):
            yield record
        self._reset()

    def _reset(self):
        """Clear internal data structure."""
        self.records = list()
        self.featsbyid = dict()
        self.featsbyparent = dict()
        self.countsbytype = dict()
