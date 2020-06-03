#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
from intervaltree import IntervalTree
import tag
import sys


class Index(defaultdict):
    """
    In-memory index for efficient interval-based queries for genome features.

    Implemented as a dictionary, with sequence IDs as keys to interval trees
    of features for the corresponding scaffold or contig sequence.

    >>> index = tag.index.Index()
    >>> index.consume_file(tag.tests.data_file('pcan-123.gff3.gz'))
    >>> list(index.seqids)
    ['scaffold_123', 'scaffold_124', 'scaffold_125']
    >>> index.extent('scaffold_124')
    (6020, 444332)
    >>> for feature in index.query('scaffold_123', 5000, 6000):
    ...     print(feature.slug)
    gene@scaffold_123[5583, 5894]
    >>> for feature in index.query('scaffold_125', 19000, 87000, strict=False):
    ...     print(feature.slug)
    gene@scaffold_125[18994, 19335]
    gene@scaffold_125[57450, 57680]
    gene@scaffold_125[86995, 87151]
    >>> reader = tag.reader.GFF3Reader(
    ...     tag.tests.data_stream('grape-cpgat.gff3')
    ... )
    >>> index = tag.index.Index()
    >>> index.consume(reader)
    >>> index.yield_inferred = False
    >>> for entry in index:
    ...     print(entry.slug)
    sequence chr8[1, 100000]
    gene@chr8[72, 5081]
    gene@chr8[10538, 11678]
    gene@chr8[22053, 23448]
    """

    def __init__(self):
        defaultdict.__init__(self, IntervalTree)
        self.declared_regions = dict()
        self.inferred_regions = dict()
        self.yield_inferred = True

    def consume_file(self, infile):
        """Load the specified GFF3 file into memory."""
        reader = tag.reader.GFF3Reader(infilename=infile)
        self.consume(reader)

    def consume_seqreg(self, seqreg):
        """Load a :code:`##sequence-region` directive into memory."""
        if not isinstance(seqreg, tag.directive.Directive) or \
                seqreg.type != 'sequence-region':
            raise ValueError('expected ##sequence-region directive')
        if seqreg.seqid in self.declared_regions:
            msg = 'duplicate sequence region "{}"'.format(seqreg.seqid)
            raise ValueError(msg)
        self.declared_regions[seqreg.seqid] = seqreg.range.copy()

    def consume_feature(self, feature):
        """Load a :code:`Feature` object into memory."""
        if not isinstance(feature, tag.feature.Feature):
            raise ValueError('expected Feature object')
        self[feature.seqid][feature.start:feature.end] = feature
        if feature.seqid not in self.inferred_regions:
            self.inferred_regions[feature.seqid] = feature._range.copy()
        newrange = self.inferred_regions[feature.seqid].merge(feature._range)
        self.inferred_regions[feature.seqid].start = newrange.start
        self.inferred_regions[feature.seqid].end = newrange.end

    def consume(self, entrystream):
        """
        Load a stream of entries into memory.

        Only Feature objects and sequence-region directives are loaded, all
        other entries are discarded.
        """
        for entry in entrystream:
            if isinstance(entry, tag.directive.Directive) and \
                    entry.type == 'sequence-region':
                self.consume_seqreg(entry)
            elif isinstance(entry, tag.feature.Feature):
                self.consume_feature(entry)

    def __iter__(self):
        regions = self.inferred_regions
        if not self.yield_inferred:
            regions = self.declared_regions

        for seqid in sorted(list(self.keys())):
            sr = regions[seqid]
            template = '##sequence-region {} {} {}'
            data = template.format(seqid, sr.start + 1, sr.end)
            yield tag.directive.Directive(data)

        for seqid in sorted(list(self.keys())):
            for interval in sorted(self[seqid]):
                yield interval.data

    def query(self, seqid, start, end=None, strict=True):
        """
        Query the index for features in the specified range.

        :param seqid: ID of the sequence to query
        :param start: start of the query interval
        :param end: end of the query interval
        :param strict: indicates whether query is strict containment or overlap
                       (:code:`True` and :code:`False`, respectively)
        """
        if end and strict:
            query = self[seqid].envelop
            args = (start, end)
        elif end and not strict:
            query = self[seqid].overlap
            args = (start, end)
        else:
            query = self[seqid].at
            args = [start]
        return sorted([intvl.data for intvl in query(*args)])

    @property
    def seqids(self):
        for seqid in sorted(self.keys()):
            yield seqid

    def extent(self, seqid):
        regions = self.inferred_regions
        if not self.yield_inferred:
            regions = self.declared_regions
        sr = regions[seqid]
        return sr.start, sr.end


class NamedIndex(object):
    """
    In-memory index for retrieving genome features by identifier.

    >>> index = tag.index.NamedIndex()
    >>> index.consume_file(tag.tests.data_file('pcan-123.gff3.gz'))
    >>> gene = index['PCAN011a001813']
    >>> print(gene.slug)
    gene@scaffold_123[200029, 201298]
    >>> exon = index['PCAN011a001859T1.exon1']
    >>> print(exon.slug)
    exon@scaffold_125[57450, 57680]
    """

    def __init__(self):
        self.data = dict()

    def consume_file(self, infile, attribute='ID'):
        reader = tag.reader.GFF3Reader(infilename=infile)
        self.consume(reader, attribute=attribute)

    def consume(self, entrystream, attribute='ID'):
        for feature in tag.select.features(entrystream):
            for subfeature in feature:
                name = subfeature.get_attribute(attribute)
                if name is None:
                    continue
                if subfeature.is_multi and subfeature.multi_rep != subfeature:
                    continue
                self.data[name] = subfeature

    def __getitem__(self, name):
        if name not in self.data:
            raise IndexError(name)
        return self.data[name]

    def __contains__(self, name):
        return name in self.data

    @property
    def names(self):
        for name in sorted(self.data.keys()):
            yield name
