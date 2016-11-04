#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import pytest
import sys
from .comment import Comment
from .directive import Directive
from .feature import Feature
from .sequence import Sequence


def parse_fasta(data):
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
                 assumesorted=False):
        assert (not instream) != (not infilename), ('provide either an '
                                                    'instream or an infile '
                                                    'name, not both')
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
    """Sanity check."""
    with open('testdata/grape-cpgat.gff3', 'r') as infile:
        reader = GFF3Reader(instream=infile)
        output = ''
        for record in reader:
            output += '%r\n' % record
        assert output == open('testdata/grape-cpgat-sorted.gff3').read()


def test_pdom():
    """Pdom GFF3 with sequence."""
    reader = GFF3Reader(infilename='testdata/pdom-withseq.gff3',
                        assumesorted=True)
    records = list()
    for record in reader:
        records.append(record)
    assert len(records) == 6

    record = records.pop()
    assert isinstance(record, Sequence) and \
        record.seqid == 'PdomSCFr1.2-0483'

    record = records.pop()
    assert isinstance(record, Feature) and \
        record.get_attribute('Name') == 'PdomGENEr1.2-04310'

    record = records.pop()
    assert isinstance(record, Feature) and \
        record.get_attribute('Name') == 'PdomGENEr1.2-04123'

    record = records.pop()
    assert isinstance(record, Comment) and \
        str(record) == 'Small scaffold from AllPathsLG assembly.'

    record = records.pop()
    assert isinstance(record, Directive) and \
        record.dirtype == 'sequence-region' and \
        record.seqid == 'PdomSCFr1.2-0483'

    record = records.pop()
    assert isinstance(record, Directive) and \
        record.dirtype == 'gff-version' and \
        record.version == '3'


def test_pbar():
    """Pbar GFF3 with sequence."""
    infile = open('testdata/pbar-withseq.gff3', 'r')
    reader = GFF3Reader(instream=infile, assumesorted=True)
    records = list()
    for record in reader:
        records.append(record)
    assert len(records) == 8

    record = records.pop()
    assert isinstance(record, Sequence) and \
        record.accession == 'NW_011929623.1'

    record = records.pop()
    assert isinstance(record, Sequence) and \
        record.seqid == 'gi|759031810|ref|NW_011929624.1|'

    record = records.pop()
    assert isinstance(record, Feature) and \
        record.get_attribute('Name') == 'LOC105423170'

    record = records.pop()
    assert isinstance(record, Feature) and \
        record.get_attribute('ID') == 'pseudogene1'

    record = records.pop()
    assert isinstance(record, Feature) and \
        record.get_attribute('Dbxref') == 'GeneID:105422795'

    record = records.pop()
    assert isinstance(record, Directive) and \
        record.dirtype == 'sequence-region' and \
        record.seqid == 'NW_011929624.1'

    record = records.pop()
    assert isinstance(record, Directive) and \
        record.dirtype == 'sequence-region' and \
        record.seqid == 'NW_011929623.1'

    record = records.pop()
    assert isinstance(record, Directive) and \
        record.dirtype == 'gff-version' and \
        record.version == '3'


def test_pbar_empty_lines():
    """Pbar GFF3 with empty lines."""
    infile = open('testdata/pbar-withspaces.gff3', 'r')
    reader = GFF3Reader(instream=infile, assumesorted=True)
    records = list()
    for record in reader:
        records.append(record)
    assert len(records) == 6


def test_child_strand():
    """Lhum bad child strand."""
    reader = GFF3Reader(infilename='testdata/lhum-cds-strand.gff3')
    with pytest.raises(AssertionError) as e:
        for record in reader:
            pass  # pragma: no cover
    assert 'child of feature LH19950-RA has a different strand' in str(e.value)


def test_parent_span():
    """Lhum child exceeds parent span."""
    reader = GFF3Reader(infilename='testdata/lhum-mrna-span.gff3')
    with pytest.raises(AssertionError) as ae:
        for record in reader:
            pass  # pragma: no cover
    testmessage = 'child of feature LH19950 is not contained within its span '
    assert testmessage in str(ae.value)


def test_id_mismatch():
    """Lhum ID mismatch."""
    reader = GFF3Reader(infilename='testdata/lhum-feat-dup.gff3')
    with pytest.raises(AssertionError) as ae:
        for record in reader:
            pass  # pragma: no cover
    testmessage = ('feature seq disagreement for ID="LH19950": scaffold2 vs '
                   'scaffold1')
    assert testmessage in str(ae.value)
