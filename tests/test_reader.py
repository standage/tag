#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import pytest
from aeneas.comment import Comment
from aeneas.directive import Directive
from aeneas.feature import Feature
from aeneas.reader import GFF3Reader
from aeneas.sequence import Sequence


def test_grape():
    """Sanity check."""
    with open('tests/testdata/grape-cpgat.gff3', 'r') as infile:
        reader = GFF3Reader(instream=infile)
        output = ''
        for record in reader:
            output += '%r\n' % record
        assert output == open('tests/testdata/grape-cpgat-sorted.gff3').read()


def test_pdom():
    """Pdom GFF3 with sequence."""
    reader = GFF3Reader(infilename='tests/testdata/pdom-withseq.gff3',
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
    infile = open('tests/testdata/pbar-withseq.gff3', 'r')
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
    infile = open('tests/testdata/pbar-withspaces.gff3', 'r')
    reader = GFF3Reader(instream=infile, assumesorted=True)
    records = list()
    for record in reader:
        records.append(record)
    assert len(records) == 6


def test_child_strand():
    """Lhum bad child strand."""
    reader = GFF3Reader(infilename='tests/testdata/lhum-cds-strand.gff3')
    with pytest.raises(AssertionError) as e:
        for record in reader:
            pass  # pragma: no cover
    assert 'child of feature LH19950-RA has a different strand' in str(e.value)


def test_parent_span():
    """Lhum child exceeds parent span."""
    reader = GFF3Reader(infilename='tests/testdata/lhum-mrna-span.gff3')
    with pytest.raises(AssertionError) as ae:
        for record in reader:
            pass  # pragma: no cover
    testmessage = 'child of feature LH19950 is not contained within its span'
    assert testmessage in str(ae.value)


def test_id_mismatch():
    """Lhum ID mismatch."""
    reader = GFF3Reader(infilename='tests/testdata/lhum-feat-dup.gff3')
    with pytest.raises(AssertionError) as ae:
        for record in reader:
            pass  # pragma: no cover
    testmessage = 'seqid mismatch for feature LH19950'
    assert testmessage in str(ae.value)
