#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import tag
from tag import Range, Comment, Directive, Feature, Sequence, GFF3Reader
from tag.reader import DuplicatedRegionError, FeatureTypeDisagreementError
from tag.tests import data_file, data_stream


def test_grape():
    """Sanity check."""
    with open(data_file('grape-cpgat.gff3'), 'r') as infile:
        reader = GFF3Reader(instream=infile)
        output = ''
        for record in reader:
            output += repr(record) + '\n'
            if isinstance(record, Feature) and record.is_complex:
                output += '###\n'
        assert output == open(data_file('grape-cpgat-sorted.gff3')).read()


def test_pdom():
    """Pdom GFF3 with sequence."""
    reader = GFF3Reader(infilename=data_file('pdom-withseq.gff3'),
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
    infile = open(data_file('pbar-withseq.gff3'), 'r')
    reader = GFF3Reader(instream=infile, assumesorted=True)
    records = list()
    for record in reader:
        records.append(record)
    assert len(records) == 9

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
    assert isinstance(record, Comment)

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
    infile = open(data_file('pbar-withspaces.gff3'), 'r')
    reader = GFF3Reader(instream=infile, assumesorted=True)
    records = list()
    for record in reader:
        records.append(record)
    assert len(records) == 6


def test_child_strand():
    """Lhum bad child strand."""
    reader = GFF3Reader(infilename=data_file('lhum-cds-strand.gff3'))
    with pytest.raises(AssertionError) as e:
        for record in reader:
            pass  # pragma: no cover
    assert 'child of feature LH19950-RA has a different strand' in str(e.value)


def test_parent_span():
    """Lhum child exceeds parent span."""
    reader = GFF3Reader(infilename=data_file('lhum-mrna-span.gff3'))
    with pytest.raises(AssertionError) as ae:
        for record in reader:
            pass  # pragma: no cover
    testmessage = 'child of feature LH19950 is not contained within its span'
    assert testmessage in str(ae.value)


def test_id_mismatch():
    """Lhum ID mismatch."""
    reader = GFF3Reader(infilename=data_file('lhum-feat-dup.gff3'))
    with pytest.raises(AssertionError) as ae:
        for record in reader:
            pass  # pragma: no cover
    testmessage = 'seqid mismatch for feature LH19950'
    assert testmessage in str(ae.value)


def test_multi():
    """Multi-feature handling"""
    reader = GFF3Reader(infilename=data_file('amel-cdna-multi.gff3'))
    records = [r for r in reader]
    assert len(records) == 4
    assert isinstance(records[0], Directive)
    assert isinstance(records[1], Directive)
    assert isinstance(records[2], Feature)
    assert records[2].type == 'gene'
    assert isinstance(records[3], Feature)
    assert records[3].type == 'cDNA_match'
    assert repr(records[3]) == \
        open(data_file('amel-cdna-multi-out.gff3')).read().strip()


def test_no_seqreg():
    """Missing sequence region pragmas"""
    reader = GFF3Reader(infilename=data_file('otau-no-seqreg.gff3'))
    records = [r for r in reader]
    assert len(records) == 8
    assert isinstance(records[0], Directive)
    assert records[0].type == 'gff-version'
    assert isinstance(records[1], Directive)
    assert records[1].type == 'sequence-region'
    assert records[1].range == Range(9779, 19804)
    for r in records[2:]:
        assert isinstance(r, Feature)
        assert r.type == 'gene'


def test_feature_out_of_range():
    """Feature out of the sequence-region-specified range."""
    reader = GFF3Reader(infilename=data_file('vcar-out-of-bounds.gff3'))
    with pytest.raises(ValueError) as e:
        records = [r for r in reader]
    assert 'feature gene@NW_003307548.1[95396, 100541] out-of-bounds' in str(e)

    reader = GFF3Reader(infilename=data_file('vcar-out-of-bounds.gff3'),
                        assumesorted=True)
    with pytest.raises(ValueError) as e:
        records = [r for r in reader]
    assert 'feature gene@NW_003307548.1[95396, 100541] out-of-bounds' in str(e)


def test_seqreg_dup():
    """Duplicated sequence-region pragma."""
    infile = data_stream('vcar-seqreg-dup.gff3.gz')
    reader = GFF3Reader(instream=infile)
    with pytest.raises(DuplicatedRegionError) as dre:
        records = [r for r in reader]
    assert 'NW_003307554.1' in str(dre)


def test_seqreg_outoforder():
    """sequence-region pragma declared after features on that sequence."""
    infile = data_stream('grape-cpgat-seqreg-after.gff3.gz')
    reader = GFF3Reader(instream=infile)
    records = [r for r in reader]
    assert len(records) == 5

    infile = data_stream('grape-cpgat-seqreg-after.gff3.gz')
    reader = GFF3Reader(instream=infile, assumesorted=True)
    with pytest.raises(ValueError) as ve:
        records = [r for r in reader]


def test_gzip_input():
    reader = GFF3Reader(infilename=data_file('gzipdata.gff3.gz'))
    records = [r for r in reader]
    assert len(records) == 5
    genes = [r for r in records if hasattr(r, 'type') and r.type == 'gene']
    assert len(genes) == 3


def test_multi_feat_type_mismatch():
    reader = GFF3Reader(infilename=data_file('eden-mismatch.gff3'))
    with pytest.raises(FeatureTypeDisagreementError) as ftde:
        records = [r for r in reader]
    assert ' CDS vs exon' in str(ftde) or 'exon vs CDS' in str(ftde)


@pytest.mark.parametrize('infile,check,positions', [
    ('grape-cpgat-unsorted.gff3', True, [72, 10538, 22053]),
    ('grape-cpgat-unsorted.gff3', False, [10538, 72, 22053]),
    ('grape-cpgat-shuffled.gff3', True, [72, 10538, 22053]),
    ('grape-cpgat-shuffled.gff3', False, [72, 10538, 22053]),
])
def test_no_sort(infile, check, positions):
    reader = GFF3Reader(infilename=data_file(infile), checkorder=check)
    genes = tag.select.features(reader, type='gene')
    testpos = [gene.start + 1 for gene in genes]
    assert testpos == positions
