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
from aeneas.range import Range
from aeneas.sequence import Sequence


def eden():
    """
    Data fixture for unit tests.

    This function is not a good example of how to parse arbitrary GFF3 files,
    it's quick and dirty for this particular case.
    """
    edenfile = open('tests/testdata/eden.gff3', 'r')
    edengff3 = edenfile.readlines()
    gene = None
    mRNAs = list()
    cdsreps = dict()
    for line in edengff3:
        if line.startswith('#'):
            continue
        feature = Feature(line.rstrip())
        if feature.type == 'gene':
            gene = feature
        elif feature.type == 'TF_binding_site':
            gene.add_child(feature)
        elif feature.type == 'mRNA':
            mRNAs.append(feature)
            gene.add_child(feature)
        else:  # pragma: no cover
            for parent in feature.get_attribute('Parent', as_list=True):
                for mRNA in mRNAs:
                    if mRNA.get_attribute('ID') == parent:
                        mRNA.add_child(feature)
                        break
            if feature.type == 'CDS':
                cdsid = feature.get_attribute('ID')
                if cdsid not in cdsreps:
                    cdsreps[cdsid] = feature
                else:
                    cdsreps[cdsid].add_sibling(feature)
    assert len(gene.children) == 4
    return gene


def test_basic():
    """Test basic constructor."""
    gff3 = ['chr', 'vim', 'gene', '1000', '2000', '.', '+', '.', '.']
    f1 = Feature('\t'.join(gff3))
    assert '%s' % f1 == '\t'.join(gff3)
    assert '%r' % f1 == '\t'.join(gff3)
    assert len(f1) == 1001

    gff3[8] = 'ID=gene1;Name=EDEN'
    f2 = Feature('\t'.join(gff3))
    assert '%s' % f2 == '\t'.join(gff3)

    gff3[5] = '0.28466'
    f3 = Feature('\t'.join(gff3))
    assert abs(f3.score - 0.28466) <= 0.00001
    gff3[5] = '0.285'
    assert '%s' % f3 == '\t'.join(gff3)


def test_seqid():
    """Test seqid handling."""
    gff3 = ['chr', 'vim', 'gene', '1000', '2000', '.', '+', '.', 'ID=g1']
    f1 = Feature('\t'.join(gff3))
    assert f1.seqid == 'chr'

    f1.seqid = 'scaffold_28'
    assert f1.seqid == 'scaffold_28'
    gff3[0] = 'scaffold_28'
    assert '%s' % f1 == '\t'.join(gff3)

    gff3 = ['scaffold_28', 'vim', 'mRNA', '1000', '2000', '.', '+', '.',
            'ID=t1;Parent=g1']
    f2 = Feature('\t'.join(gff3))
    gff3 = ['scaffold_28', 'vim', 'exon', '1000', '2000', '.', '+', '.',
            'Parent=t1']
    f3 = Feature('\t'.join(gff3))
    f2.add_child(f3)
    f1.add_child(f2)
    for feature in f1:
        assert feature.seqid == 'scaffold_28'
    f1.seqid = 'ctg123'
    for feature in f1:
        assert feature.seqid == 'ctg123'


def test_source():
    """Test source handling."""
    gff3 = ['chr', 'vim', 'gene', '1000', '2000', '.', '+', '.', 'ID=g1']
    f1 = Feature('\t'.join(gff3))
    assert f1.source == 'vim'

    f1.source = 'nano'
    assert f1.source == 'nano'
    gff3[1] = 'nano'
    assert '%s' % f1 == '\t'.join(gff3)

    gff3 = ['chr', 'emacs', 'mRNA', '1000', '2000', '.', '+', '.',
            'ID=t1;Parent=g1']
    f2 = Feature('\t'.join(gff3))
    gff3 = ['chr', 'nano', 'exon', '1000', '2000', '.', '+', '.',
            'Parent=t1']
    f3 = Feature('\t'.join(gff3))
    f2.add_child(f3)
    f1.add_child(f2)
    for feature in f1:
        if feature.type == 'mRNA':
            assert feature.source == 'emacs'
        else:
            assert feature.source == 'nano'
    f1.source = 'gedit'
    for feature in f1:
        if feature.type == 'mRNA':
            assert feature.source == 'emacs'
        else:
            assert feature.source == 'gedit', '%s=%s' % (feature.type,
                                                         feature.source)


def test_type():
    """Test type handling."""
    gff3 = ['chr', 'vim', 'mRNA', '1000', '2000', '.', '+', '.', 'ID=mRNA1']
    f1 = Feature('\t'.join(gff3))
    assert f1.type == 'mRNA'

    gff3[2] = 'messenger RNA'
    f2 = Feature('\t'.join(gff3))
    assert f2.type == 'messenger RNA'


def test_region():
    """Test coordinate handling."""
    gff3 = ['chr', 'vim', 'gene', '1000', '2000', '.', '+', '.', 'ID=g1']
    f1 = Feature('\t'.join(gff3))
    assert f1._region == Range(999, 2000)
    assert f1.start == 999 and f1.end == 2000

    gff3 = ['contig5', 'vim', 'mRNA', '500', '2500', '.', '+', '.',
            'ID=t1;Parent=g1']
    f2 = Feature('\t'.join(gff3))
    with pytest.raises(AssertionError) as ae:
        f1.add_child(f2, rangecheck=True)
    assert 'seqid mismatch for feature g1' in str(ae)

    f2.seqid = 'chr'
    with pytest.raises(AssertionError) as ae:
        f1.add_child(f2, rangecheck=True)
    assert 'is not contained within its span' in str(ae)

    f1.add_child(f2)
    assert len(f1.children) == 1
    assert f2._region == Range(499, 2500)

    f2.set_coord(999, 2000)
    assert f2.start == 999 and f2.end == 2000

    f1.transform(100000)
    for feature in f1:
        assert feature._region == Range(100999, 102000)

    f1.transform(100000, newseqid='scf89')
    for feature in f1:
        assert feature.seqid == 'scf89'
        assert feature._region == Range(200999, 202000)


def test_score():
    """Test score handling."""
    gff3 = ['chr', 'vim', 'EST_match', '57229', '57404', '.', '+', '.', '.']
    f1 = Feature('\t'.join(gff3))
    assert f1.score is None

    gff3[5] = '0.97'
    f2 = Feature('\t'.join(gff3))
    assert abs(f2.score - 0.97) <= 0.00001
    gff3[5] = '0.970'
    assert str(f2) == '\t'.join(gff3)

    gff3[5] = '-1.8332'
    f3 = Feature('\t'.join(gff3))
    assert abs(f3.score + 1.8332) <= 0.00001
    gff3[5] = '-1.833'
    assert str(f3) == '\t'.join(gff3)


def test_strand():
    """Test strand handling."""
    gff3 = ['chr', 'vim', 'EST_match', '57229', '57404', '.', '.', '.', '.']
    f1 = Feature('\t'.join(gff3))
    assert f1.strand == '.'
    assert str(f1) == '\t'.join(gff3)

    gff3[6] = '-'
    f2 = Feature('\t'.join(gff3))
    assert f2.strand == '-'
    assert str(f2) == '\t'.join(gff3)

    gff3[2:5] = 'match_part', '57229', '57298'
    f3 = Feature('\t'.join(gff3))
    assert f2.strand == '-'

    gff3[3:5] = '57377', '57404'
    gff3[6] = '+'
    f4 = Feature('\t'.join(gff3))
    assert f4.strand == '+'

    f2.add_child(f3)
    with pytest.raises(AssertionError) as ae:
        f2.add_child(f4, rangecheck=True)
    assert 'has a different strand' in str(ae)
    assert len(f2.children) == 1
    f4._strand = '-'
    f2.add_child(f4)
    assert len(f2.children) == 2

    gff3 = ['chr', 'vim', 'EST_match', '57229', '57404', '.', '~', '.', '.']
    with pytest.raises(AssertionError) as ae:
        f5 = Feature('\t'.join(gff3))
    assert 'invalid strand' in str(ae)


def test_phase():
    """Test phase handling."""
    gff3 = ['chr', 'vim', 'mRNA', '1001', '1420', '.', '+', '.', 'ID=t1']
    m1 = Feature('\t'.join(gff3))
    assert m1.phase is None

    gff3 = ['chr', 'vim', 'CDS', '1001', '1100', '.', '+', '0', 'ID=CDS1']
    c1 = Feature('\t'.join(gff3))
    assert c1.phase == 0

    gff3[3:5] = '1201', '1236'
    gff3[7] = '2'
    c2 = Feature('\t'.join(gff3))
    c1.add_sibling(c2)
    assert c2.phase == 2

    gff3[3:5] = '1301', '1364'
    c3 = Feature('\t'.join(gff3))
    c1.add_sibling(c3)
    assert c3.phase == 2

    gff3[3:5] = '1401', '1420'
    gff3[7] = '1'
    c4 = Feature('\t'.join(gff3))
    c1.add_sibling(c4)
    assert c4.phase == 1


def test_attributes():
    """Test attribute handling."""
    gene = eden()
    assert gene.get_attribute('ID') == 'gene00001'
    assert gene.get_attribute('ID', as_list=True) == 'gene00001'
    assert gene.get_attribute('Name') == 'EDEN'
    assert gene.get_attribute('Name', as_list=True) == ['EDEN']
    assert gene.get_attribute_keys() == ['ID', 'Name']

    gene.add_attribute('Name', 'Gandalf')
    assert gene.get_attribute('Name') == 'Gandalf'
    gene.add_attribute('Name', 'Aragorn', append=True)
    assert gene.get_attribute('Name') == ['Aragorn', 'Gandalf']
    assert gene.get_attribute('Name', as_string=True) == 'Aragorn,Gandalf'

    gene.add_attribute('ID', 'g1')
    for child in gene.children:
        assert child.get_attribute('Parent') == 'g1'

    for feature in gene:
        if feature.type == 'exon' and feature._region == Range(4999, 5500):
            assert feature.get_attribute('Parent') == [
                'mRNA00001',
                'mRNA00002',
                'mRNA00003']
            assert feature.get_attribute('Parent', as_string=True) == \
                'mRNA00001,mRNA00002,mRNA00003'

    for feature in gene:
        if feature.get_attribute('ID') == 'mRNA00003':
            feature.add_attribute('ID', 'mRNA3')
        # elif feature.get_attribute('ID') in ['cds00003', 'cds00004']:
        #     assert feature.get_attribute('Parent') == 'mRNA3'

    assert repr(gene) == \
        open('tests/testdata/eden-mod.gff3', 'r').read().rstrip()

    gene.add_attribute('Note', 'I need to test access of other attributes')
    assert gene.attributes == ('ID=g1;Name=Aragorn,Gandalf;Note='
                               'I need to test access of other attributes')

    gff3 = ['chr', 'vim', 'mRNA', '1001', '1420', '.', '+', '.', '.']
    m1 = Feature('\t'.join(gff3))
    assert m1.phase is None
    m1.add_attribute('ID', 't1')
    assert m1.get_attribute('ID') == 't1'

    gff3 = ['chr', 'vim', 'CDS', '1001', '1420', '.', '+', '.', '.']
    c1 = Feature('\t'.join(gff3))
    m1.add_child(c1)
    m1.add_attribute('ID', 'mRNA1')
    assert c1.get_attribute('Parent') == 'mRNA1'
    c1.type = 'exon'
    assert c1.type == 'exon'


def test_multi():
    """
    Test handling of multi-features.

    In GFF3 discontiguous sequence features are often encoded using multi-
    features: that is, a single feature described across mulitple entries
    (lines). This unit test validates the handling of multi-features.
    """
    gene = eden()
    for feature in gene:
        if feature.type == 'CDS':
            assert feature.is_multi
            if feature.get_attribute('ID') == 'cds00001':
                assert feature.multi_rep._region == Range(1200, 1500)
            elif feature.get_attribute('ID') == 'cds00002':
                assert feature.multi_rep._region == Range(1200, 1500)
            elif feature.get_attribute('ID') == 'cds00003':
                assert feature.multi_rep._region == Range(3300, 3902)
            elif feature.get_attribute('ID') == 'cds00004':  # pragma no cover
                assert feature.multi_rep._region == Range(3390, 3902)
        else:
            assert not feature.is_multi and feature.multi_rep is None

    for feature in gene:
        if feature.get_attribute('ID') == 'cds00004':
            if feature.multi_rep == feature:
                feature.type = 'coding sequence'
            assert feature.type == 'coding sequence'


def test_compare():
    """Test comparison and sorting of features."""
    gff3 = ['chr1', 'vim', 'gene', '1000', '2000', '.', '+', '.', '.']
    g1 = Feature('\t'.join(gff3))
    gff3 = ['chr1', 'vim', 'gene', '3000', '4000', '.', '+', '.', '.']
    g2 = Feature('\t'.join(gff3))
    gff3 = ['chr1', 'vim', 'gene', '3000', '4500', '.', '+', '.', '.']
    g3 = Feature('\t'.join(gff3))

    assert g1 < g2
    assert g1 <= g3
    assert g3 >= g1
    assert not g1 > g3
    assert g3 >= g2
    assert not g3 <= g2
    assert sorted([g3, g2, g1]) == [g1, g2, g3]

    gff3 = ['chr10', 'vim', 'gene', '100', '400', '.', '-', '.', '.']
    g4 = Feature('\t'.join(gff3))
    gff3 = ['chr2', 'vim', 'gene', '2000', '2500', '.', '-', '.', '.']
    g5 = Feature('\t'.join(gff3))

    assert g2 <= g4
    assert not g4 <= g2
    assert g4 < g5
    assert not g5 < g4
    assert g5 > g1
    assert not g1 > g5
    assert g5 >= g1
    assert not g1 >= g5
    assert g5 >= g5
    assert g5 <= g5
    assert sorted([g3, g5, g1, g4, g2]) == [g1, g2, g3, g4, g5]

    d = Directive('##gff-version')
    c = Comment('# Cool story, bro!')
    assert g1 > c
    assert g1 >= c
    assert g2 > d
    assert g2 >= d
    assert not g3 < c
    assert not g3 <= c
    assert not g4 < d
    assert not g4 <= d

    gff3 = ['chr', 'vim', 'mRNA', '1000', '2000', '.', '+', '.', 'ID=mRNA1']
    f1 = Feature('\t'.join(gff3))
    gff3 = ['chr', 'vim', 'tRNA', '1000', '2000', '.', '+', '.', 'ID=tRNA1']
    f2 = Feature('\t'.join(gff3))

    assert f1 < Sequence('>contig1', 'GATTACA')
    assert f1 <= Sequence('>contig1', 'GATTACA')
    assert not f1 > Sequence('>contig1', 'GATTACA')
    assert not f1 >= Sequence('>contig1', 'GATTACA')
    assert f1 > f2
    assert f1 >= f2


def test_cyclic():
    """Test handling of cyclic features."""
    gene = eden()
    for feature in gene:
        if feature.get_attribute('ID') == 'exon00002':
            feature.add_child(gene, rangecheck=False)

    try:
        gff3string = repr(gene)
    except Exception:
        pass