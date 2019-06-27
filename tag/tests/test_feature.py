#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
from tag import Comment
from tag import Directive
from tag import Feature
from tag import Range
from tag import Score
from tag import Sequence
from tag import select
from tag.tests import data_file, data_stream


@pytest.fixture
def eden():
    """Data fixture for unit tests.

    This function is not a good example of how to parse arbitrary GFF3 files,
    it's quick and dirty for this particular case.
    """
    edenfile = open(data_file('eden.gff3'), 'r')
    edengff3 = edenfile.readlines()
    gene = None
    mRNAs = list()
    cdsreps = dict()
    for line in edengff3:
        if line.startswith('#'):
            continue
        feature = Feature.from_gff3(line.rstrip())
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
    """Test basic constructor and GFF3 factory constructor."""
    f0 = Feature('chr', 'gene', 1000, 2000, strand='+')
    gff3 = ['chr', 'tag', 'gene', '1001', '2000', '.', '+', '.', '.']
    f1 = Feature.from_gff3('\t'.join(gff3))
    assert f0.like(f1)
    assert str(f1) == '\t'.join(gff3)
    assert repr(f1) == '\t'.join(gff3)
    assert len(f1) == 1000
    assert f1.slug == 'gene@chr[1001, 2000]'

    gff3[7] = '-1'
    with pytest.raises(ValueError) as ve:
        f999 = Feature.from_gff3('\t'.join(gff3))
    assert 'invalid phase "-1"' in str(ve)
    gff3[6] = '$'
    with pytest.raises(ValueError) as ve:
        f999 = Feature.from_gff3('\t'.join(gff3))
    assert 'invalid strand "$"' in str(ve)

    f2 = Feature(
        'chr', 'gene', 1000, 2000, strand='+', attrstr='ID=gene1;Name=EDEN'
    )
    f2.add_attribute('perc', 0.75)
    assert str(f2) == (
        'chr\ttag\tgene\t1001\t2000\t.\t+\t.\tID=gene1;Name=EDEN;perc=0.7500'
    )

    gff3[5] = 0.28466
    f3 = Feature('chr', 'gene', 1000, 2000, strand='+', score=0.28466)
    assert f3.score == pytest.approx(0.28466)
    assert '\t0.285\t' in str(f3)


def test_seqid():
    """Test seqid handling."""
    f1 = Feature('chr', 'gene', 999, 2000, strand='+')
    assert f1.seqid == 'chr'
    f1.seqid = 'scaffold_28'
    assert f1.seqid == 'scaffold_28'

    f2 = Feature('scaffold_28', 'mRNA', 999, 2000, strand='+')
    f3 = Feature('scaffold_28', 'exon', 999, 2000, strand='+')
    f2.add_child(f3)
    f1.add_child(f2)
    for feature in f1:
        assert feature.seqid == 'scaffold_28'
    f1.seqid = 'ctg123'
    for feature in f1:
        assert feature.seqid == 'ctg123'


def test_source():
    """Test source handling."""
    f0 = Feature('chr', 'gene', 999, 2000, strand='+')
    assert f0.source == 'tag'

    f1 = Feature('chr', 'gene', 999, 2000, strand='+', source='vim')
    assert f1.source == 'vim'
    f1.source = 'nano'
    assert f1.source == 'nano'

    f2 = Feature('chr', 'mRNA', 999, 2000, strand='+', source='emacs')
    f3 = Feature('chr', 'exon', 999, 2000, strand='+', source='vim')
    f2.add_child(f3)
    f1.add_child(f2)

    # Test cascading when child sources are mixed
    f1.source = 'gedit'
    assert f1.source == 'gedit'
    for feature in select.features(f1, type='mRNA'):
        assert feature.source == 'emacs'
    for feature in select.features(f1, type='exon'):
        assert feature.source == 'vim'

    # Test cascading when child sources match parent
    for f in f1:
        f.source = 'atom'
    f1.source = 'Sublime'
    for feature in f1:
        assert feature.source == 'Sublime'


def test_type():
    """Test type handling."""
    f1 = Feature('chr', 'mRNA', 999, 2000, strand='+')
    assert f1.type == 'mRNA'
    f2 = Feature('chr', 'messenger RNA', 999, 2000, strand='+')
    assert f2.type == 'messenger RNA'


def test_region():
    """Test coordinate handling."""
    f1 = Feature('chr', 'gene', 999, 2000, strand='+', attrstr='ID=g1')
    assert f1._range == Range(999, 2000)
    assert f1.start == 999 and f1.end == 2000
    assert f1.contains(Range(1500, 1600))
    assert f1.contains_point(5) is False
    assert f1.overlap(Range(2000, 2001)) is False
    assert f1.overlap(Range(1999, 2001)) is True

    f2 = Feature('contig5', 'mRNA', 499, 2500, strand='+', attrstr='Parent=g1')
    with pytest.raises(AssertionError) as ae:
        f1.add_child(f2, rangecheck=True)
    assert 'seqid mismatch for feature g1' in str(ae)

    f2.seqid = 'chr'
    with pytest.raises(AssertionError) as ae:
        f1.add_child(f2, rangecheck=True)
    assert 'is not contained within its span' in str(ae)

    f1.add_child(f2)
    assert len(f1.children) == 1
    assert f2._range == Range(499, 2500)

    f2.set_coord(999, 2000)
    assert f2.start == 999 and f2.end == 2000

    f1.transform(100000)
    for feature in f1:
        assert feature._range == Range(100999, 102000)

    f1.transform(100000, newseqid='scf89')
    for feature in f1:
        assert feature.seqid == 'scf89'
        assert feature._range == Range(200999, 202000)


def test_score():
    """Test score handling."""
    gff3 = ['chr', 'vim', 'EST_match', '57229', '57404', '.', '+', '.', '.']
    f1 = Feature.from_gff3('\t'.join(gff3))
    assert f1.score is None

    f2 = Feature('chr', 'EST_match', 57228, 57404, strand='+', score=0.97)
    assert f2.score == pytest.approx(0.97)
    assert '\t0.970\t' in str(f2)

    f3 = Feature(
        'chr', 'EST_match', 57228, 57404, strand='+', score=Score(-1.8332),
    )
    assert f3.score == pytest.approx(-1.8332)
    assert '\t-1.833\t' in str(f3)
    f3.score = 42
    assert f3.score == 42
    with pytest.raises(TypeError) as te:
        f3.score = 'BogusScore'
    assert 'please convert score to a numeric type' in str(te)


def test_strand():
    """Test strand handling."""
    gff3 = ['chr', 'vim', 'EST_match', '57229', '57404', '.', '.', '.', '.']
    f1 = Feature('chr', 'EST_match', 57228, 57404)
    assert f1.strand == '.'

    f2 = Feature('chr', 'EST_match', 57228, 57404, strand='-')
    f3 = Feature('chr', 'match_part', 57228, 57298, strand='-')
    assert f2.strand == '-'
    assert f3.strand == '-'

    f4 = Feature('chr', 'match_part', 57376, 57404, strand='+')
    assert f4.strand == '+'

    f2.add_child(f3)
    with pytest.raises(AssertionError) as ae:
        f2.add_child(f4, rangecheck=True)
    assert 'has a different strand' in str(ae)
    assert len(f2.children) == 1
    f4._strand = '-'
    f2.add_child(f4)
    assert len(f2.children) == 2

    with pytest.raises(ValueError) as ve:
        f5 = Feature('chr', 'EST_match', 57228, 57404, strand='~')
    assert 'invalid strand "~"' in str(ve)
    f5 = Feature('chr', 'EST_match', 57228, 57404)
    with pytest.raises(ValueError) as ve:
        f5.strand = '~'
    assert 'invalid strand "~"' in str(ve)
    f5.strand = '-'
    assert f5.strand == '-'


def test_phase():
    """Test phase handling."""
    m1 = Feature('chr', 'mRNA', 1000, 1420, strand='+')
    assert m1.phase is None
    c1 = Feature('chr', 'CDS', 1000, 1100, strand='+', phase=0)
    assert c1.phase == 0
    c2 = Feature('chr', 'CDS', 1201, 1236, strand='+', phase=2)
    c1.add_sibling(c2)
    assert c2.phase == 2

    c3 = Feature('chr', 'CDS', 1300, 1364, strand='+', phase=2)
    c1.add_sibling(c3)
    assert c3.phase == 2
    c4 = Feature('chr', 'CDS', 1400, 1420, strand='+', phase=1)
    c1.add_sibling(c4)
    assert c4.phase == 1

    with pytest.raises(ValueError) as ve:
        c5 = Feature('chr', 'CDS', 1400, 1420, strand='+', phase='4')
    assert 'invalid phase "4"' in str(ve)


def test_attributes(eden):
    """Test attribute handling."""
    assert eden.get_attribute('ID') == 'gene00001'
    assert eden.get_attribute('ID', as_list=True) == 'gene00001'
    assert eden.get_attribute('Name') == 'EDEN'
    assert eden.get_attribute('Name', as_list=True) == ['EDEN']
    assert eden.get_attribute_keys() == ['ID', 'Name']

    eden.add_attribute('Name', 'Gandalf')
    assert eden.get_attribute('Name') == 'Gandalf'
    eden.add_attribute('Name', 'Aragorn', append=True)
    assert eden.get_attribute('Name') == ['Aragorn', 'Gandalf']
    assert eden.get_attribute('Name', as_string=True) == 'Aragorn,Gandalf'

    eden.add_attribute('ID', 'g1')
    for child in eden.children:
        assert child.get_attribute('Parent') == 'g1'

    for feature in eden:
        if feature.type == 'exon' and feature._range == Range(4999, 5500):
            assert feature.get_attribute('Parent') == [
                'mRNA00001',
                'mRNA00002',
                'mRNA00003']
            assert feature.get_attribute('Parent', as_string=True) == \
                'mRNA00001,mRNA00002,mRNA00003'

    for feature in eden:
        if feature.get_attribute('ID') == 'mRNA00003':
            feature.add_attribute('ID', 'mRNA3')

    assert repr(eden) == \
        open(data_file('eden-mod.gff3'), 'r').read().rstrip()

    eden.add_attribute('Note', 'I need to test access of other attributes')
    assert eden.attributes == ('ID=g1;Name=Aragorn,Gandalf;Note='
                               'I need to test access of other attributes')
    eden.drop_attribute('BogusAttributeButThatIsOkayNoHarmDone')
    eden.drop_attribute('Name')
    assert eden.get_attribute('Name') is None

    gff3 = ['chr', 'vim', 'mRNA', '1001', '1420', '.', '+', '.', '.']
    m1 = Feature.from_gff3('\t'.join(gff3))
    assert m1.phase is None
    m1.add_attribute('ID', 't1')
    assert m1.get_attribute('ID') == 't1'

    gff3 = ['chr', 'vim', 'CDS', '1001', '1420', '.', '+', '.', '.']
    c1 = Feature.from_gff3('\t'.join(gff3))
    m1.add_child(c1)
    m1.add_attribute('ID', 'mRNA1')
    assert c1.get_attribute('Parent') == 'mRNA1'
    c1.type = 'exon'
    assert c1.type == 'exon'


def test_multi(eden):
    """Test handling of multi-features.

    In GFF3 discontiguous sequence features are often encoded using multi-
    features: that is, a single feature described across mulitple entries
    (lines). This unit test validates the handling of multi-features.
    """
    for feature in eden:
        if feature.type == 'CDS':
            assert feature.is_multi
            if feature.get_attribute('ID') == 'cds00001':
                assert feature.multi_rep._range == Range(1200, 1500)
            elif feature.get_attribute('ID') == 'cds00002':
                assert feature.multi_rep._range == Range(1200, 1500)
            elif feature.get_attribute('ID') == 'cds00003':
                assert feature.multi_rep._range == Range(3300, 3902)
            elif feature.get_attribute('ID') == 'cds00004':  # pragma no cover
                assert feature.multi_rep._range == Range(3390, 3902)
        else:
            assert not feature.is_multi and feature.multi_rep is None

    for feature in eden:
        if feature.get_attribute('ID') == 'cds00004':
            if feature.multi_rep == feature:
                feature.type = 'coding sequence'
            assert feature.type == 'coding sequence'


def test_compare():
    """Test comparison and sorting of features."""
    g1 = Feature('chr1', 'gene', 999, 2000, strand='+')
    g2 = Feature('chr1', 'gene', 2999, 4000, strand='+')
    g3 = Feature('chr1', 'gene', 2999, 4500, strand='+')

    assert g1 < g2
    assert g2 > g1
    assert g1 <= g3
    assert g3 >= g1
    assert not g1 > g3
    assert g3 >= g2
    assert not g3 <= g2
    assert sorted([g3, g2, g1]) == [g1, g2, g3]

    g4 = Feature('chr10', 'gene', 99, 400, strand='-')
    g5 = Feature('chr2', 'gene', 1999, 2500, strand='-')
    g6 = Feature('chr2', 'exon', 1999, 2500, strand='-')

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
    assert g5 < g6
    assert g6 >= g5
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
    f1 = Feature('chr', 'mRNA', 999, 2000, strand='+')
    gff3 = ['chr', 'vim', 'tRNA', '1000', '2000', '.', '+', '.', 'ID=tRNA1']
    f2 = Feature('chr', 'tRNA', 999, 2000, strand='+')

    assert f1 < Sequence('>contig1', 'GATTACA')
    assert f1 <= Sequence('>contig1', 'GATTACA')
    assert not f1 > Sequence('>contig1', 'GATTACA')
    assert not f1 >= Sequence('>contig1', 'GATTACA')
    assert f1 > f2
    assert f1 >= f2


def test_cyclic(eden):
    """Test handling of cyclic features."""
    for feature in eden:
        if feature.get_attribute('ID') == 'exon00002':
            feature.add_child(eden, rangecheck=False)
    with pytest.raises(ValueError) as ve:
        gff3string = repr(eden)
    assert 'feature graph is cyclic' in str(ve)


def test_pseudo_1():
    feat_x = Feature('1', 'cDNA_match', 999, 2000, source='atom', strand='+')
    feat_y = Feature('1', 'cDNA_match', 2999, 4000, source='atom', strand='+')
    feat_z = Feature('1', 'cDNA_match', 4999, 6000, source='atom', strand='+')
    feat_x.add_sibling(feat_y)
    feat_x.add_sibling(feat_z)
    parent = feat_x.pseudoify()
    assert str(parent) == ''
    assert parent.is_toplevel is True
    assert parent.slug == 'cDNA_match@1[1000, 6000]'
    assert repr(parent) == repr(feat_x)
    assert repr(parent) == repr(feat_z.pseudoify())


def test_ncbi_geneid():
    gff3 = ('NW_007377440.1	Gnomon	gene	63775	73670	.	-	.	'
            'ID=gene2;Name=LOC103504972;Dbxref=GeneID:103504972;gbkey=Gene;'
            'gene=LOC103504972;gene_biotype=protein_coding')
    gene = Feature.from_gff3(gff3)
    assert gene.ncbi_geneid == '103504972'

    gff3 = ('NW_007378253.1	Gnomon	mRNA	103380	167368	.	-	.	'
            'ID=mRNA10000;Parent=gene9300;Name=XM_008477076.2;'
            'Dbxref=Genbank:XM_008477076.2,GeneID:103512317;gbkey=mRNA;'
            'gene=LOC103512317;model_evidence=Supporting evidence includes '
            'similarity to: 5 Proteins%2C and 84%25 coverage of the annotated '
            'genomic feature by RNAseq alignments;'
            'product=glutamine--fructose-6-phosphate aminotransferase '
            '[isomerizing] 2-like;transcript_id=XM_008477076.2')
    mrna = Feature.from_gff3(gff3)
    assert mrna.ncbi_geneid == '103512317'

    gff3 = ('NW_007377513.1	RefSeq	cDNA_match	271974	274535	.	+	.	'
            'ID=cDNA_match42;Gap=M2086 D2 M474;Target=XM_008486908.2 1 2560 +;'
            'for_remapping=2;gap_count=1;num_ident=2840;num_mismatch=0;'
            'pct_coverage=100;pct_coverage_hiqual=100;'
            'pct_identity_gap=99.9296;pct_identity_ungap=100;rank=1')
    match = Feature.from_gff3(gff3)
    assert match.ncbi_geneid is None

    gff3 = ('chr	atom	region	1000	2000	.	.	.	'
            'Dbxref=MyDB:ID12345')
    region = Feature.from_gff3(gff3)
    assert region.ncbi_geneid is None


def test_attribute_crawl():
    from tag import GFF3Reader
    from tag import select

    reader = GFF3Reader(data_stream('pbar-withseq.gff3'))
    features = select.features(reader)
    gene = next(features)
    names = ['XP_011630649.1', 'LOC105422795', 'XM_011632347.1']
    assert gene.attribute_crawl('Name') == set(names)


def test_feature_like():
    f1 = Feature('chr1', 'mRNA', 999, 2000, source='vim', strand='+')
    f2 = Feature('chr2', 'tRNA', 2999, 4000, source='emacs', strand='+')
    f3 = Feature('chr1', 'tRNA', 2999, 4000, source='emacs', strand='+')
    f4 = Feature('chr1', 'mRNA', 999, 2000, source='emacs', strand='+')
    f5 = Feature('chr1', 'mRNA', 999, 2000, source='vim', strand='+')
    d1 = Directive('##gff-version   3')
    assert f1.like(f2) is False
    assert f1.like(f3) is False
    assert f1.like(f4) is False
    assert f1.like(f5) is True
    assert f1.like(d1) is False
