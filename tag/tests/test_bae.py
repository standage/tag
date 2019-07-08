#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2019 Battelle National Biodefense Institute.
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import tag
from tag.tests import data_file


def test_baesic():
    instream = tag.GFF3Reader(infilename=data_file('Ye.region01.gff3'))
    locusstream = tag.locus.loci(instream)
    baestream = tag.bae.eval_stream(locusstream)
    features = list(tag.select.features(baestream, type='CDS'))
    assert len(features) == 5
    assert [f.get_attribute('locus_orfs') for f in features] == [4, 4, 4, 4, 1]
    cov = [f.get_attribute('protein_coverage') for f in features]
    assert cov == [1.0, 1.0, 1.0, 1.0, 0.0]


def test_bae_cli(capsys):
    arglist = ['bae', data_file('Ye.region02.gff3')]
    args = tag.cli.parser().parse_args(arglist)
    tag.cli.bae.main(args)
    terminal = capsys.readouterr()
    with open(data_file('Ye.region02.bae.gff3'), 'r') as fh:
        testout = fh.read()
    assert terminal.out.strip() == testout.strip()


def test_collapse():
    instream = tag.GFF3Reader(infilename=data_file('bork.gff3'))
    locusstream = tag.locus.loci(instream)
    loci = [f for s, c, f in locusstream]

    feats = list(tag.bae.collapse_locus(loci[0]))
    assert len(feats) == 1
    assert feats[0].get_attribute('support_from') == 'BBTools'

    feats = list(tag.bae.collapse_locus(loci[1]))
    assert len(feats) == 5
    supports = [
        f.get_attribute('support_from', as_string=True)
        for f in tag.select.features(feats, type='CDS')
    ]
    test_supports = [
        'AUGUSTUS', 'GeneMark.hmm,MetaGeneAnnotator,Prodigal_v2.6.3',
        'Glimmer3', 'BBTools'
    ]
    print(supports, test_supports)
    assert sorted(supports) == sorted(test_supports)

    feats = list(tag.bae.collapse_locus(loci[2]))
    assert len(feats) == 4
    supports = [
        f.get_attribute('support_from', as_string=True)
        for f in tag.select.features(feats, type='CDS')
    ]
    test_supports = [
        'GeneMark.hmm,Glimmer3', 'AUGUSTUS,MetaGeneAnnotator',
        'BBTools,Prodigal_v2.6.3'
    ]
    assert sorted(supports) == sorted(test_supports)

    feats = list(tag.bae.collapse_locus(loci[3]))
    assert len(feats) == 2
    supports = [
        f.get_attribute('support_from', as_string=True)
        for f in tag.select.features(feats, type='CDS')
    ]
    test_supports = ['AUGUSTUS,BBTools,GeneMark.hmm,Glimmer3,Prodigal_v2.6.3']
    assert sorted(supports) == sorted(test_supports)
    assert feats[1].score == 5


def test_bcollapse_cli(capsys):
    arglist = ['bcollapse', data_file('bork.gff3')]
    args = tag.cli.parser().parse_args(arglist)
    tag.cli.bcollapse.main(args)
    terminal = capsys.readouterr()
    testout = open(data_file('bork-out.gff3'), 'r').read()
    assert terminal.out.strip() == testout.strip()
