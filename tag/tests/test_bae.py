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


def test_cli(capsys):
    arglist = ['bae', data_file('Ye.region02.gff3')]
    args = tag.cli.parser().parse_args(arglist)
    tag.cli.bae.main(args)
    terminal = capsys.readouterr()
    with open(data_file('Ye.region02.bae.gff3'), 'r') as fh:
        testout = fh.read()
    assert terminal.out.strip() == testout.strip()
