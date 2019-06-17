#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import pytest
import tag
from tag import GFF3Reader
from tag.tests import data_file, data_stream


def test_features():
    reader = GFF3Reader(infilename=data_file('pbar-withseq.gff3'))
    features = [f for f in tag.select.features(reader)]
    assert len(features) == 3

    reader = GFF3Reader(infilename=data_file('pbar-withseq.gff3'))
    features = [f for f in tag.select.features(reader, type='pseudogene')]
    assert len(features) == 1

    reader = GFF3Reader(infilename=data_file('pbar-withseq.gff3'))
    with pytest.raises(ValueError) as ve:
        features = [f for f in tag.select.features(reader, traverse=True)]
    assert 'cannot traverse without a specific feature type' in str(ve)


def test_window_contain():
    reader = GFF3Reader(infilename=data_file('pdom-withseq.gff3'))
    features = tag.select.features(reader)
    window = tag.select.window(features, 'PdomSCFr1.2-0483', 3000, 7000)
    selected = [f for f in window]
    assert len(selected) == 1
    assert selected[0].get_attribute('Name') == 'PdomGENEr1.2-04310'


def test_window_overlap():
    reader = GFF3Reader(infilename=data_file('osat-twoscaf.gff3.gz'))
    features = tag.select.features(reader)
    window = tag.select.window(features, 'NW_015379189.1', 15000, 20000,
                               strict=False)
    selected = [f for f in window]
    assert len(selected) == 3
    assert [f.type for f in selected] == ['region', 'gene', 'gene']


def test_window_seqid():
    reader = GFF3Reader(infilename=data_file('osat-twoscaf.gff3.gz'))
    features = tag.select.features(reader)
    window = tag.select.window(features, 'NW_015379208.1')
    selected = [f for f in window]
    assert len(selected) == 2
    assert [f.type for f in selected] == ['region', 'gene']


def test_directives():
    reader = GFF3Reader(infilename=data_file('pbar-withseq.gff3'))
    directives = [d for d in tag.select.directives(reader)]
    assert len(directives) == 3

    reader = GFF3Reader(infilename=data_file('pbar-withseq.gff3'))
    dirs = [d for d in tag.select.directives(reader, type='sequence-region')]
    assert len(dirs) == 2
    assert [d.seqid for d in dirs] == ['NW_011929623.1', 'NW_011929624.1']


def test_sequences():
    reader = GFF3Reader(infilename=data_file('pbar-withseq.gff3'))
    seqs = [s for s in tag.select.sequences(reader)]
    assert len(seqs) == 2
    assert sorted([s.seq[0:5] for s in seqs]) == ['GCTAA', 'TAGAC']


def test_merge():
    merger = tag.select.merge(
        GFF3Reader(infilename=data_file('ex-red-1.gff3')),
        GFF3Reader(infilename=data_file('ex-red-2.gff3')),
        GFF3Reader(infilename=data_file('ex-red-3.gff3')),
    )
    merge_stream = tag.select.features(merger)
    test_stream = tag.select.features(
        GFF3Reader(infilename=data_file('ex-red-merged.gff3'))
    )
    for mergefeat, testfeat in zip(merge_stream, test_stream):
        print(mergefeat, testfeat, sep='\n', end='\n\n')
        assert mergefeat.like(testfeat)
