#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import tag


def test_features():
    ifn = 'tests/testdata/pbar-withseq.gff3'

    reader = tag.reader.GFF3Reader(infilename=ifn)
    features = [f for f in tag.select.features(reader)]
    assert len(features) == 3

    reader = tag.reader.GFF3Reader(infilename=ifn)
    features = [f for f in tag.select.features(reader, type='pseudogene')]
    assert len(features) == 1

    reader = tag.reader.GFF3Reader(infilename=ifn)
    with pytest.raises(ValueError) as ve:
        features = [f for f in tag.select.features(reader, traverse=True)]
    assert 'cannot traverse without a specific feature type' in str(ve)


def test_directives():
    ifn = 'tests/testdata/pbar-withseq.gff3'

    reader = tag.reader.GFF3Reader(infilename=ifn)
    directives = [d for d in tag.select.directives(reader)]
    assert len(directives) == 3

    reader = tag.reader.GFF3Reader(infilename=ifn)
    dirs = [d for d in tag.select.directives(reader, type='sequence-region')]
    assert len(dirs) == 2
    assert [d.seqid for d in dirs] == ['NW_011929623.1', 'NW_011929624.1']


def test_sequences():
    ifn = 'tests/testdata/pbar-withseq.gff3'

    reader = tag.reader.GFF3Reader(infilename=ifn)
    seqs = [s for s in tag.select.sequences(reader)]
    assert len(seqs) == 2
    assert sorted([s.seq[0:5] for s in seqs]) == ['GCTAA', 'TAGAC']
