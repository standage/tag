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


def test_entry():
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
