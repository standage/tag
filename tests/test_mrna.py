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


def test_primary():
    reader = tag.reader.GFF3Reader(tag.pkgdata('nanosplice.gff3'))
    gene = next(tag.select.features(tag.mrna.primary(reader), types=['gene']))
    assert gene.cdslen is None
    assert gene.num_children == 1
    assert gene.children[0].get_attribute('ID') == 'mRNAsecond'

    reader = tag.reader.GFF3Reader(tag.pkgdata('pdom-withseq.gff3'))
    for gene in tag.select.features(tag.mrna.primary(reader), types=['gene']):
        assert gene.num_children == 1

    reader = tag.reader.GFF3Reader(tag.pkgdata('psyllid-100k.gff3'))
    for gene in tag.select.features(tag.mrna.primary(reader), types=['gene']):
        mrnas = [f for f in gene if f.type == 'mRNA']
        assert len(mrnas) <= 1, mrnas
