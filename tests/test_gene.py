#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2017 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import pytest
import tag


def test_exon_dedup():
    output = StringIO()
    reader = tag.GFF3Reader(tag.pkgdata('turt-iso.gff3'))
    dedup = tag.gene.exon_dedup(reader)
    exons = [e for e in tag.select.features(dedup, type='exon')]
    assert len(exons) == 7
