#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import sys
import tag
from tag.tests import data_file, data_stream


def test_tag_open_read():
    fh = tag.open(data_file('boring-genes.gff3'), 'r')
    line = next(fh)
    expout = 'chr1\tnano\tgene\t6615\t11733\t0.72\t+\t.\tID=gEnE1;Name=Adh\n'
    assert line == expout

    with pytest.raises(IOError):
        fh = tag.open('tag/bogus.py', 'r')

    fh = data_stream('gzipdata.gff3.gz')
    line = next(fh)
    assert line == '##gff-version 3\n'

    with pytest.raises(IOError):
        fh = tag.open('tag/bogus.gz', 'r')


def test_tag_open_badmode():
    with pytest.raises(ValueError) as ae:
        fh = tag.open('tag/__init__.py', 'z')
    assert 'invalid mode' in str(ae)
