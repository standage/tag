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


def test_tag_open_read():
    fh = tag.open('tag/__init__.py', 'r')
    assert hasattr(fh, 'mode')

    with pytest.raises(IOError):
        fh = tag.open('tag/bogus.py', 'r')

    fh = tag.open('tests/testdata/gzipdata.gff3.gz', 'r')
    assert not hasattr(fh, 'mode')

    with pytest.raises(IOError):
        fh = tag.open('tag/bogus.gz', 'r')


def test_tag_open_badmode():
    with pytest.raises(ValueError) as ae:
        fh = tag.open('tag/__init__.py', 'z')
    assert 'invalid mode' in str(ae)
