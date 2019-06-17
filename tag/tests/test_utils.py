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
    fh = tag.open('tag/__init__.py', 'r')
    line = next(fh)
    assert line == '#!/usr/bin/env python\n'

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
