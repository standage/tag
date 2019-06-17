#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2017 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from tag import Score


def test_basic():
    for score in ['.', '10', '12.345', '-11', '-98765.432', '-4.555E+09']:
        assert str(Score(score)) == score
    assert str(Score('10.0')) == '10.000'
    assert str(Score('1.32e12')) == '1.320E+12'
    assert str(Score('1.2e-16')) == '1.200E-16'
