#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2017 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
from tag import Score


@pytest.mark.parametrize('score,scorestr', [
    (10, '10'),
    (12.345, '12.345'),
    (-11, '-11'),
    (-98765.432, '-98765.432'),
    (-4.555E+09, '-4.555E+09'),
    (10.0, '10.000'),
    (1.32e12, '1.320E+12'),
    (1.2e-16, '1.200E-16'),
])
def test_basic(score, scorestr):
    scoreobj = Score(score)
    assert scoreobj.value == score
    assert str(scoreobj) == scorestr


def test_bad_score():
    with pytest.raises(TypeError) as te:
        score = Score('ZaphodBeeblebrox')
    assert 'please convert score to a numeric type' in str(te)


def test_from_str():
    for score in ['.', '10', '12.345', '-11', '-98765.432', '-4.555E+09']:
        assert str(Score.from_str(score)) == score
    assert str(Score.from_str('10.0')) == '10.000'
    assert str(Score.from_str('1.32e12')) == '1.320E+12'
    assert str(Score.from_str('1.2e-16')) == '1.200E-16'
