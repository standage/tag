#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import pytest
from aeneas.comment import Comment
from aeneas.directive import Directive
from aeneas.feature import Feature


def test_init():
    """Test constructor."""
    c1 = Comment('# A strange game. The only winning move is not to play. '
                 'How about a nice game of chess?')
    assert c1._rawdata == ('# A strange game. The only winning move is not to '
                           'play. How about a nice game of chess?')
    with pytest.raises(AssertionError):
        c2 = Comment('This is not a valid comment')


def test_repr():
    """Test default representation."""
    c1 = Comment('# This gene model is a fragment')
    c2 = Comment('############## Ignore below this point.')

    assert repr(c1) == '# This gene model is a fragment'
    assert repr(c2) == '############## Ignore below this point.'


def test_str():
    """Test string representation."""
    c1 = Comment('# This gene model is a fragment')
    c2 = Comment('############## Ignore below this point.')

    assert str(c1) == 'This gene model is a fragment'
    assert str(c2) == 'Ignore below this point.'


def test_sort():
    """Test sorting and comparison."""
    c1 = Comment('# This gene model is a fragment')
    c2 = Comment('############## Ignore below this point.')
    d = Directive('##sequence-region chr 1 1000')
    gff3 = ['chr', 'vim', 'mRNA', '1001', '1420', '.', '+', '.', 'ID=t1']
    f = Feature('\t'.join(gff3))

    assert c1 > d
    assert not c1 < d
    assert c1 < f
    assert not c1 > f
    assert c1 < c2
    assert c2 > c1
