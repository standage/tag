#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the ISC license: see LICENSE.txt.
# -----------------------------------------------------------------------------


class Comment():
    """
    Represents a comment in an annotation (GFF3) file.

    Any GFF3 entry starting with >= 1 '#' characters is treated as a comment,
    with two exceptions: first, the separator directive, a line containing
    '###' and nothing more; second, any entry beginning with just two '#'
    characters is treated as a directive.
    """

    def __init__(self, data):
        assert data.startswith('#')
        self._rawdata = data

    def __repr__(self):
        return self._rawdata

    def __str__(self):
        i = 0
        while self._rawdata[i] in '# ':
            i += 1
        return self._rawdata[i:]

    def __lt__(self, other):
        from .directive import Directive
        from .feature import Feature
        if isinstance(other, Directive):
            return False
        if isinstance(other, Feature):
            return True
        return self._rawdata < other._rawdata

    def __gt__(self, other):
        from .directive import Directive
        from .feature import Feature
        if isinstance(other, Directive):
            return True
        if isinstance(other, Feature):
            return False
        return self._rawdata > other._rawdata


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------

def test_init():
    """[aeneas::Comment] Test constructor."""
    c1 = Comment('# A strange game. The only winning move is not to play. '
                 'How about a nice game of chess?')
    try:
        c2 = Comment('This is not a valid comment')
    except AssertionError:
        pass
    assert c1._rawdata == ('# A strange game. The only winning move is not to '
                           'play. How about a nice game of chess?')


def test_repr():
    """[aeneas::Comment] Test default representation."""
    c1 = Comment('# This gene model is a fragment')
    c2 = Comment('############## Ignore below this point.')

    assert "%r" % c1 == '# This gene model is a fragment'
    assert "%r" % c2 == '############## Ignore below this point.'


def test_str():
    """[aeneas::Comment] Test string representation."""
    c1 = Comment('# This gene model is a fragment')
    c2 = Comment('############## Ignore below this point.')

    assert "%s" % c1 == 'This gene model is a fragment'
    assert "%s" % c2 == 'Ignore below this point.'


def test_sort():
    """[aeneas::Comment] Test sorting and comparison."""
    from .directive import Directive
    from .feature import Feature

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
