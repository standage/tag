#!/usr/bin/env python
#
# ------------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# ------------------------------------------------------------------------------


class Comment(object):
    """
    Represents a comment in an annotation (GFF3) file.

    Any GFF3 entry starting with >= 1 '#' characters is treated as a comment,
    with two exceptions:

    - the separator directive, a line containing '###' and nothing more
    - any entry beginning with just two '#' characters is treated as a
      directive.
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
