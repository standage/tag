#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2017 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import re


class Score(object):
    @staticmethod
    def from_str(datastr):
        if datastr == '.':
            return Score(None)
        elif re.search(r'^-*\d+$', datastr):
            return Score(int(datastr))
        else:
            return Score(float(datastr))

    def __init__(self, data):
        if isinstance(data, str):
            raise TypeError(
                'please convert score to a numeric type or instantiate the '
                'object from the `tag.Score.from_str` function'
            )
        self.value = data
        self._type = type(data)

    def __str__(self):
        if self.value is None:
            return '.'
        elif self._type == int:
            return '{:d}'.format(self.value)
        elif abs(self.value) < 1e6 and abs(self.value) > 1e-4:
            return '{:1.3f}'.format(self.value)
        else:
            return '{:1.3E}'.format(self.value)
