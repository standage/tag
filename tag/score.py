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

    def __init__(self, datastr):
        self._type = None
        if datastr == '.':
            self.score = None
        elif re.search('^-*\d+$', datastr):
            self.score = int(datastr)
            self._type = int
        else:
            self.score = float(datastr)
            self._type = float

    def __str__(self):
        if self.score is None:
            return '.'
        elif self._type == int:
            return '{:d}'.format(self.score)
        elif abs(self.score) < 1e6 and abs(self.score) > 1e-4:
            return '{:1.3f}'.format(self.score)
        else:
            return '{:1.3E}'.format(self.score)
