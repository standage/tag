#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------
"""Package-wide configuration"""

from . import comment
from . import feature
from . import range
from . import reader
from gzip import open as gzopen

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


def file_handle(filename, mode):
    if mode not in ['r', 'w']:
        raise ArgumentError('invalid mode "{}"'.format(mode))
    openfunc = open
    if filename.endswith('.gz'):
        openfunc = gzopen
        mode += 't'
    return openfunc(filename, mode)
