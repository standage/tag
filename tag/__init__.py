#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------
"""Package-wide configuration"""

try:
    import __builtin__ as builtins
except ImportError:  # pragma: no cover
    import builtins
from tag.comment import Comment
from tag.directive import Directive
from tag.feature import Feature
from tag.sequence import Sequence
from tag.range import Range
from tag.reader import GFF3Reader
from tag.writer import GFF3Writer
from tag.score import Score
from tag import bae
from tag import cli
from tag import index
from tag import locus
from tag import select
from tag import transcript
from gzip import open as gzopen
import sys

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


def open(filename, mode):
    if mode not in ['r', 'w']:
        raise ValueError('invalid mode "{}"'.format(mode))
    if filename in ['-', None]:  # pragma: no cover
        filehandle = sys.stdin if mode == 'r' else sys.stdout
        return filehandle
    openfunc = builtins.open
    if filename.endswith('.gz'):
        openfunc = gzopen
        mode += 't'
    return openfunc(filename, mode)
