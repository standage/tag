#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2019 Battelle National Biodefense Institute.
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import os
from pkg_resources import resource_filename
import tag


def data_file(path):
    pathparts = path.split('/')
    relpath = os.path.join('tests', 'data', *pathparts)
    return resource_filename('tag', relpath)


def data_stream(path):
    return tag.open(data_file(path), 'r')
