#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import tag


def features(entrystream, type=None, traverse=False):
    for feature in entry_type_filter(entrystream, tag.feature.Feature):
        if traverse:
            assert type
            if type == feature.type:
                yield feature
            else:
                for subfeature in feature:
                    if type == subfeature.type:
                        yield subfeature
        else:
            if not type or type == feature.type:
                yield feature


def entry_type_filter(entrystream, entryclass):
    for entry in entrystream:
        if isinstance(entry, entryclass):
            yield entry
