#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import tag


def primary(entrystream, parenttype='gene'):
    for entry in entrystream:
        if not isinstance(entry, tag.feature.Feature):
            yield entry
            continue

        for feature in tag.select.features(entry, parenttype, traverse=True):
            mrnas = [m for m in tag.select.features(feature.children, 'mRNA')]
            mrnas.sort(key=lambda m: (m.cdslen, m.get_attribute('ID')))
            mrnas.pop()
            feature.children = [c for c in feature.children if c not in mrnas]
        yield entry
