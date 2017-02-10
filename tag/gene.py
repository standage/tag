#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2017 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from collections import defaultdict
import tag


def exon_dedup(entrystream, parenttype='gene'):
    for entry in entrystream:
        if not isinstance(entry, tag.feature.Feature):
            yield entry
            continue

        exons_by_location = defaultdict(list)
        feature_by_id = dict()
        for feature in tag.select.features(entry, parenttype, traverse=True):
            for subfeature in feature:
                fid = subfeature.get_attribute('ID')
                if fid is not None:
                    feature_by_id[fid] = subfeature

                if subfeature.type == 'exon':
                    exons_by_location[subfeature._range].append(subfeature)

            # resolve
        yield entry
