#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import tag


def primary_mrna(entrystream, parenttype='gene'):
    """
    Select a single mRNA as a representative for each protein-coding gene.

    The primary mRNA is the one with the longest translation product. In cases
    where multiple isoforms have the same translated length, the feature ID is
    used for sorting.

    >>> reader = tag.reader.GFF3Reader(tag.pkgdata('pdom-withseq.gff3'))
    >>> filter = tag.transcript.primary_mrna(reader)
    >>> for gene in tag.select.features(filter, type='gene'):
    ...    assert gene.num_children == 1
    """
    for entry in entrystream:
        if not isinstance(entry, tag.feature.Feature):
            yield entry
            continue

        for feature in tag.select.features(entry, parenttype, traverse=True):
            mrnas = [f for f in feature.children if f.type == 'mRNA']
            if len(mrnas) == 0:
                continue
            mrnas.sort(key=lambda m: (m.cdslen, m.get_attribute('ID')))
            mrnas.pop()
            feature.children = [c for c in feature.children if c not in mrnas]
        yield entry
