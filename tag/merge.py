#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2019 Battelle National Biodefense Institute.
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import heapq
import itertools
import tag


def merge(*sorted_streams):
    heap = list()
    heapq.heapify(heap)
    streams = itertools.chain(*sorted_streams)
    for record in heapq.merge(heap, *streams):
        yield record
