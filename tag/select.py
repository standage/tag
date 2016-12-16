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
    """
    Pull features out of the specified entry stream.

    :param entrystream: a stream of entries
    :param type: retrieve only features of the specified type; set to
                 :code:`None` to retrieve all features
    :param traverse: by default, only top-level features are selected; set
                     to :code:`True` to search each feature graph for the
                     specified feature type
    """
    for feature in entry_type_filter(entrystream, tag.feature.Feature):
        if traverse:
            if type is None:
                message = 'cannot traverse without a specific feature type'
                raise ValueError(message)
            if type == feature.type:
                yield feature
            else:
                for subfeature in feature:
                    if type == subfeature.type:
                        yield subfeature
        else:
            if not type or type == feature.type:
                yield feature


def directives(entrystream, type=None):
    """
    Pull directives out of the specified entry stream.

    :param entrystream: a stream of entries
    :param type: retrieve only directives of the specified type; set to
                 :code:`None` to retrieve all directives
    """
    for directive in entry_type_filter(entrystream, tag.directive.Directive):
        if not type or type == directive.type:
            yield directive


def sequences(entrystream):
    """Pull sequences out of the specified entry stream."""
    for sequence in entry_type_filter(entrystream, tag.sequence.Sequence):
        yield sequence


def entry_type_filter(entrystream, entryclass):
    """
    Generic entry filter.

    :param entrystream: a stream of entries
    :param entryclass: specify the type of entry upon which to filter
    """
    for entry in entrystream:
        if isinstance(entry, entryclass):
            yield entry
