#!/usr/bin/env python
#
# ------------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# ------------------------------------------------------------------------------

import re
from tag.range import Range


dirtypes = ['gff-version', 'sequence-region', 'feature-ontology', 'species',
            'attribute-ontology', 'source-ontology', 'genome-build']


class Directive(object):
    """
    Represents a directive from a GFF3 file.

    This class is primarily for error checking and data access. Once created,
    `Directive` objects should be treated as read-only: modify at your peril!
    Also, separator directives (`###`) and the `##FASTA` directive are handled
    directly by parsers and not by this class.

    Directives not explicitly declared in the GFF3 spec are application
    specific: they will be parsed without complaint, but no guarantees can be
    made about accessing their attributes.

    >>> sr = Directive('##sequence-region chr1 5000 10000')
    >>> sr.type
    'sequence-region'
    >>> sr.seqid
    'chr1'
    >>> gb = Directive('##genome-build   BeeBase 4.5')
    >>> gb.type
    'genome-build'
    >>> gb.source
    'BeeBase'
    """

    def __init__(self, data):
        assert data.startswith('##')
        self._rawdata = data

        formatmatch = re.match('##gff-version\s+(\d+)', data)
        if formatmatch:
            self.dirtype = 'gff-version'
            self.version = formatmatch.group(1)
            assert self.version == '3', 'Only GFF version 3 is supported'
            return

        formatmatch = re.match('##sequence-region\s+(\S+) (\d+) (\d+)', data)
        if formatmatch:
            self.dirtype = 'sequence-region'
            self.seqid = formatmatch.group(1)
            self.range = Range(int(formatmatch.group(2)) - 1,
                               int(formatmatch.group(3)))
            return

        formatmatch = re.match('##((feature|attribute|source)-ontology)'
                               '\s+(\S+)', data)
        if formatmatch:
            self.dirtype = formatmatch.group(1)
            self.uri = formatmatch.group(3)
            return

        formatmatch = re.match('##species\s+(\S+)', data)
        if formatmatch:
            self.dirtype = 'species'
            self.uri = formatmatch.group(1)
            return

        formatmatch = re.match('##genome-build\s+(\S+)\s+(\S+)', data)
        if formatmatch:
            self.dirtype = 'genome-build'
            self.source = formatmatch.group(1)
            self.build_name = formatmatch.group(2)
            return

        formatmatch = re.match('##(\S+)(\s+(.+))*', data)
        assert formatmatch
        self.dirtype = formatmatch.group(1)
        self.data = formatmatch.group(3)

        assert self.dirtype is not None

    @property
    def type(self):
        if self.dirtype in dirtypes:
            return self.dirtype
        return None

    def __repr__(self):
        return self._rawdata

    def __lt__(self, other):
        if self.type == 'gff-version':
            return True

        if self.type == 'sequence-region':
            if isinstance(other, Directive):
                if other.type == 'gff-version':
                    return False
                elif other.type == 'sequence-region':
                    if self.seqid == other.seqid:
                        return self.range.__lt__(other.range)
                    else:
                        return self.seqid < other.seqid
                else:
                    return True
            else:
                return True

        if isinstance(other, Directive):
            return self._rawdata < other._rawdata
        else:
            return True

    def __le__(self, other):
        if self.type == 'gff-version':
            return True

        if self.type == 'sequence-region':
            if isinstance(other, Directive):
                if other.type == 'gff-version':
                    return False
                elif other.type == 'sequence-region':
                    if self.seqid == other.seqid:
                        return self.range.__le__(other.range)
                    else:
                        return self.seqid <= other.seqid
                else:
                    return True
            else:
                return True

        if isinstance(other, Directive):
            return self._rawdata <= other._rawdata
        else:
            return True

    def __gt__(self, other):
        return not self.__le__(other)

    def __ge__(self, other):
        return not self.__lt__(other)
