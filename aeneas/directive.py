#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# -----------------------------------------------------------------------------

import re
from .range import Range


dirtypes = ['gff-version', 'sequence-region', 'feature-ontology', 'species',
            'attribute-ontology', 'source-ontology', 'genome-build']


class Directive():
    """
    Represents a directive from a GFF3 file.

    This class is primarily for error checking and data access. Once created,
    `Directive` objects should be treated as read-only: modify at your peril!
    Also, separator directives (`###`) and the `##FASTA` directive are handled
    directly by parsers and not by this class.
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
            self.region = Range(int(formatmatch.group(2)) - 1,
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
        """
        Directives not following one of the explicitly described formats in the
        GFF3 spec are application specific and not supported.
        """
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
                        return self.region.__lt__(other.region)
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
                        return self.region.__le__(other.region)
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
