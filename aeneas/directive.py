#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the ISC license: see LICENSE.txt.
# -----------------------------------------------------------------------------
"""Represents a directive from a GFF3 file."""

import re
from region import Region


dirtypes = ['gff-version', 'sequence-region', 'feature-ontology', 'species',
            'attribute-ontology', 'source-ontology', 'genome-build']


class Directive():
    """
    This class is primarily for error checking and data access.

    Once created, `Directive` objects should be treated as read-only: modify at
    your peril! Also, separator directives (`###`) and the `##FASTA` directive
    are handled directly by parsers and not by this class.
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
            self.region = Region(int(formatmatch.group(2)),
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


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------

def test_basic():
    """[aeneas::Directive] Test basic object construction."""
    d = Directive('##gff-version 3')
    assert d.type == 'gff-version' and d.version == '3'
    d = Directive('##gff-version    3')  # 3 spaces
    assert d.type == 'gff-version' and d.version == '3'
    d = Directive('##gff-version	3')  # tab
    assert d.type == 'gff-version' and d.version == '3'

    try:
        d = Directive('')  # No data
    except AssertionError:
        pass
    try:
        d = Directive('##gff-version   2.2')  # Only version 3 supported
    except AssertionError:
        pass
    try:
        d = Directive('not a directive')
    except AssertionError:
        pass
    try:
        d = Directive('# Still not a directive')
    except AssertionError:
        pass


def test_custom_directive():
    """[aeneas::Directive] Test custom directive type."""
    d1 = Directive('##bogus-directive')
    d2 = Directive('##bonus-directive   abc 1 2 3')
    d3 = Directive('##Type DNA NC_005213.1')

    assert d1.type is None and d1.dirtype == 'bogus-directive'
    assert d1.data is None
    assert d2.type is None and d2.dirtype == 'bonus-directive'
    assert d2.data == 'abc 1 2 3'
    assert d3.type is None and d3.dirtype == 'Type'
    assert d3.data == 'DNA NC_005213.1'


def test_sequence_region():
    """[aeneas::Directive] Test sequence-region directive type."""
    r1 = Directive('##sequence-region ctg123 1 1497228')
    r2 = Directive('##sequence-region   ctg123 1 1497228')  # 3 spaces
    r3 = Directive('##sequence-region	ctg123 1 1497228')  # tab
    r4 = Directive('##sequence-region 1 1 1000')

    assert r1.type == 'sequence-region' and r1.seqid == 'ctg123' and \
        r1.region == Region(1, 1497228)
    assert r2.type == 'sequence-region' and r2.seqid == 'ctg123' and \
        r2.region == Region(1, 1497228)
    assert r3.type == 'sequence-region' and r3.seqid == 'ctg123' and \
        r3.region == Region(1, 1497228)
    assert r4.type == 'sequence-region' and r4.seqid == '1' and \
        r4.region == Region(1, 1000)

    try:
        # Invalid region
        r5 = Directive('##sequence-region   BoGuScHr 123456 4321')
    except AssertionError:
        pass


def test_ontology_directives():
    """[aeneas::Directive] Test ontology directives."""
    so_uri = ('http://song.cvs.sourceforge.net/viewvc/song/ontology/'
              'so.obo?revision=1.263')
    attr_uri = 'http://www.bogus.edu/attr-o.obo'
    src_uri = 'http://www.bogus.edu/src-o.obo'

    o1 = Directive('##feature-ontology ' + so_uri)
    o2 = Directive('##attribute-ontology   ' + attr_uri)
    o3 = Directive('##source-ontology	' + src_uri)

    assert o1.type == 'feature-ontology' and o1.uri == so_uri
    assert o2.type == 'attribute-ontology' and o2.uri == attr_uri
    assert o3.type == 'source-ontology' and o3.uri == src_uri


def test_species_directive():
    """[aeneas::Directive] Test species directive."""
    amel = 'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=7460'
    pdom = 'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=743375'
    sinv = 'http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=13686'

    s1 = Directive('##species ' + amel)
    s2 = Directive('##species   ' + pdom)
    s3 = Directive('##species	' + sinv)

    assert s1.type == 'species' and s1.uri == amel
    assert s2.type == 'species' and s2.uri == pdom
    assert s3.type == 'species' and s3.uri == sinv


def test_genome_build_directive():
    """[aeneas::Directive] Test genome-build directive."""
    b1 = Directive('##genome-build NCBI B36')
    b2 = Directive('##genome-build   WormBase ws110')
    b3 = Directive('##genome-build	FlyBase  r4.1')

    assert b1.type == 'genome-build' and b1.source == 'NCBI' and \
        b1.build_name == 'B36'
    assert b2.type == 'genome-build' and b2.source == 'WormBase' and \
        b2.build_name == 'ws110'
    assert b3.type == 'genome-build' and b3.source == 'FlyBase' and \
        b3.build_name == 'r4.1'
