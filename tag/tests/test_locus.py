#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2019 Battelle National Biodefense Institute.
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import tag
from tag.tests import data_file


@pytest.fixture
def ye_in():
    return [
        tag.GFF3Reader(
            infilename=data_file('Ye.callgenes.min.gff3.gz'), assumesorted=True
        ),
        tag.GFF3Reader(
            infilename=data_file('Ye.glimmer.min.gff3.gz'), assumesorted=True
        ),
        tag.GFF3Reader(
            infilename=data_file('Ye.prodigal.min.gff3.gz'), assumesorted=True
        ),
    ]


def test_locus_parsing_basic(ye_in):
    test_ranges = (
        ('NC_008791.1', tag.Range(591, 1287)),
        ('NC_008791.1', tag.Range(1283, 1550)),
        ('NC_008791.1', tag.Range(1598, 2147)),
        ('NC_008791.1', tag.Range(2290, 2839)),
        ('NC_008791.1', tag.Range(3351, 4320)),
        ('NC_008800.1', tag.Range(269, 710)),
        ('NC_008800.1', tag.Range(801, 1263)),
        ('NC_008800.1', tag.Range(1448, 2441)),
        ('NC_008800.1', tag.Range(2594, 3266)),
        ('NC_008800.1', tag.Range(3275, 4742)),
    )
    compare_stream = zip(
        tag.locus.loci(*ye_in),
        test_ranges
    )
    for (seqid, rng, features), (testseqid, testrng) in compare_stream:
        assert seqid == seqid
        assert rng == testrng


def test_locus_parsing_strict_overlap(ye_in):
    test_ranges = (
        ('NC_008791.1', tag.Range(591, 1550)),
        ('NC_008791.1', tag.Range(1598, 2147)),
        ('NC_008791.1', tag.Range(2290, 2839)),
        ('NC_008791.1', tag.Range(3351, 4320)),
        ('NC_008800.1', tag.Range(269, 710)),
        ('NC_008800.1', tag.Range(801, 1263)),
        ('NC_008800.1', tag.Range(1448, 2441)),
        ('NC_008800.1', tag.Range(2594, 3266)),
        ('NC_008800.1', tag.Range(3275, 4742)),
    )
    compare_stream = zip(
        tag.locus.loci(*ye_in, minbp=1, minperc=0.0),
        test_ranges
    )
    for (seqid, rng, features), (testseqid, testrng) in compare_stream:
        assert seqid == seqid
        assert rng == testrng


def test_locus_parsing_ye():
    instreams = [
        tag.GFF3Reader(
            infilename=data_file('Ye.callgenes.gff3.gz'), assumesorted=True
        ),
        tag.GFF3Reader(
            infilename=data_file('Ye.glimmer.gff3.gz'), assumesorted=True
        ),
        tag.GFF3Reader(
            infilename=data_file('Ye.prodigal.gff3.gz'), assumesorted=True
        ),
    ]
    test_ranges = list()
    with tag.open(data_file('Ye.loci.tsv.gz'), 'r') as fh:
        for line in fh:
            seqid, start, end = line.strip().split('\t')
            test_ranges.append((seqid, tag.Range(int(start), int(end))))
    compare_stream = zip(
        tag.locus.loci(*instreams),
        test_ranges
    )
    for (seqid, rng, features), (testseqid, testrng) in compare_stream:
        assert seqid == seqid
        assert rng == testrng


@pytest.mark.parametrize('ftype,minperc,numloci', [
    (None, 0.0, 11),
    (None, 0.25, 12),
    ('gene', 0.0, 11),
    ('gene', 0.25, 12),
    ('cDNA_match', 0.0, 1),
    ('cDNA_match', 0.25, 1),
])
def test_apis_mellifera(ftype, minperc, numloci):
    instream = tag.GFF3Reader(infilename=data_file('honeybee-100kb.gff3.gz'))
    locusstream = tag.locus.loci(instream, featuretype=ftype, minperc=minperc)
    ranges = [r for s, r, f in locusstream]
    assert len(ranges) == numloci
