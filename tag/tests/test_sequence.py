#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
from tag import Feature
from tag import Sequence


def test_repr():
    """Test string representation of sequences."""
    try:
        from StringIO import StringIO
    except ImportError:  # pragma: no cover
        from io import StringIO

    s1 = Sequence('>seq1', 'ACGT')
    assert str(s1) == '>seq1\nACGT'
    assert repr(s1) == '>seq1\nACGT'
    sio = StringIO()
    s1.format_seq(linewidth=0, outstream=sio)
    assert sio.getvalue() == 'ACGT\n'
    sio.close()

    s2 = Sequence('>contig2', 'AAAAACCCCCGGGGGNNNNNTTTTT')
    sio = StringIO()
    s2.format_seq(linewidth=5, outstream=sio)
    assert sio.getvalue() == 'AAAAA\nCCCCC\nGGGGG\nNNNNN\nTTTTT\n'
    sio.close()

    s3 = Sequence('>scf3', ('A'*70) + ('C'*70))
    assert str(s3) == '>scf3\n{}\n{}\n'.format('A'*70, 'C'*70)


def test_seqid():
    """Test sequence ID parsing and retrieval."""
    s1 = Sequence('>scaffold_789', 'ACGT')
    assert s1.seqid == 'scaffold_789'

    s2 = Sequence('>gnl|aeneas|ABC123 [Bogus vulgaris]', 'ACGT')
    assert s2.seqid == 'gnl|aeneas|ABC123'

    with pytest.raises(AssertionError):
        s3 = Sequence('> no_space_allowed', 'ACGT')


def test_len():
    """Test length."""
    assert len(Sequence('>chr1', 'ACGT')) == 4
    assert len(Sequence('>contig2', 'AAAAACCCCCGGGGGNNNNNTTTTT')) == 25
    assert len(Sequence('>lcl|TheAccession', 'AACCGGNNTT')) == 10


def test_compare():
    """Test sorting and comparison"""
    s1 = Sequence('>chr1', 'ACGT')
    s2 = Sequence('>contig2', 'AAAAACCCCCGGGGGNNNNNTTTTT')
    s3 = Sequence('>gnl|aeneas|ABC123 [Bogus vulgaris]', 'ACGT')
    f1 = Feature('chr', 'gene', 1000, 2000, strand='+')

    assert s1 < s2
    assert s2 > s1
    assert s1 <= s3
    assert s3 >= s2
    assert s1 > f1
    assert s2 >= f1
    assert not s3 < f1
    assert not s3 <= f1
    assert sorted([s1, s2, s3, f1]) == [f1, s1, s2, s3], \
        sorted([s1, s2, s3, f1])


def test_accession():
    """Test accession parsing."""
    s1 = Sequence('>gi|572257426|ref|XP_006607122.1|', 'ACGT')
    s2 = Sequence('>gnl|Tcas|XP_008191512.1', 'ACGT')
    s3 = Sequence('>lcl|PdomMRNAr1.2-10981.1 some description here', 'ACGT')
    s4 = Sequence('>not|standard|format', 'ACGT')
    s5 = Sequence('>rawid', 'ACGT')
    s6 = Sequence('>gi|bogus|fmt|doesntparse', 'ACGT')
    s7 = Sequence('>gnl|bogus', 'ACGT')
    s8 = Sequence('>lcl|', 'ACGT')

    assert s1.accession == 'XP_006607122.1'
    assert s2.accession == 'XP_008191512.1'
    assert s3.accession == 'PdomMRNAr1.2-10981.1' and \
        s3.seqid == 'lcl|PdomMRNAr1.2-10981.1'
    for seq in [s4, s5, s6, s7, s8]:
        assert seq.accession is None

    with pytest.raises(AssertionError):
        s = Sequence('gi|572257426|ref|XP_006607122.1|', 'ACGT')
