#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the ISC license: see LICENSE.txt.
# -----------------------------------------------------------------------------


class Sequence():
    """Represents a biological sequence."""

    def __init__(self, defline, seq):
        assert defline.startswith('>') and defline[1] != ' '
        self.defline = defline
        self.seq = seq

    def __str__(self):
        return self.defline + '\n' + self.format_seq()

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.seq)

    def __lt__(self, other):
        """Rich comparison operator for Python 3 support."""
        if isinstance(other, Sequence):
            return self.seqid < other.seqid
        return False

    def __le__(self, other):
        """Rich comparison operator for Python 3 support."""
        if isinstance(other, Sequence):
            return self.seqid <= other.seqid
        return False

    def __gt__(self, other):
        """Rich comparison operator for Python 3 support."""
        if isinstance(other, Sequence):
            return self.seqid > other.seqid
        return True

    def __ge__(self, other):
        """Rich comparison operator for Python 3 support."""
        if isinstance(other, Sequence):
            return self.seqid >= other.seqid
        return True

    @property
    def seqid(self):
        return self.defline[1:].split(' ')[0]

    def format_seq(self, outstream=None, linewidth=70):
        """Print a sequence in a readable format."""
        if linewidth == 0 or len(self.seq) <= linewidth:
            if outstream is None:
                return self.seq
            else:
                print >> outstream, self.seq
                return

        i = 0
        seq = ''
        while i < len(self.seq):
            if outstream is None:
                seq += self.seq[i:i+linewidth] + '\n'
            else:
                print >> outstream, self.seq[i:i+linewidth]
            i += linewidth
        if outstream is None:
            return seq


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------

def test_repr():
    """[aeneas::Sequence] Test string representation of sequences."""
    import StringIO

    s1 = Sequence('>seq1', 'ACGT')
    assert '%s' % s1 == '>seq1\nACGT'
    assert '%r' % s1 == '>seq1\nACGT'

    s2 = Sequence('>contig2', 'AAAAACCCCCGGGGGNNNNNTTTTT')
    sio = StringIO.StringIO()
    s2.format_seq(linewidth=5, outstream=sio)
    assert sio.getvalue() == 'AAAAA\nCCCCC\nGGGGG\nNNNNN\nTTTTT\n'
    sio.close()

    s3 = Sequence('>scf3', ('A'*70) + ('C'*70))
    assert '%s' % s3 == '>scf3\n%s\n%s\n' % ('A'*70, 'C'*70)


def test_seqid():
    """[aeneas::Sequence] Test sequence ID parsing and retrieval."""
    s1 = Sequence('>scaffold_789', 'ACGT')
    assert s1.seqid == 'scaffold_789'

    s2 = Sequence('>gnl|aeneas|ABC123 [Bogus vulgaris]', 'ACGT')
    assert s2.seqid == 'gnl|aeneas|ABC123'

    try:
        s3 = Sequence('> no_space_allowed', 'ACGT')
    except AssertionError:
        pass


def test_len():
    """[aeneas::Sequence] Test length."""
    assert len(Sequence('>chr1', 'ACGT')) == 4
    assert len(Sequence('>contig2', 'AAAAACCCCCGGGGGNNNNNTTTTT')) == 25
    assert len(Sequence('>lcl|TheAccession', 'AACCGGNNTT')) == 10


def test_compare():
    """[aeneas::Sequence] Test sorting and comparison"""
    from .feature import Feature
    s1 = Sequence('>chr1', 'ACGT')
    s2 = Sequence('>contig2', 'AAAAACCCCCGGGGGNNNNNTTTTT')
    s3 = Sequence('>gnl|aeneas|ABC123 [Bogus vulgaris]', 'ACGT')

    gff3 = ['chr', 'vim', 'gene', '1000', '2000', '.', '+', '.', 'ID=g1']
    f1 = Feature('\t'.join(gff3))

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
