#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the ISC license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from __future__ import print_function


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
                print(self.seq, file=outstream)
                return

        i = 0
        while i < len(self.seq):
            if outstream is None:
                return self.seq[i:i+linewidth]
            else:
                print(self.seq[i:i+linewidth], file=outstream)
            i += linewidth


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------

def test_repr():
    """[aeneas::Sequence] Test string representation of sequences."""
    try:
        from StringIO import StringIO
    except ImportError:
        from io import StringIO

    s1 = Sequence('>seq1', 'ACGT')
    assert '%s' % s1 == '>seq1\nACGT'
    assert '%r' % s1 == '>seq1\nACGT'

    s2 = Sequence('>contig2', 'AAAAACCCCCGGGGGNNNNNTTTTT')
    sio = StringIO()
    s2.format_seq(linewidth=5, outstream=sio)
    assert sio.getvalue() == 'AAAAA\nCCCCC\nGGGGG\nNNNNN\nTTTTT\n'
    sio.close()


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
