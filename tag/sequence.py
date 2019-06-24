#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import re


class Sequence(object):
    """
    Represents a sequence (almost always a nucleotide sequence).

    We do not encourage putting FASTA sequences in your GFF3 files, but since
    the specification explicitly allows it we have to handle it. :-(

    >>> s = Sequence('>contig1 description', 'ACGT')
    >>> s.defline
    '>contig1 description'
    >>> s.seq
    'ACGT'
    >>> s.seqid
    'contig1'
    >>> len(s)
    4
    >>> s = Sequence('>gi|12345|gb|BOGUSSEQ', 'GATTACA')
    >>> s.seq
    'GATTACA'
    >>> s.accession
    'BOGUSSEQ'
    """

    def __init__(self, defline, seq):
        assert defline.startswith('>') and defline[1] != ' '
        self.defline = defline
        self.seq = seq.strip()

    def __str__(self):
        return self.defline + '\n' + self.format_seq()

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.seq)

    def __lt__(self, other):
        if isinstance(other, Sequence):
            return self.seqid < other.seqid
        return False

    def __le__(self, other):
        if isinstance(other, Sequence):
            return self.seqid <= other.seqid
        return False

    def __gt__(self, other):
        if isinstance(other, Sequence):
            return self.seqid > other.seqid
        return True

    def __ge__(self, other):
        if isinstance(other, Sequence):
            return self.seqid >= other.seqid
        return True

    @property
    def seqid(self):
        return self.defline[1:].split(' ')[0]

    @property
    def accession(self):
        """
        Parse accession number from commonly supported formats.

        If the defline does not match one of the following formats, the entire
        description (sans leading caret) will be returned.

        * >gi|572257426|ref|XP_006607122.1|
        * >gnl|Tcas|XP_008191512.1
        * >lcl|PdomMRNAr1.2-10981.1
        """
        accession = None
        if self.defline.startswith('>gi|'):
            match = re.match(r'>gi\|\d+\|[^\|]+\|([^\|\n ]+)', self.defline)
            if match:
                accession = match.group(1)
        elif self.defline.startswith('>gnl|'):
            match = re.match(r'>gnl\|[^\|]+\|([^\|\n ]+)', self.defline)
            if match:
                accession = match.group(1)
        elif self.defline.startswith('>lcl|'):
            match = re.match(r'>lcl\|([^\|\n ]+)', self.defline)
            if match:
                accession = match.group(1)
        return accession

    def format_seq(self, outstream=None, linewidth=70):
        """
        Print a sequence in a readable format.

        :param outstream: if `None`, formatted sequence is returned as a
                          string; otherwise, it is treated as a file-like
                          object and the formatted sequence is printed to the
                          outstream
        :param linewidth: width for wrapping sequences over multiple lines; set
                          to 0 for no wrapping
        """
        if linewidth == 0 or len(self.seq) <= linewidth:
            if outstream is None:
                return self.seq
            else:
                print(self.seq, file=outstream)
                return

        i = 0
        seq = ''
        while i < len(self.seq):
            if outstream is None:
                seq += self.seq[i:i+linewidth] + '\n'
            else:
                print(self.seq[i:i+linewidth], file=outstream)
            i += linewidth
        if outstream is None:
            return seq
