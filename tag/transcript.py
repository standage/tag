#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
from collections import defaultdict
from sys import stderr
import tag


# One day I may write a program to crawl the Sequence Ontology and build a more
# comprehensive version of this set of terms. For now, this should handle >95%
# of GFF3 files in the wild.
type_terms = set([
    'mRNA',
    'transcript',
    'primary_transcript',
    'lnc_RNA',
    'tRNA',
    'rRNA',
])


def _emplace_pmrna(mrnas, parent, strict=False):
    """Retrieve the primary mRNA and discard all others."""
    mrnas.sort(key=lambda m: (m.cdslen, m.get_attribute('ID')))
    pmrna = mrnas.pop()
    if strict:
        parent.children = [pmrna]
    else:
        parent.children = [c for c in parent.children if c not in mrnas]


def _emplace_transcript(transcripts, parent):
    """Retrieve the primary transcript and discard all others."""
    transcripts.sort(key=lambda t: (len(t), t.get_attribute('ID')))
    pt = transcripts.pop()
    parent.children = [pt]


def primary_mrna(entrystream, parenttype='gene'):
    """
    Select a single mRNA as a representative for each protein-coding gene.

    The primary mRNA is the one with the longest translation product. In cases
    where multiple isoforms have the same translated length, the feature ID is
    used for sorting.

    This function **does not** return only mRNA features, it returns all GFF3
    entry types (pragmas, features, sequences, etc). The function **does**
    modify the gene features that pass through to ensure that they have at most
    a single mRNA feature.

    >>> reader = tag.GFF3Reader(tag.tests.data_stream('pdom-withseq.gff3'))
    >>> filter = tag.transcript.primary_mrna(reader)
    >>> for gene in tag.select.features(filter, type='gene'):
    ...    assert gene.num_children == 1
    """
    for entry in entrystream:
        if not isinstance(entry, tag.Feature):
            yield entry
            continue

        for parent in tag.select.features(entry, parenttype, traverse=True):
            mrnas = [f for f in parent.children if f.type == 'mRNA']
            if len(mrnas) == 0:
                continue
            _emplace_pmrna(mrnas, parent)
        yield entry


def _get_primary_type(ttypes, parent, logstream=stderr):
    """Check for multiple transcript types and, if possible, select one."""
    if len(ttypes) > 1:
        if logstream:  # pragma: no branch
            message = '[tag::transcript::primary_transcript]'
            message += ' WARNING: feature {:s}'.format(parent.slug)
            message += ' has multiple associated transcript types'
            message += ' {}'.format(ttypes)
            print(message, file=logstream)
        if 'mRNA' not in ttypes:
            message = (
                'cannot resolve multiple transcript types if "mRNA" is'
                ' not one of those types {}'.format(ttypes)
            )
            raise Exception(message)
        ttypes = ['mRNA']
    return ttypes[0]


def primary_transcript(entrystream, parenttype='gene', logstream=stderr):
    """
    Select a single transcript as a representative for each gene.

    This function is a generalization of the `primary_mrna` function that
    attempts, under certain conditions, to select a single transcript as a
    representative for each gene. If a gene encodes multiple transcript types,
    one of those types must be **mRNA** or the function will complain loudly
    and fail.

    For mRNAs, the primary transcript is selected according to translated
    length. For all other transcript types, the length of the transcript
    feature itself is used. I'd be eager to hear suggestions for alternative
    selection criteria.

    Like the `primary_mrna` function, this function **does not** return only
    transcript features. It **does** modify gene features to ensure that each
    has at most one transcript feature.

    >>> reader = tag.GFF3Reader(
    ...     tag.tests.data_stream('psyllid-mixed-gene.gff3.gz')
    ... )
    >>> gene_filter = tag.select.features(reader, type='gene')
    >>> trans_filter = tag.transcript.primary_transcript(gene_filter)
    >>> for gene in trans_filter:
    ...     assert gene.num_children == 1

    In cases where the direct children of a gene feature have heterogenous
    types, the `primary_mrna` function will only discard mRNA features. This
    function, however, will discard all direct children of the gene that are
    not the primary transcript, including non-transcript children. This is a
    retty subtle distinction, and anecdotal experience suggests that cases in
    which the distinction actually matters are extremely rare.
    """
    for entry in entrystream:
        if not isinstance(entry, tag.Feature):
            yield entry
            continue

        for parent in tag.select.features(entry, parenttype, traverse=True):
            if parent.num_children == 0:
                continue

            transcripts = defaultdict(list)
            for child in parent.children:
                if child.type in type_terms:
                    transcripts[child.type].append(child)

            if len(transcripts) == 0:
                continue

            ttypes = list(transcripts.keys())
            ttype = _get_primary_type(ttypes, parent)
            transcript_list = transcripts[ttype]

            if ttype == 'mRNA':
                _emplace_pmrna(transcript_list, parent, strict=True)
            else:
                _emplace_transcript(transcript_list, parent)

        yield entry
