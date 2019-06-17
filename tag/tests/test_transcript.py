#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import pytest
import tag
from tag.tests import data_file, data_stream
from tag.transcript import primary_mrna, primary_transcript


def test_primary_mrna():
    reader = tag.GFF3Reader(data_stream('nanosplice.gff3'))
    gene = next(tag.select.features(primary_mrna(reader), type='gene'))
    assert gene.cdslen is None
    assert gene.num_children == 1
    assert gene.children[0].get_attribute('ID') == 'mRNAsecond'

    reader = tag.GFF3Reader(data_stream('pdom-withseq.gff3'))
    for gene in tag.select.features(primary_mrna(reader), type='gene'):
        assert gene.num_children == 1

    reader = tag.GFF3Reader(data_stream('psyllid-100k.gff3'))
    for gene in tag.select.features(primary_mrna(reader), type='gene'):
        mrnas = [f for f in gene if f.type == 'mRNA']
        assert len(mrnas) <= 1, mrnas


def test_primary_transcript_mixed():
    reader = tag.GFF3Reader(data_stream('psyllid-mixed-gene.gff3.gz'))
    trans_filter = primary_transcript(reader)
    gene_filter = tag.select.features(trans_filter, type='gene')
    gene = next(gene_filter)

    assert gene.num_children == 1
    assert gene.children[0].type == 'mRNA'


def test_primary_transcript_mixed_unresolvable():
    reader = tag.GFF3Reader(data_stream('unresolvable-mixed.gff3'))
    trans_filter = primary_transcript(reader)
    gene_filter = tag.select.features(trans_filter, type='gene')
    with pytest.raises(Exception) as e:
        genes = list(gene_filter)
    assert 'cannot resolve multiple transcript types' in str(e)


def test_primary_transcript_multi():
    reader = tag.GFF3Reader(data_stream('psyllid-multi-trans.gff3.gz'))
    trans_filter = primary_transcript(reader)
    gene_filter = tag.select.features(trans_filter, type='gene')
    gene = next(gene_filter)

    assert gene.num_children == 1
    assert gene.children[0].type == 'transcript'
    assert gene.children[0].get_attribute('Name') == 'XR_541450.2'


def test_primary_transcript_boring_genes():
    reader = tag.GFF3Reader(data_stream('boring-genes.gff3'))
    trans_filter = primary_transcript(reader)
    gene_filter = tag.select.features(trans_filter, type='gene')
    for gene in gene_filter:
        assert gene.num_children == 0


def test_primary_transcript_boring_exons():
    reader = tag.GFF3Reader(data_stream('boring-exons.gff3'))
    trans_filter = primary_transcript(reader)
    gene_filter = tag.select.features(trans_filter, type='gene')
    for gene in gene_filter:
        t = [c for c in gene.children if c.type in tag.transcript.type_terms]
        assert len(t) == 0
