#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2016 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import glob
import pytest
import tag
import tag.__main__
from tag.tests import data_file, data_stream
import sys


def test_cli_args(capsys):
    with pytest.raises(SystemExit) as se:
        args = tag.cli.parser().parse_args(['-v'])
    terminal = capsys.readouterr()
    assert tag.__version__ in terminal.out or tag.__version__ in terminal.err

    with pytest.raises(SystemExit) as se:
        args = tag.cli.parser().parse_args(['gff3', '-h'])
    terminal = capsys.readouterr()
    assert 'input file in GFF3 format' in terminal.out

    arglist = ['gff3', '-r', data_file('mito-trna.gff3')]
    args = tag.cli.parser().parse_args(arglist)
    tag.__main__.main(args)
    terminal = capsys.readouterr()
    testout = data_stream('mito-trna-out.gff3').read()
    assert terminal.out == testout

    with pytest.raises(SystemExit):
        tag.__main__.main()


def test_gff3_strict(capsys):
    arglist = ['gff3', data_file('mito-trna.gff3')]
    args = tag.cli.parser().parse_args(arglist)
    with pytest.raises(AssertionError) as ae:
        tag.cli.gff3.main(args)
    assert 'not contained within its span' in str(ae)
    terminal = capsys.readouterr()

    args.strict = False
    tag.cli.gff3.main(args)
    terminal = capsys.readouterr()
    testout = data_stream('mito-trna-out.gff3').read()
    assert terminal.out == testout


@pytest.mark.parametrize('gff3,ftype,expected_output', [
    ('oluc-20kb.gff3', 'CDS', '14100\n'),
    ('bogus-aligns.gff3', 'cDNA_match', '7006\n'),
    ('bogus-genes.gff3', 'gene', '18000\n'),
    ('bogus-genes.gff3', 'mRNA', '18000\n'),
    ('bogus-genes.gff3', 'exon', '11000\n'),
])
def test_occ(gff3, ftype, expected_output, capsys):
    arglist = ['occ', data_file(gff3), ftype]
    args = tag.cli.parser().parse_args(arglist)
    tag.cli.occ.main(args)
    terminal = capsys.readouterr()
    assert terminal.out == expected_output


def test_pmrna(capsys):
    arglist = ['pmrna', data_file('nanosplice.gff3')]
    args = tag.cli.parser().parse_args(arglist)
    tag.cli.pmrna.main(args)
    terminal = capsys.readouterr()
    exp_out = data_stream('nanosplice-primary.gff3').read()
    assert terminal.out.strip() == exp_out.strip()


def test_sum(capsys):
    infile = data_file('GCF_001639295.1_ASM163929v1_genomic.gff.gz')
    args = tag.cli.parser().parse_args(['sum', infile])
    tag.cli.sum.main(args)
    terminal = capsys.readouterr()
    out = terminal.out.strip().split('\n')[1:]
    exp_out = data_stream('sum-test-out.txt').read().strip().split('\n')[1:]
    assert out == exp_out


def test_merge(capsys):
    infiles = glob.glob(data_file('ex-red-?.gff3'))
    arglist = ['merge'] + infiles
    args = tag.cli.parser().parse_args(arglist)
    tag.cli.merge.main(args)
    terminal = capsys.readouterr()
    exp_out = data_stream('ex-red-merged.gff3').read()
    assert terminal.out.strip() == exp_out.strip()


def test_locuspocus(capsys):
    infiles = glob.glob(data_file('Ye.*.min.gff3.gz'))
    arglist = ['locuspocus'] + infiles
    args = tag.cli.parser().parse_args(arglist)
    tag.cli.locuspocus.main(args)
    terminal = capsys.readouterr()
    exp_out = data_stream('Ye.loci.gff3').read()
    assert terminal.out.strip() == exp_out.strip()


def test_pep2nuc(capsys):
    arglist = [
        'pep2nuc', '-k', 'protein', data_file('Ypes-abinit.gff3.gz'),
        data_file('Ypes-signalp-prot.gff3.gz')
    ]
    args = tag.cli.parser().parse_args(arglist)
    tag.cli.pep2nuc.main(args)
    terminal = capsys.readouterr()
    exp_out = data_stream('Ypes-signalp-nucl.gff3.gz').read()
    assert terminal.out.strip() == exp_out.strip()


def test_pep2nuc_nokeep(capsys):
    arglist = [
        'pep2nuc', data_file('Ypes-abinit.gff3.gz'),
        data_file('Ypes-signalp-prot.gff3.gz')
    ]
    args = tag.cli.parser().parse_args(arglist)
    tag.cli.pep2nuc.main(args)
    terminal = capsys.readouterr()
    exp_out = data_stream('Ypes-signalp-nucl-nk.gff3.gz').read()
    assert terminal.out.strip() == exp_out.strip()


def test_pep2nuc_missing_id(capsys):
    arglist = [
        'pep2nuc', data_file('Ypes-abinit-minus1.gff3.gz'),
        data_file('Ypes-signalp-prot.gff3.gz')
    ]
    args = tag.cli.parser().parse_args(arglist)
    tag.cli.pep2nuc.main(args)
    terminal = capsys.readouterr()
    msg = '[tag::pep2nuc] WARNING: protein identifier "cds000008" not defined'
    assert msg in terminal.err
