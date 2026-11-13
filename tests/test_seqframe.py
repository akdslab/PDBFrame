"""Smoke tests for SeqFrame on small hardcoded input."""

import pandas as pd
import pytest

from pdbframe import SeqFrame

P1 = 'MKVLWAALLVTFLAGCQAK'
P2 = 'XXMKVLWAALLVTFLAGCQAKXX'


def test_single_sequence_table():
    s = SeqFrame('mkv', name='a')
    assert s.df['a'].to_list() == ['M', 'K', 'V']
    assert s.df.index.to_list() == [1, 2, 3]
    assert s.statistics['length'] == 3


def test_both_call_styles_agree():
    a = SeqFrame('ATCGATCGATCG').align('ATCGATGGATCG')
    b = SeqFrame.align('ATCGATCGATCG', 'ATCGATGGATCG')
    pd.testing.assert_frame_equal(a.df, b.df)


def test_alignment_table_and_statistics():
    out = SeqFrame.align('ATCGATCGATCG', 'ATCGATGGATCG')
    assert out.df.columns.to_list() == ['seq1', 'seq2', 'seq1_pos', 'seq2_pos', 'match']
    assert out.df.loc[7, 'match'] == 'mismatch'
    st = out.statistics
    assert (st['matches'], st['mismatches'], st['gap_columns']) == (11, 1, 0)
    assert st['identity'] == pytest.approx(11 / 12)


def test_gaps_have_na_positions():
    out = SeqFrame.align('MKVAAL', 'MKVL', names=('long', 'short'))
    assert out.df['short_pos'].isna().sum() == 2
    assert out.statistics['short_gaps'] == 2
    assert out.sequences == {'long': 'MKVAAL', 'short': 'MKVL'}


def test_similarity_never_exceeds_one():
    # The source counted identical residues twice in "similarity"
    assert SeqFrame.align('AAAA', 'AAAA').statistics['similarity'] == 1.0


def test_local_alignment_keeps_aligned_region():
    st = SeqFrame.align(P1, P2, mode='local').statistics
    assert (st['seq2_start'], st['seq2_end']) == (3, 21)
    assert st['coverage'] == 1.0


def test_names_from_seqframes():
    out = SeqFrame(P1, name='hla').align(SeqFrame(P2, name='ref'))
    assert out.names == ['hla', 'ref']


def test_show_one_line(capsys):
    SeqFrame.align('MKVAAL', 'MKVL').show(one_line=True)
    lines = [l for l in capsys.readouterr().out.splitlines() if l]
    assert len(lines) == 3


def test_drop_gaps():
    out = SeqFrame.align('MKVAAL', 'MKVL', names=('long', 'short'))
    assert ''.join(out.drop_gaps().df['long']) == 'MKVL'
    assert len(out.drop_gaps('long')) == 6          # long has no gaps: nothing dropped
    assert out.drop_gaps().statistics['score'] is None


def test_multiple_alignment():
    seqs = ['CASSIRSSYEQYF', 'CASSIRSAYEQYF', 'CASSLGQAYEQYF', 'CASSIRSSGYEQYF']
    out = SeqFrame.align(*seqs)
    assert out.names == ['seq1', 'seq2', 'seq3', 'seq4']
    assert len({len(s) for s in out.aligned.values()}) == 1
    assert list(out.sequences.values()) == seqs


def test_multiple_alignment_keeps_each_pairwise_alignment():
    """Star merge: ref vs seq_i, read off the MSA, is the pairwise alignment."""
    ref, others = 'GSHSMRYFYTSVSRPG', ['GSHSMRYFFTSVSRPGRG', 'HSMRYFYTAVSRPG', 'GSHSWMRYFYTSVSRPG']
    msa = SeqFrame.align(ref, *others)
    for i, o in enumerate(others, 2):
        pair = SeqFrame.align(ref, o).aligned
        sub = msa.df[['seq1', f'seq{i}']]
        sub = sub[(sub != '-').any(axis=1)]
        assert ''.join(sub['seq1']) == pair['seq1']
        assert ''.join(sub[f'seq{i}']) == pair['seq2']


def test_local_multiple_rejected():
    with pytest.raises(ValueError):
        SeqFrame.align('MKV', 'MKA', 'MKL', mode='local')


CDR3 = ['CASSIRSSYEQYF', 'CASSIRSAYEQYF', 'CASSLGQAYEQYF', 'CASSIRSSGYEQYF']
CDR3_NAMES = ['tcr1', 'tcr2', 'tcr3', 'tcr4']


def test_pairwise_against_first():
    aln = SeqFrame.align(*CDR3, names=CDR3_NAMES)
    pairs = aln.pairwise()
    assert list(pairs) == ['tcr2', 'tcr3', 'tcr4']
    for n, p in pairs.items():
        assert p.names == ['tcr1', n]
        assert p.sequences == {'tcr1': CDR3[0], n: CDR3[CDR3_NAMES.index(n)]}
    # tcr2 has no insertion, so the column tcr4's G needed is gone
    assert len(pairs['tcr2']) == 13 and pairs['tcr2'].statistics['matches'] == 12
    assert pairs['tcr4'].statistics['gap_columns'] == 1
    assert pairs['tcr2'].score is None


def test_pairwise_other_reference_and_two_sequences():
    aln = SeqFrame.align(*CDR3, names=CDR3_NAMES)
    assert list(aln.pairwise('tcr3')) == ['tcr1', 'tcr2', 'tcr4']
    two = SeqFrame.align(P1, P2)
    (p,) = two.pairwise().values()
    assert p.aligned == two.aligned and p.score == two.score
    with pytest.raises(KeyError):
        aln.pairwise('nope')
    with pytest.raises(ValueError):
        SeqFrame(P1).pairwise()


def test_show_against_reference(capsys):
    aln = SeqFrame.align(*CDR3, names=CDR3_NAMES)
    aln.show(reference='tcr1')
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].split()[2] == 'CASSIRSS-YEQYF'
    assert lines[1].split()[2] == '.......A-.....'
    assert lines[3].split()[2] == '........G.....'
    aln.show(reference=True)
    assert capsys.readouterr().out.splitlines() == lines
    aln.show()
    assert capsys.readouterr().out.splitlines()[1].split()[2] == 'CASSIRSA-YEQYF'
