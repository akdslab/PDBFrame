"""Tests for PDBFrame on small structures in tests/data.

3tfk_rcsb.pdb / 3tfk_af2.pdb are an experimental structure and an AlphaFold2 model
of the same 4-chain complex; chain A is numbered one higher in the model. The
pinned metric values were checked against independent implementations: DockQ
against the official DockQ program (v2), lDDT and the no-superposition TM-score
against an independent implementation.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pdbframe import PDBFrame, SeqFrame
from pdbframe.structure.geometry import kabsch

DATA = Path(__file__).parent / 'data'


@pytest.fixture(scope='module')
def ref():
    return PDBFrame.read(DATA / '3tfk_rcsb.pdb')


@pytest.fixture(scope='module')
def model():
    return PDBFrame.read(DATA / '3tfk_af2.pdb')


def _rotated(pf, seed=0):
    rng = np.random.default_rng(seed)
    q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    q *= np.sign(np.linalg.det(q))           # a proper rotation
    return pf.transform(q, rng.normal(size=3) * 10)


# -------------------------------------------------------------------- reading

def test_read_pdb(ref):
    assert ref.chains == ['A', 'C', 'D', 'E']
    assert len(ref) == 4972
    assert {c: len(s) for c, s in ref.sequences.items()} == {'A': 174, 'C': 9, 'D': 201, 'E': 239}
    assert ref.sequences['C'] == 'QLSDVPMDL'
    assert ref.df['observed_index'].min() == 1


def test_read_cif():
    pf = PDBFrame.read(DATA / 'aurka_tpx2_af3.cif')
    assert pf.chains == ['A', 'B'] and len(pf) == 2585


@pytest.mark.parametrize('suffix', ['.pdb', '.cif', '.csv'])
def test_write_read_round_trip(ref, tmp_path, suffix):
    back = PDBFrame.read(ref.write(tmp_path / f'x{suffix}'))
    assert back.sequences == ref.sequences
    np.testing.assert_allclose(back.coords, ref.coords, atol=1e-3)


def test_unknown_format(tmp_path):
    with pytest.raises(ValueError, match='unknown structure format'):
        PDBFrame.read(tmp_path / 'x.xyz')


# ------------------------------------------------------------------ selecting

def test_select(ref):
    assert ref.select(chains='C').chains == ['C']
    assert len(ref.select(atoms='ca')) == 623
    assert set(ref.select(atoms='backbone').df['atom_name']) == {'N', 'CA', 'C', 'O'}
    stretch = ref.select(chains='A', residues=(10, 19))
    assert stretch.residues['residue_number'].tolist() == list(range(10, 20))
    assert stretch.df['observed_index'].min() == 1          # recomputed per frame
    assert len(ref.select(chains='C', invert=True)) == len(ref) - 70


def test_residues_and_seqframe(ref):
    res = ref.residues
    assert len(res) == 623 and res['atoms'].sum() == len(ref)
    sf = ref.seqframe('C')
    assert isinstance(sf, SeqFrame) and sf.names == ['C']


def test_seqframes_all_and_some_chains(ref):
    sfs = ref.seqframes()
    assert list(sfs) == ref.chains
    assert all(sf.sequences[c] == ref.sequences[c] for c, sf in sfs.items())
    assert list(ref.seqframes('C')) == ['C']
    assert list(ref.seqframes(['D', 'C'])) == ['D', 'C']
    with pytest.raises(KeyError):
        ref.seqframes(['Z'])


def test_contacts_and_near(ref):
    peptide = ref.select(chains='C')
    shell = ref.near(peptide, cutoff=5)
    assert 'C' in shell.chains and len(shell.residues) == 55
    pairs = ref.contacts(cutoff=4)
    assert (pairs['chain_1'] != pairs['chain_2']).all()
    assert (pairs['distance'] <= 4).all()


# -------------------------------------------------------------------- editing

def test_rename_renumber_add(ref):
    pep = ref.select(chains='C').rename_chains({'C': 'P'})
    assert pep.chains == ['P']
    assert pep.renumber(100).residues['residue_number'].tolist() == list(range(100, 109))
    combined = ref.select(chains='A') + pep
    assert combined.chains == ['A', 'P']
    with pytest.raises(ValueError, match='rename_chains'):
        ref + ref


# ------------------------------------------------------------------- matching

def test_match_handles_offset_numbering(ref, model):
    mine, theirs = model.match(ref)
    assert len(mine) == len(theirs) == 4972
    assert (mine.df['atom_name'].to_numpy() == theirs.df['atom_name'].to_numpy()).all()
    offset = mine.df['residue_number'] - theirs.df['residue_number']
    assert set(offset[mine.df['chain'] == 'A']) == {1}
    # by number, chain A pairs the wrong residues, so far fewer atoms match by name
    assert len(model.match(ref, by='number')[0]) < len(mine)


def test_match_explicit_chain_pairs(ref):
    renamed = ref.rename_chains({'A': 'X'})
    mine, theirs = renamed.match(ref, chains={'X': 'A'})
    assert mine.chains == ['X'] and theirs.chains == ['A']


# -------------------------------------------------------------- chain mapping

RENAME = {'A': 'H', 'C': 'P', 'D': 'X', 'E': 'Y'}


def test_map_chains_by_sequence(ref, model):
    shuffled = model.rename_chains(RENAME)
    assert shuffled.map_chains(ref) == {v: k for k, v in RENAME.items()}


def test_metrics_do_not_need_matching_chain_ids(ref, model):
    shuffled = model.rename_chains(RENAME)
    assert shuffled.metrics.lddt(ref) == pytest.approx(model.metrics.lddt(ref))
    assert shuffled.metrics.rmsd(ref) == pytest.approx(model.metrics.rmsd(ref))
    # DockQ sides are named with the reference's chain IDs
    assert shuffled.metrics.dockq(ref, 'A', 'C')['dockq'] == pytest.approx(0.914, abs=1e-3)
    # 'same' pairs by ID only, so nothing pairs up
    with pytest.raises(ValueError, match='no atoms in common'):
        shuffled.match(ref, chains='same')


def test_rename_chains_like_reference(ref, model):
    back = model.rename_chains(RENAME).rename_chains(ref)
    assert back.chains == ref.chains


def test_homodimer_keeps_chain_ids(ref):
    a = ref.select(chains='A')
    dimer = a + a.rename_chains({'A': 'B'})
    assert dimer.map_chains(dimer) == {'A': 'A', 'B': 'B'}


def test_rename_refuses_to_merge_chains(ref):
    with pytest.raises(ValueError, match='merge'):
        ref.rename_chains({'A': 'C'})


def test_summary(ref, model):
    s = model.rename_chains(RENAME).metrics.summary(ref)
    assert s['chain_map'] == {v: k for k, v in RENAME.items()}
    assert s['lddt'] == pytest.approx(model.metrics.lddt(ref))


def test_low_identity_chain_is_skipped(ref):
    with pytest.warns(UserWarning, match='identity'):
        with pytest.raises(ValueError, match='no atoms in common'):
            ref.select(chains='A').match(ref.select(chains='D'), chains={'A': 'D'})


# ---------------------------------------------------------------- superposing

def test_superpose_recovers_rigid_motion(ref):
    moved = _rotated(ref)
    assert moved.metrics.rmsd(ref, superpose=False) > 1
    assert moved.superpose(ref).metrics.rmsd(ref, superpose=False) < 1e-6
    assert moved.metrics.rmsd(ref) < 1e-6


def test_kabsch_is_a_rotation():
    rng = np.random.default_rng(1)
    R, _ = kabsch(rng.normal(size=(20, 3)), rng.normal(size=(20, 3)))
    assert np.isclose(np.linalg.det(R), 1)


# -------------------------------------------------------------------- metrics

def test_self_scores(ref):
    assert ref.metrics.rmsd(ref) < 1e-6
    assert ref.metrics.tm_score(ref) == pytest.approx(1)
    assert ref.metrics.lddt(ref) == pytest.approx(1)
    assert ref.metrics.dockq(ref, 'A', 'C')['dockq'] == pytest.approx(1)


def test_lddt_matches_source(ref, model):
    assert model.metrics.lddt(ref, atoms='all') == pytest.approx(0.756971, abs=1e-5)
    per_res = model.metrics.lddt(ref, per_residue=True)
    assert len(per_res) == 623 and per_res['lddt'].between(0, 1).all()


def test_tm_score(ref, model):
    assert model.metrics.tm_score(ref, superpose=False) == pytest.approx(0.500526, abs=1e-5)
    # chain A on its own, against TM-align's 0.9775
    assert model.metrics.tm_score(ref.select(chains='A'), chains=['A']) == pytest.approx(0.977, abs=2e-3)


@pytest.mark.parametrize('side1, side2, dockq, lrmsd', [
    ('A', 'C', 0.914, 1.383),        # official DockQ: 0.914 / 1.383
    ('A', 'D', 0.278, 44.477),       # D is larger, so A is the ligand
    ('D', 'E', 0.252, 33.280),
])
def test_dockq_matches_official(ref, model, side1, side2, dockq, lrmsd):
    out = model.metrics.dockq(ref, side1, side2)
    assert out['dockq'] == pytest.approx(dockq, abs=1e-3)
    assert out['lrmsd'] == pytest.approx(lrmsd, abs=1e-3)


def test_plddt(model):
    assert model.metrics.plddt() == pytest.approx(89.86, abs=0.01)
    assert len(model.metrics.plddt(per_residue=True)) == 623


def test_view(ref, model):
    pytest.importorskip('py3Dmol')
    assert ref.py3Dmol(model) is not None


# --------------------------------------------------- Biopython, formats, copying

def test_from_biopython(ref):
    from Bio.PDB import PDBParser
    structure = PDBParser(QUIET=True).get_structure('x', str(DATA / '3tfk_rcsb.pdb'))
    pf = PDBFrame.from_biopython(structure)
    assert pf.name == 'x' and pf.sequences == ref.sequences
    np.testing.assert_allclose(pf.coords, ref.coords)
    assert PDBFrame.from_biopython(structure[0]).sequences == ref.sequences   # one Model
    assert PDBFrame.ensure(structure).sequences == ref.sequences
    assert PDBFrame.from_biopython(ref.to_biopython()).sequences == ref.sequences  # round trip
    with pytest.raises(TypeError, match='Structure or Model'):
        PDBFrame.from_biopython(structure[0]['A'])                           # a Chain


def test_cif_to_pdb(tmp_path):
    cif = PDBFrame.read(DATA / 'aurka_tpx2_af3.cif')
    # its two phosphothreonines are HETATM records, so not read: PDB cannot list them
    with pytest.warns(UserWarning, match='chains A \\(2\\)'):
        back = PDBFrame.read(cif.write(tmp_path / 'x.pdb'))
    assert back.sequences == cif.sequences
    np.testing.assert_allclose(back.coords, cif.coords, atol=1e-3)


def test_derived_frames_are_independent(ref):
    sub = ref.select(chains='A')
    moved = ref.transform(np.eye(3), [1, 0, 0])
    x0 = ref.df.loc[0, 'x']
    sub.df.loc[0, 'x'] = 999.0             # changing a child ...
    moved.df.loc[0, 'y'] = 999.0
    assert ref.df.loc[0, 'x'] == x0        # ... never reaches the parent
    assert ref.df.loc[0, 'y'] != 999.0
    assert moved.df.loc[0, 'x'] == pytest.approx(x0 + 1)


def test_parent_changes_do_not_reach_children(ref):
    parent = ref.copy()
    child = parent.rename_chains({'A': 'Z'})
    before = child.df.loc[0, 'x']
    parent.df.loc[0, 'x'] = -1.0
    assert child.df.loc[0, 'x'] == before


def test_copy_is_deep(ref):
    c = ref.copy()
    assert not np.shares_memory(c.df['x'].to_numpy(), ref.df['x'].to_numpy())


def test_derived_frames_keep_the_subclass(ref):
    class Sub(PDBFrame):
        pass

    sub = Sub.read(DATA / '3tfk_rcsb.pdb')
    a, b = sub.chains[:2]
    derived = [sub.select(chains=a), sub.copy(), sub.renumber(),
               sub.transform(np.eye(3), np.zeros(3)), sub.rename_chains({a: 'Z'}),
               sub.select(chains=a) + sub.select(chains=b)]
    assert all(type(d) is Sub for d in derived)
    assert type(ref.select(chains=a)) is PDBFrame

def test_observed_index_matches_reference_implementation(ref):
    from pdbframe.structure.io import observed_index

    def reference(df):      # the original per-chain factorize version
        key = df['residue_number'].astype(str) + df['insertion_code']
        return (df.assign(_key=key).groupby('chain', sort=False)['_key']
                .transform(lambda s: pd.factorize(s)[0] + 1).to_numpy())

    odd = ref.df.copy()
    odd.loc[odd.index[:30], 'insertion_code'] = 'A'
    odd.loc[odd.index[30:60], 'residue_number'] = -5
    for df in (ref.df, odd):
        assert (observed_index(df) == reference(df)).all()


# ---------------------------------------------------- chain matching at scale

def test_assembly_maps_each_copy_to_itself(ref):
    # 10 copies of the complex: 40 chains, only 4 distinct sequences
    asm = PDBFrame(pd.concat([ref.df.assign(chain=ref.df['chain'] + str(i)) for i in range(10)],
                             ignore_index=True))
    m = asm.map_chains(asm)
    assert len(m) == 40 and all(k == v for k, v in m.items())


def test_kmer_screen_keeps_true_matches(ref, model):
    from pdbframe.structure.matching import chain_map
    # screened (default 0.7) and unscreened (0.6) mapping agree on a real model
    assert chain_map(model.df, ref.df) == chain_map(model.df, ref.df, min_identity=0.6)


# ------------------------------------------------------------ RCSB (needs internet)

def _online():
    import urllib.request
    try:
        urllib.request.urlopen('https://files.rcsb.org', timeout=5)
        return True
    except Exception:
        return False


@pytest.mark.skipif(not _online(), reason='no internet access to files.rcsb.org')
def test_from_rcsb(tmp_path):
    pf = PDBFrame.from_rcsb('3tfk', save=tmp_path)
    assert pf.name == '3TFK' and 'QLSDVPMDL' in pf.sequences.values()
    assert (tmp_path / '3TFK.cif').exists()
    assert PDBFrame.from_rcsb('3TFK', fmt='pdb').sequences == pf.sequences
    with pytest.raises(ValueError, match='no cif file'):
        PDBFrame.from_rcsb('3tfk', assembly=9)


@pytest.mark.skipif(not _online(), reason='no internet access to files.rcsb.org')
def test_rcsb_info():
    table = PDBFrame.rcsb_info(['3tfk', '1AO7'])
    assert list(table['pdb_id']) == ['3TFK', '1AO7']
    assert table['method'].eq('X-ray').all() and table['resolution'].between(1, 4).all()
    assert table['release_date'].dt.year.tolist() == [2011, 1997]
    with pytest.raises(ValueError, match='no entry'):
        PDBFrame.rcsb_info('0AAA')      # IDs never start with 0
    with pytest.raises(ValueError, match='not a PDB ID'):
        PDBFrame.rcsb_info('3tfk_')
    pf = PDBFrame.from_rcsb('3tfk', metadata=True)
    assert pf.metadata['resolution'] == table['resolution'][0]
    assert pf.select(chains='A').metadata == pf.metadata         # kept by derived frames
    assert PDBFrame.from_rcsb('3tfk').metadata == {}             # off by default


def test_metadata_is_kept_by_derived_frames(ref):
    pf = PDBFrame(ref.df, name='x', metadata={'resolution': 2.1})
    for derived in [pf.select(chains='A'), pf.transform(np.eye(3), np.ones(3)), pf.copy(),
                    pf.rename_chains({'A': 'Z'})]:
        assert derived.metadata == {'resolution': 2.1}
    assert ref.metadata == {} and (pf.select(chains='A') + pf.select(chains='B')).metadata == {}
    pf.select(chains='A').metadata['resolution'] = 9     # each frame has its own dict
    assert pf.metadata == {'resolution': 2.1}


def test_rcsb_urls():
    from pdbframe.structure import rcsb
    assert rcsb.url('3tfk') == 'https://files.rcsb.org/download/3TFK.cif'
    assert rcsb.url('3tfk', assembly=1) == 'https://files.rcsb.org/download/3TFK-assembly1.cif'
    assert rcsb.url('3tfk', assembly=1, fmt='pdb') == 'https://files.rcsb.org/download/3TFK.pdb1'
    with pytest.raises(ValueError, match='not a PDB ID'):
        rcsb.url('3tf')


# ------------------------------------------------------ sequences without coordinates

RCSB_FASTA = """>3TFK_2|Chain B|p4B10 peptide|null
QLSDVPMDL
>3TFK_1|Chains A, C[auth X]|H2-Ld|Mus musculus (10090)
MGPHSMRY
"""


def test_parse_fasta():
    from pdbframe.structure.io import parse_fasta
    assert parse_fasta(RCSB_FASTA) == {'A': 'MGPHSMRY', 'X': 'MGPHSMRY', 'B': 'QLSDVPMDL'}
    assert parse_fasta('>first desc\nMKV\nLW\n>second\nAA*\n') == {'first': 'MKVLW', 'second': 'AA'}
    with pytest.raises(ValueError, match='twice'):
        parse_fasta('>a\nM\n>a\nK\n')


def test_sequence_frame(tmp_path):
    path = tmp_path / 'x.fasta'
    path.write_text(RCSB_FASTA)
    pf = PDBFrame.read(path)
    assert not pf.has_coordinates and pf.chains == ['A', 'X', 'B']
    assert pf.sequences['B'] == 'QLSDVPMDL'
    sub = pf.select(chains='A', residues=(2, 4))
    assert sub.sequences == {'A': 'GPH'} and sub.residues['residue_number'].tolist() == [2, 3, 4]
    assert pf.rename_chains({'X': 'C'}).chains == ['A', 'C', 'B']
    back = PDBFrame.read(pf.write(tmp_path / 'y.fasta'))           # round trip keeps chains
    assert back.sequences == pf.sequences
    assert PDBFrame.from_sequences({'H': 'MKVX'}).residues['residue_name'].tolist() == \
        ['MET', 'LYS', 'VAL', 'UNK']
    assert 'sequences only' in repr(pf)


def test_sequence_frame_maps_onto_structure(ref):
    seq = PDBFrame.from_sequences({'A': 'MG' + ref.sequences['A']})   # two extra leading residues
    mine, theirs = seq.match(ref, chains={'A': 'A'})
    assert len(mine) == len(ref.sequences['A'])                         # every modeled residue paired
    assert mine.residues['residue_number'].iloc[0] == 3


def test_sequence_frame_needs_coordinates_errors(ref):
    seq = PDBFrame.from_sequences(ref.sequences)
    for call in (lambda: seq.metrics.rmsd(ref), lambda: ref.metrics.lddt(seq),
                 lambda: seq.view(), lambda: seq.py3Dmol(), lambda: seq.to_string(),
                 lambda: seq.near(ref), lambda: seq.contacts()):
        with pytest.raises(ValueError, match='no coordinates'):
            call()


@pytest.mark.skipif(not _online(), reason='no internet access to files.rcsb.org')
def test_from_rcsb_fasta(tmp_path):
    seq = PDBFrame.from_rcsb('3tfk', fmt='fasta', save=tmp_path)
    assert not seq.has_coordinates and seq.chains == ['A', 'B', 'C', 'D']
    assert seq.sequences['B'] == 'QLSDVPMDL' and (tmp_path / '3TFK.fasta').exists()
    with pytest.raises(ValueError, match='assembly'):
        PDBFrame.from_rcsb('3tfk', fmt='fasta', assembly=1)


def _gapped(tmp_path, n=3):
    """aurka_tpx2_af3.cif with the atoms of `n` residues (from the 11th) of its first
    chain removed: its sequence table still lists them, as unobserved."""
    lines = (DATA / 'aurka_tpx2_af3.cif').read_text().splitlines(keepends=True)
    full = PDBFrame.read(DATA / 'aurka_tpx2_af3.cif')
    chain = full.chains[0]
    gone = list(full.residues.query('chain == @chain')['residue_number'].iloc[10:10 + n])

    def dropped(line):
        parts = line.split()
        return (line.startswith('ATOM') and parts[6] == chain     # label_asym_id
                and int(parts[8]) in gone)                       # label_seq_id
    path = tmp_path / 'gapped.cif'
    path.write_text(''.join(line for line in lines if not dropped(line)))
    return full, PDBFrame.read(path), chain, gone


def test_unobserved_residues_are_in_the_residue_table(tmp_path):
    full, pf, chain, gone = _gapped(tmp_path)
    res, before = pf.residues, full.residues
    assert {'observed', 'full_index'} <= set(res.columns)
    # the same residues in the same order; the three have no atoms (the file's
    # phosphothreonines are HETATM records, so unobserved in both without hetero=True)
    assert res[['chain', 'residue_number', 'aa']].equals(before[['chain', 'residue_number', 'aa']])
    newly = res['observed'].to_numpy() != before['observed'].to_numpy()
    assert res.loc[newly, 'residue_number'].tolist() == gone
    assert (res.loc[newly, 'atoms'] == 0).all() and res.loc[newly, 'observed_index'].isna().all()
    assert pf.unobserved['residue_number'].isin(gone).sum() == 3
    assert pf.full_sequences == full.full_sequences
    assert str(pf.seqframe(chain, full=True)) == str(full.seqframe(chain, full=True))
    assert len(pf.seqframes(full=True)[chain]) == len(pf.seqframe(chain)) + 3 + 2
    assert len(pf.sequences[chain]) == len(full.sequences[chain]) - 3
    assert len(pf.df) < len(full.df)            # the atom table never has them
    plddt = pf.metrics.plddt(per_residue=True)
    assert len(plddt) == len(res) and plddt.loc[newly, 'plddt'].isna().all()
    # formats without a deposited sequence: every residue is observed
    other = PDBFrame.read(DATA / '3tfk_rcsb.pdb')
    assert other.residues['observed'].all() and other.full_sequences == other.sequences


def test_selections_keep_the_unobserved_residues_inside_them(tmp_path):
    _, pf, chain, gone = _gapped(tmp_path)
    lo, hi = gone[0], gone[-1]

    def kept(frame):            # the removed residues a selection keeps
        unseen = frame.unobserved
        return [n for n in unseen.loc[unseen['chain'] == chain, 'residue_number'] if n in gone]

    num = pf.df['residue_number']
    on_chain = (pf.df['chain'] == chain).to_numpy()
    assert kept(pf.select(chains=chain)) == gone
    assert kept(pf.select(atoms='ca')) == gone
    assert kept(pf.select(chains=chain, residues=(lo - 5, lo))) == [lo]   # by number
    assert kept(pf.select(chains=chain, residues=(lo - 5, hi + 5))) == gone
    assert kept(pf.select(chains=chain, invert=True)) == []
    # masks: a gap is kept when both its sides are
    assert kept(pf.select(mask=on_chain & num.between(lo - 5, hi + 5).to_numpy())) == gone
    assert kept(pf.select(mask=on_chain & num.between(lo - 5, lo - 1).to_numpy())) == []
    assert kept(pf.near(pf.select(chains=chain, residues=(lo - 2, hi + 2)), 3.0)) == gone
    # selecting again from a selection, and pairing atoms with another structure
    part = pf.select(chains=chain, residues=(lo - 5, hi + 5))
    assert kept(part.select(residues=(lo - 5, lo - 1))) == []
    mine, _ = pf.match(pf)
    assert kept(mine) == gone


def test_unobserved_residues_follow_renaming_renumbering_and_joining(tmp_path):
    _, pf, chain, gone = _gapped(tmp_path)
    renamed = pf.rename_chains({chain: 'Z'})
    assert renamed.unobserved.query("chain == 'Z'")['residue_number'].tolist()[:3] == gone
    numbered = pf.select(chains=chain).renumber(start=1)
    res = numbered.residues
    assert res['residue_number'].tolist() == list(range(1, len(res) + 1))  # gap counted
    other = PDBFrame.read(DATA / '3tfk_rcsb.pdb')
    joined = pf + other.rename_chains({c: c.lower() for c in other.chains})
    assert joined.unobserved.equals(pf.unobserved)
    assert pf.copy().residues.equals(pf.residues)
    moved = pf.transform(np.eye(3), np.ones(3))
    assert moved.residues.equals(pf.residues)


def test_mmcif_writes_the_unobserved_residues(tmp_path):
    _, pf, chain, gone = _gapped(tmp_path)
    for frame in (pf, pf.rename_chains({chain: 'Z'}).renumber(start=5)):
        frame.write(tmp_path / 'out.cif')
        back = PDBFrame.read(tmp_path / 'out.cif')
        assert back.residues.equals(frame.residues)
        assert back.full_sequences == frame.full_sequences
    with pytest.warns(UserWarning, match='unobserved residues'):
        pf.write(tmp_path / 'out.pdb')
    with pytest.warns(UserWarning, match='write .cif'):
        pf.to_string('pdb')
    assert PDBFrame.read(tmp_path / 'out.pdb').residues['observed'].all()
