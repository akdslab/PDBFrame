"""The deposited sequence: residues an experiment did not resolve, kept beside the atoms.

A structure file has atoms only for the residues that were seen (in a crystal
structure, the ones with electron density). mmCIF files also list each chain's whole
deposited sequence (``_pdbx_poly_seq_scheme``), the unobserved residues included.
:class:`PDBFrame` keeps that list, so :attr:`PDBFrame.residues` has every residue of
the sequence, with ``observed`` False for those without atoms, while
:attr:`PDBFrame.df` keeps only atoms.

The list follows the frame: a selection keeps the unobserved residues that belong to
it (see :func:`subset`), renaming chains and renumbering carry over (see
:func:`remap`), and ``+`` joins the two lists. :func:`cif_tables` writes it back to
an mmCIF file.

The table (columns :data:`COLUMNS`) holds the deposited residues of the frame, in
chain and sequence order; every residue with ``observed`` True has atoms in the frame.
Extra columns (a package's labels for the residues, say) are carried along.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ['COLUMNS', 'KEYS', 'observed_residues', 'residue_table', 'subset', 'inside', 'remap',
           'renumbered', 'concat', 'cif_tables']

#: The deposited-sequence table, in column order: ``full_index`` is the residue's
#: 1-based position in the chain's deposited sequence, ``residue_number`` /
#: ``insertion_code`` its number in the file, ``observed`` whether it has atoms.
COLUMNS = ['chain', 'full_index', 'residue_name', 'aa', 'residue_number',
           'insertion_code', 'observed']

#: How a deposited residue is found in the atom table.
KEYS = ['chain', 'residue_number', 'insertion_code']


def _keys(df: pd.DataFrame) -> pd.MultiIndex:
    return pd.MultiIndex.from_frame(df[KEYS])


def observed_residues(df: pd.DataFrame) -> pd.DataFrame:
    """One row per residue of an atom table, in file order: numbering, name,
    one-letter code and atom count."""
    keys = ['chain', 'observed_index']
    return (df.groupby(keys, sort=False)
            .agg(residue_number=('residue_number', 'first'),
                 insertion_code=('insertion_code', 'first'),
                 residue_name=('residue_name', 'first'),
                 aa=('aa', 'first'),
                 atoms=('atom_name', 'size'))
            .reset_index())


def residue_table(df: pd.DataFrame, seqres: pd.DataFrame | None) -> pd.DataFrame:
    """Every residue: the atom table's, and the unobserved ones of `seqres` between them.

    Observed residues keep their order; each run of unobserved residues goes before
    the observed residue that follows it in the deposited sequence (after the last one
    at a chain's end). ``observed_index`` is empty and ``atoms`` 0 for unobserved
    residues; ``full_index`` is empty for residues not in a deposited sequence
    (ligands, or every residue of a file without one).
    """
    res = observed_residues(df)
    res['observed'] = True
    if seqres is None:
        res['full_index'] = pd.array([pd.NA] * len(res), dtype='Int64')
        return res
    seen = seqres.loc[seqres['observed'], KEYS + ['full_index']]
    res = res.merge(seen, on=KEYS, how='left')          # a left merge keeps the order
    res['full_index'] = res['full_index'].astype('Int64')
    unseen = seqres[~seqres['observed'].to_numpy()
                    & seqres['chain'].isin(res['chain']).to_numpy()].reset_index(drop=True)
    if unseen.empty:
        return res
    # sort key: observed residues their position; unobserved ones just before the next
    # observed residue (by full_index) of their chain
    res['_key'] = np.arange(len(res), dtype=float)
    keys = np.empty(len(unseen))
    for chain, rows in unseen.groupby('chain', sort=False):
        mine = res[res['chain'] == chain]
        index = mine['full_index'].to_numpy(dtype=float, na_value=np.nan)
        pos = mine['_key'].to_numpy()
        targets = rows['full_index'].to_numpy(float)
        later = index[None, :] > targets[:, None]                 # unobserved x observed
        after = np.where(later.any(axis=1),
                         np.where(later, pos[None, :], np.inf).min(axis=1), pos[-1] + 1)
        # within one gap, in sequence order, all between the two observed residues
        order = pd.Series(targets).groupby(after).rank(method='first').to_numpy()
        size = pd.Series(after).map(pd.Series(after).value_counts()).to_numpy()
        keys[rows.index.to_numpy()] = after - 1 + order / (size + 1)
    add = pd.DataFrame({
        'chain': unseen['chain'].to_numpy(),
        'observed_index': pd.NA,
        'residue_number': unseen['residue_number'].to_numpy(),
        'insertion_code': unseen['insertion_code'].to_numpy(),
        'residue_name': unseen['residue_name'].to_numpy(),
        'aa': unseen['aa'].to_numpy(),
        'atoms': 0,
        'observed': False,
        'full_index': unseen['full_index'].to_numpy(),
        '_key': keys,
    })
    rank = {c: i for i, c in enumerate(dict.fromkeys(res['chain']))}
    out = pd.concat([res, add], ignore_index=True)
    out['_rank'] = out['chain'].map(rank)
    out = out.sort_values(['_rank', '_key'], kind='stable').drop(columns=['_rank', '_key'])
    out['observed_index'] = out['observed_index'].astype('Int64')
    out['full_index'] = out['full_index'].astype('Int64')
    out['atoms'] = out['atoms'].astype(int)
    out['observed'] = out['observed'].astype(bool)
    return out.reset_index(drop=True)


def subset(seqres: pd.DataFrame | None, df: pd.DataFrame,
           unobserved: np.ndarray | None = None) -> pd.DataFrame | None:
    """The deposited residues of a selection `df` of the frame `seqres` belongs to.

    Observed residues stay if `df` has them. Unobserved ones have no atoms to test, so:

    - `unobserved` given (one bool per row of `seqres`): those it marks, on chains
      `df` keeps. :meth:`PDBFrame.select` passes this for its residue-level filters
      (``chains``, ``residues``), which unobserved residues have.
    - otherwise (atom masks, distances, row positions): a run of unobserved residues
      stays if the selection keeps every observed residue next to it in the sequence
      (one side at a chain's end), i.e. the run lies inside what was kept.
    """
    if seqres is None:
        return None
    observed = seqres['observed'].to_numpy()
    kept = _keys(seqres).isin(_keys(df.drop_duplicates(KEYS)))
    on_chain = seqres['chain'].isin(set(df['chain'])).to_numpy()
    chosen = inside(seqres, df) if unobserved is None else np.asarray(unobserved, bool)
    keep = (observed & kept) | (~observed & chosen & on_chain)
    return seqres if keep.all() else seqres[keep].reset_index(drop=True)


def inside(seqres: pd.DataFrame, df: pd.DataFrame) -> np.ndarray:
    """True for the unobserved residues of `seqres` whose observed neighbours in the
    sequence (one at a chain's end) are all in `df`: the runs inside a selection."""
    observed = seqres['observed'].to_numpy()
    kept = _keys(seqres).isin(_keys(df.drop_duplicates(KEYS)))
    state = pd.Series(np.where(observed, kept, np.nan))
    chain = seqres['chain'].to_numpy()
    before = state.groupby(chain).ffill().fillna(1).to_numpy()
    after = state.groupby(chain).bfill().fillna(1).to_numpy()
    return ~observed & (before == 1) & (after == 1)


def remap(seqres: pd.DataFrame | None, old: pd.DataFrame, new: pd.DataFrame
          ) -> pd.DataFrame | None:
    """`seqres` for `new`, the atom table `old` changed row by row (chains renamed,
    residues renumbered, coordinates moved): observed residues take their new chain
    and number, unobserved ones the new ID of their chain."""
    if seqres is None:
        return None
    if len(old) != len(new):
        raise ValueError(f'expected the same atoms ({len(old)}), got {len(new)}')
    if old[KEYS].reset_index(drop=True).equals(new[KEYS].reset_index(drop=True)):
        return seqres
    pairs = pd.concat([old[KEYS].reset_index(drop=True),
                       new[KEYS].add_prefix('new_').reset_index(drop=True)], axis=1)
    pairs = pairs.drop_duplicates()
    if pairs.duplicated(KEYS).any():
        raise ValueError('a residue was split in two; cannot carry the deposited sequence')
    chains = pairs[['chain', 'new_chain']].drop_duplicates()
    if chains['chain'].duplicated().any():
        raise ValueError('a chain was split in two; cannot carry the deposited sequence')
    out = seqres.merge(pairs, on=KEYS, how='left')
    seen = out['observed'].to_numpy()
    if out.loc[seen, 'new_chain'].isna().any():
        raise ValueError('an observed residue is missing from the atom table')
    renamed = out['chain'].map(dict(zip(chains['chain'], chains['new_chain'])))
    out['chain'] = renamed.fillna(out['chain'])
    for col in ['residue_number', 'insertion_code']:
        out.loc[seen, col] = out.loc[seen, 'new_' + col]
    out['residue_number'] = out['residue_number'].astype(int)
    return out.drop(columns=['new_' + c for c in KEYS])


def renumbered(seqres: pd.DataFrame, numbers: pd.DataFrame) -> pd.DataFrame:
    """`seqres` with new numbers: `numbers` has :data:`KEYS` (old) and ``new_number``,
    ``new_insertion_code`` for each residue."""
    out = seqres.merge(numbers, on=KEYS, how='left')
    if out['new_number'].isna().any():
        raise ValueError('a deposited residue has no new number')
    out['residue_number'] = out.pop('new_number').astype(int)
    out['insertion_code'] = out.pop('new_insertion_code').astype(str)
    return out


def concat(*tables: pd.DataFrame | None) -> pd.DataFrame | None:
    """The deposited residues of frames joined with ``+`` (their chains never overlap)."""
    given = [t for t in tables if t is not None]
    return pd.concat(given, ignore_index=True) if given else None


def _listed(value) -> list[str]:
    """An mmCIF dictionary value as a list (a loop of one row is a plain string)."""
    return [value] if isinstance(value, str) else list(value)


def cif_tables(mmcif: dict, seqres: pd.DataFrame) -> dict:
    """`mmcif` (an mmCIF dictionary of a written structure) with the deposited sequence
    added, so the unobserved residues are written too.

    Each chain becomes its own entity: ``_entity_poly_seq`` lists its sequence and
    ``_pdbx_poly_seq_scheme`` places every residue, the unobserved ones with ``?``
    for their atoms' number and name (as RCSB writes them). The atoms of observed
    residues get their ``label_seq_id`` (the position in the sequence) and
    ``label_entity_id``.
    """
    out = dict(mmcif)
    site = '_atom_site.'
    chain, number, codes, label_seq, entity, asym = (
        _listed(out[site + f]) for f in ['auth_asym_id', 'auth_seq_id', 'pdbx_PDB_ins_code',
                                         'label_seq_id', 'label_entity_id', 'label_asym_id'])
    icode = ['' if x in ('?', '.') else x for x in codes]
    chains = list(dict.fromkeys(seqres['chain']))
    entities = {c: str(i) for i, c in enumerate(chains, 1)}
    where = {(c, int(n), i): int(s) for c, n, i, s in
             zip(seqres['chain'], seqres['residue_number'], seqres['insertion_code'],
                 seqres['full_index'])}
    label_asym = {}
    for k in range(len(chain)):
        s = where.get((chain[k], int(number[k]), icode[k]))
        if s is None:
            continue
        label_seq[k], entity[k] = str(s), entities[chain[k]]
        label_asym.setdefault(chain[k], asym[k])
    out[site + 'label_seq_id'] = label_seq
    out[site + 'label_entity_id'] = entity
    seq = pd.concat([g.sort_values('full_index', kind='stable')
                     for _, g in seqres.groupby('chain', sort=False)])
    one_letter = seq['aa'].tolist()
    seen = seq['observed'].to_numpy()
    first = seq.drop_duplicates('chain')
    out['_entity_poly.entity_id'] = [entities[c] for c in first['chain']]
    out['_entity_poly.type'] = ['polypeptide(L)'] * len(first)
    out['_entity_poly.pdbx_seq_one_letter_code_can'] = [
        ''.join(a for a, c in zip(one_letter, seq['chain']) if c == x) for x in first['chain']]
    out['_entity_poly.pdbx_strand_id'] = list(first['chain'])
    out['_entity_poly_seq.entity_id'] = [entities[c] for c in seq['chain']]
    out['_entity_poly_seq.num'] = [str(s) for s in seq['full_index']]
    out['_entity_poly_seq.mon_id'] = list(seq['residue_name'])
    out['_entity_poly_seq.hetero'] = ['n'] * len(seq)
    scheme = '_pdbx_poly_seq_scheme.'
    numbers = [str(n) for n in seq['residue_number']]
    codes = [i or '.' for i in seq['insertion_code']]
    names = list(seq['residue_name'])
    out[scheme + 'asym_id'] = [label_asym.get(c, c) for c in seq['chain']]
    out[scheme + 'entity_id'] = [entities[c] for c in seq['chain']]
    out[scheme + 'seq_id'] = [str(s) for s in seq['full_index']]
    out[scheme + 'mon_id'] = names
    out[scheme + 'ndb_seq_num'] = [str(s) for s in seq['full_index']]
    out[scheme + 'pdb_seq_num'] = numbers
    out[scheme + 'auth_seq_num'] = [n if s else '?' for n, s in zip(numbers, seen)]
    out[scheme + 'pdb_mon_id'] = [m if s else '?' for m, s in zip(names, seen)]
    out[scheme + 'auth_mon_id'] = [m if s else '?' for m, s in zip(names, seen)]
    out[scheme + 'pdb_strand_id'] = list(seq['chain'])
    out[scheme + 'pdb_ins_code'] = codes
    out[scheme + 'hetero'] = ['n'] * len(seq)
    return out
