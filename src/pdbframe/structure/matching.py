"""Pairing up the atoms of two structures of the same protein.

Two structures of one protein rarely line up row for row: numbering differs, one
has missing loops, one has hydrogens. :func:`match` finds the atoms present in both
and returns them as two atom tables in the same order, so row i of one is row i of
the other. Every structure comparison (superposition, RMSD, lDDT, ...) builds on it.
"""

from __future__ import annotations

import warnings
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment

from ..sequence import SeqFrame

__all__ = ['match', 'chain_map', 'BACKBONE']

#: Protein backbone atoms - the only atoms paired between two different residue types.
BACKBONE = ('N', 'CA', 'C', 'O')


def _sequences(df: pd.DataFrame) -> dict[str, str]:
    res = df.drop_duplicates(['chain', 'observed_index'])
    return {c: ''.join(g['aa']) for c, g in res.groupby('chain', sort=False)}


@lru_cache(maxsize=4096)
def _align(seq1: str, seq2: str) -> tuple[pd.DataFrame, int, float]:
    """Aligned residue pairs, how many are identical, and that as a share of the pairs.

    Cached by sequence: assemblies repeat the same chains many times, and the
    residue pairing after :func:`chain_map` reuses its alignments.
    """
    aln = SeqFrame(seq1).align(seq2).df
    both = aln.dropna(subset=['seq1_pos', 'seq2_pos'])
    same = int((both['match'] == 'identical').sum())
    return both, same, same / len(both) if len(both) else 0.0


def _kmers(seq: str, k: int = 3) -> frozenset:
    return frozenset(seq[i:i + k] for i in range(len(seq) - k + 1))


def chain_map(df1: pd.DataFrame, df2: pd.DataFrame, min_identity: float = 0.7,
              chains=None) -> dict[str, str]:
    """Which chain of `df2` each chain of `df1` corresponds to, by sequence.

    Every chain pair is aligned and scored by its number of identical residues.
    Chains are then assigned one-to-one to maximise the total, so each chain is
    used once even for homo-oligomers; between equally good choices, equal chain
    IDs win. Pairs below `min_identity` (share of aligned residues that are
    identical) are left out.

    Large assemblies stay fast: each distinct pair of sequences is aligned once,
    and pairs sharing no 3-residue stretch are skipped without aligning. That
    skip is safe whenever `min_identity` is above 2/3 - without three identical
    aligned residues in a row, identity can be at most 2/3 - so it is turned off
    for lower thresholds.

    Args:
        chains: only map these chains of `df1` (default: all).

    Returns:
        ``{chain_in_1: chain_in_2}``; chains without a counterpart are absent.
    """
    s1, s2 = _sequences(df1), _sequences(df2)
    if chains is not None:
        chains = [chains] if isinstance(chains, str) else list(chains)
        missing = set(chains) - set(s1)
        if missing:
            raise KeyError(f'chains not found: {sorted(missing)}')
        s1 = {c: s1[c] for c in chains}
    c1, c2 = list(s1), list(s2)
    screen = min_identity > 2 / 3
    k1 = {c: _kmers(s1[c]) for c in c1} if screen else {}
    k2 = {c: _kmers(s2[c]) for c in c2} if screen else {}
    score = np.zeros((len(c1), len(c2)))
    identity = np.zeros_like(score)
    for i, a in enumerate(c1):
        for j, b in enumerate(c2):
            # sequences shorter than 3 residues have no 3-mers, so always align those
            if screen and k1[a] and k2[b] and k1[a].isdisjoint(k2[b]):
                continue
            _, same, identity[i, j] = _align(s1[a], s2[b])
            # half a residue: breaks ties toward equal IDs, never outweighs a residue
            score[i, j] = same + 0.5 * (a == b)
    score[identity < min_identity] = 0
    rows, cols = linear_sum_assignment(score, maximize=True)
    return {c1[i]: c2[j] for i, j in zip(rows, cols) if identity[i, j] >= min_identity}


def _chain_pairs(df1, df2, chains, min_identity) -> list[tuple[str, str]]:
    if isinstance(chains, dict):
        return list(chains.items())
    if chains == 'same':
        common = set(df2['chain'])
        return [(c, c) for c in df1['chain'].unique() if c in common]
    if chains is None or chains == 'auto':
        return list(chain_map(df1, df2, min_identity).items())
    return list(chain_map(df1, df2, min_identity, chains=chains).items())


def _residue_pairs(r1: pd.DataFrame, r2: pd.DataFrame, by: str, min_identity: float,
                   label: str) -> pd.DataFrame:
    """Pair the residues of one chain: observed_index in 1 -> observed_index in 2."""
    if by == 'number':
        return r1.merge(r2, on=['residue_number', 'insertion_code'],
                        suffixes=('_1', '_2'))[['observed_index_1', 'observed_index_2']]
    if by != 'sequence':
        raise ValueError(f"by must be 'sequence' or 'number', not {by!r}")

    both, _, identity = _align(''.join(r1['aa']), ''.join(r2['aa']))
    if identity < min_identity:
        warnings.warn(f'chain {label}: sequence identity {identity:.0%} is below '
                      f'min_identity={min_identity:.0%}; chain skipped', stacklevel=4)
        return pd.DataFrame(columns=['observed_index_1', 'observed_index_2'])
    # alignment positions are 1-based positions in each chain's residue list
    return pd.DataFrame({
        'observed_index_1': r1['observed_index'].to_numpy()[both['seq1_pos'].to_numpy(int) - 1],
        'observed_index_2': r2['observed_index'].to_numpy()[both['seq2_pos'].to_numpy(int) - 1],
    })


def match(df1: pd.DataFrame, df2: pd.DataFrame, by: str = 'sequence', chains='auto',
          min_identity: float = 0.7) -> tuple[np.ndarray, np.ndarray]:
    """Row positions of the atoms the two tables have in common, paired up.

    Residues are paired chain by chain, then atoms by name within each residue
    pair. When the two residues differ (a mutation), only backbone atoms are
    paired, since side-chain atom names mean different things.

    Args:
        df1, df2: atom tables.
        by: ``'sequence'`` pairs residues by aligning each chain's sequence, so
            numbering may differ and residues may be missing; ``'number'`` pairs
            residues with equal residue number and insertion code.
        chains: which chains to pair.

            - ``'auto'`` (default): find each chain's counterpart by sequence with
              :func:`chain_map`, so chain IDs need not agree.
            - a list of chains of `df1`: only those, mapped automatically.
            - ``'same'``: pair equal chain IDs only.
            - a dict ``{chain_in_1: chain_in_2}``: exactly these pairs.
        min_identity: with ``by='sequence'``, skip (with a warning) chains whose
            aligned residues are less identical than this.

    Returns:
        ``(rows1, rows2)``: integer positions into `df1` and `df2`, same length,
        with ``df1.iloc[rows1]`` and ``df2.iloc[rows2]`` the paired atoms.
    """
    res_cols = ['observed_index', 'residue_number', 'insertion_code', 'aa']
    all1 = df1.reset_index(drop=True).assign(_row=np.arange(len(df1)))
    all2 = df2.reset_index(drop=True).assign(_row=np.arange(len(df2)))
    # split by chain once: filtering the whole table per chain is slow for big assemblies
    by1 = dict(tuple(all1.groupby('chain', sort=False)))
    by2 = dict(tuple(all2.groupby('chain', sort=False)))
    rows1, rows2 = [], []
    for c1, c2 in _chain_pairs(df1, df2, chains, min_identity):
        if c1 not in by1 or c2 not in by2:
            raise KeyError(f'chain {c1 if c1 not in by1 else c2!r} not found')
        a1, a2 = by1[c1], by2[c2]
        r1 = a1.drop_duplicates('observed_index')[res_cols]
        r2 = a2.drop_duplicates('observed_index')[res_cols]
        pairs = _residue_pairs(r1, r2, by, min_identity, c1 if c1 == c2 else f'{c1}/{c2}')

        pairs = pairs.set_axis(['res1', 'res2'], axis=1)
        cols = ['observed_index', 'residue_name', 'atom_name', '_row']
        left = a1[cols].set_axis(['res1', 'name1', 'atom_name', 'row1'], axis=1)
        right = a2[cols].set_axis(['res2', 'name2', 'atom_name', 'row2'], axis=1)
        m = left.merge(pairs, on='res1').merge(right, on=['res2', 'atom_name'])
        m = m[(m['name1'] == m['name2']) | m['atom_name'].isin(BACKBONE)].sort_values('row1')
        rows1.append(m['row1'].to_numpy())
        rows2.append(m['row2'].to_numpy())
    if not rows1:
        return np.array([], int), np.array([], int)
    return np.concatenate(rows1), np.concatenate(rows2)
