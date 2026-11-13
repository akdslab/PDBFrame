"""Picking atoms and residues out of an atom table, and finding contacts.

The functions work on plain atom tables (see :mod:`.io`); :class:`PDBFrame`'s
``select``, ``near`` and ``contacts`` methods call them.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .geometry import atom_pairs_within
from .matching import BACKBONE

__all__ = ['ATOM_SETS', 'atom_mask', 'residue_labels', 'select_mask', 'near_mask', 'contacts']

#: Named atom sets accepted by ``atoms=``; a list of atom names also works.
ATOM_SETS = ('all', 'heavy', 'backbone', 'ca')

_XYZ = ['x', 'y', 'z']


def atom_mask(df: pd.DataFrame, atoms) -> np.ndarray:
    """True for the atoms in a named set (:data:`ATOM_SETS`) or a list of atom names."""
    if isinstance(atoms, str):
        if atoms == 'all':
            return np.ones(len(df), bool)
        if atoms == 'heavy':
            return ~df['element'].isin(['H', 'D']).to_numpy()
        if atoms == 'backbone':
            return df['atom_name'].isin(BACKBONE).to_numpy()
        if atoms == 'ca':
            return (df['atom_name'] == 'CA').to_numpy()
        raise ValueError(f'atoms must be one of {ATOM_SETS} or a list of atom names, '
                         f'not {atoms!r}')
    return df['atom_name'].isin(list(atoms)).to_numpy()


def residue_labels(df: pd.DataFrame) -> np.ndarray:
    """One label per atom naming its residue, unique across chains."""
    return (df['chain'] + ':' + df['observed_index'].astype(str)).to_numpy()


def select_mask(df: pd.DataFrame, chains=None, residues=None, atoms=None, mask=None) -> np.ndarray:
    """True for the atoms matching every filter given; see :meth:`PDBFrame.select`."""
    keep = np.ones(len(df), bool)
    if chains is not None:
        keep &= df['chain'].isin([chains] if isinstance(chains, str) else chains).to_numpy()
    if residues is not None:
        num = df['residue_number']
        if isinstance(residues, tuple):
            start, end = residues
            keep &= num.between(start, end).to_numpy()
        else:
            keep &= num.isin([residues] if np.isscalar(residues) else list(residues)).to_numpy()
    if atoms is not None:
        keep &= atom_mask(df, atoms)
    if mask is not None:
        keep &= np.asarray(mask, bool)
    return keep


def near_mask(df: pd.DataFrame, other: pd.DataFrame, cutoff: float) -> np.ndarray:
    """True for every atom of `df`'s residues that have any atom within `cutoff` Å of `other`."""
    i, _, _ = atom_pairs_within(df[_XYZ].to_numpy(float), other[_XYZ].to_numpy(float), cutoff)
    labels = residue_labels(df)
    return np.isin(labels, labels[i])


def contacts(df: pd.DataFrame, other: pd.DataFrame | None = None,
             cutoff: float = 5.0) -> pd.DataFrame:
    """Residue pairs in contact, with their closest atom-atom distance.

    With `other`, pairs between `df` and `other`; without it, pairs between
    different chains of `df` (its interfaces), each listed once.
    """
    a = df
    b = a if other is None else other
    i, j, d = atom_pairs_within(a[_XYZ].to_numpy(float), b[_XYZ].to_numpy(float), cutoff)
    if other is None:
        # each interface pair once, from the earlier chain to the later one
        order = {c: k for k, c in enumerate(a['chain'].unique())}
        rank = a['chain'].map(order).to_numpy()
        keep = rank[i] < rank[j]
        i, j, d = i[keep], j[keep], d[keep]
    cols = ['chain', 'residue_number', 'insertion_code', 'residue_name']
    left = a.iloc[i][cols].add_suffix('_1').reset_index(drop=True)
    right = b.iloc[j][cols].add_suffix('_2').reset_index(drop=True)
    pairs = pd.concat([left, right], axis=1).assign(distance=d)
    keys = [c for c in pairs.columns if c != 'distance']
    return (pairs.groupby(keys, sort=False, as_index=False)['distance'].min()
            .sort_values(['chain_1', 'residue_number_1', 'chain_2', 'residue_number_2'])
            .reset_index(drop=True))
