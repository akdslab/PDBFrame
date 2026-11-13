"""Structure scores: RMSD, TM-score, lDDT, DockQ and pLDDT.

All metric code lives here, in two layers:

- :class:`Metrics`, reached as ``pf.metrics``: compares a structure with a
  reference. It pairs up their atoms (chains and residues matched by sequence; see
  :mod:`.matching`), then calls the functions below.
- ``*_coords`` functions: the maths, on ``(N, 3)`` arrays whose rows are already
  paired. Use them directly to score coordinates from elsewhere.

    >>> model.metrics.rmsd(ref)
    >>> model.metrics.summary(ref)        # every score in one Series
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from .geometry import atom_pairs_within, kabsch, rmsd, transform
from .matching import BACKBONE
from .pdbframe import PDBFrame
from .selection import residue_labels

__all__ = ['Metrics', 'rmsd', 'tm_score_coords', 'lddt_coords', 'dockq_coords']

_XYZ = ['x', 'y', 'z']


class Metrics:
    """Scores of one structure (the model) against a reference: ``pf.metrics``.

    Every comparison takes the reference first, then options. The ``match_kw``
    (``chains``, ``by``, ``min_identity``) go to :meth:`PDBFrame.match`: by
    default each model chain is paired with its counterpart in the reference by
    sequence, so chain IDs and residue numbering do not need to agree.
    """

    def __init__(self, frame: PDBFrame):
        self.frame = frame

    def __repr__(self) -> str:
        return (f'Metrics({self.frame.name!r}): rmsd, tm_score, lddt, dockq, plddt, '
                'summary')

    def _xyz(self, reference, atoms, match_kw):
        a, b = self.frame._paired(PDBFrame.ensure(reference), atoms, match_kw)
        return a[_XYZ].to_numpy(float), b[_XYZ].to_numpy(float), a, b

    def rmsd(self, reference: PDBFrame, atoms='ca', superpose: bool = True,
             **match_kw) -> float:
        """RMSD (Å) over the common `atoms`, after superposing them unless ``superpose=False``."""
        a, b, _, _ = self._xyz(reference, atoms, match_kw)
        return rmsd(a, b, superpose=superpose)

    def tm_score(self, reference: PDBFrame, superpose: bool = True, **match_kw) -> float:
        """TM-score (0-1) against `reference`, normalised by the reference's length.

        See :func:`tm_score_coords`; ``superpose=False`` scores the current pose.
        """
        reference = PDBFrame.ensure(reference)
        a, b, _, _ = self._xyz(reference, 'ca', match_kw)
        length = int((reference.df['atom_name'] == 'CA').sum())
        return tm_score_coords(a, b, length=length, superpose=superpose)

    def lddt(self, reference: PDBFrame, atoms='heavy', cutoff: float = 15.0,
             per_residue: bool = False, **match_kw):
        """lDDT (0-1) against `reference`; no superposition needed.

        With ``per_residue=True``, returns the reference's residue table with an
        ``lddt`` column instead of the global score.
        """
        a, b, _, ref_atoms = self._xyz(reference, atoms, match_kw)
        score = lddt_coords(a, b, residue_labels(ref_atoms), cutoff=cutoff,
                            per_residue=per_residue)
        if not per_residue:
            return score
        return PDBFrame(ref_atoms).residues.assign(lddt=score)

    def dockq(self, reference: PDBFrame, side1, side2=None, **match_kw) -> dict:
        """DockQ of the interface between `side1` and `side2`, against `reference`.

        Args:
            side1: chain ID(s) of the reference on one side of the interface.
            side2: chain ID(s) of the reference on the other side; default all
                its other chains.

        As in the DockQ program, the side with fewer residues is the ligand for
        LRMSD. Returns ``dockq``, ``fnat``, ``irmsd``, ``lrmsd``, ``fnonnat`` and
        the contact counts; see :func:`dockq_coords`.
        """
        a, b, _, ref_atoms = self._xyz(reference, 'heavy', match_kw)
        side1 = [side1] if isinstance(side1, str) else list(side1)
        side2 = ([side2] if isinstance(side2, str) else list(side2)) if side2 \
            else [c for c in ref_atoms['chain'].unique() if c not in side1]
        keep = ref_atoms['chain'].isin(side1 + side2).to_numpy()
        a, b, ref_atoms = a[keep], b[keep], ref_atoms[keep]
        n_res = ref_atoms.drop_duplicates(['chain', 'observed_index'])['chain'].value_counts()
        ligand = side2 if n_res.reindex(side2).sum() <= n_res.reindex(side1).sum() else side1
        return dockq_coords(a, b, residue=residue_labels(ref_atoms),
                            ligand=ref_atoms['chain'].isin(ligand).to_numpy(),
                            backbone=ref_atoms['atom_name'].isin(BACKBONE).to_numpy())

    def plddt(self, per_residue: bool = False):
        """Mean B-factor per residue - the pLDDT of AlphaFold-style predicted models.

        Needs no reference. Returns the mean over residues, or the residue table
        with a ``plddt`` column when ``per_residue=True`` (NaN for unobserved residues).
        """
        df = self.frame.df
        per_res = df.groupby(['chain', 'observed_index'], sort=False)['bfactor'].mean()
        if not per_residue:
            return float(per_res.mean())
        res = self.frame.residues
        res['plddt'] = np.nan
        res.loc[res['observed'], 'plddt'] = per_res.to_numpy()    # observed keep their order
        return res

    def summary(self, reference: PDBFrame, **match_kw) -> pd.Series:
        """RMSD (Cα), TM-score, lDDT and pLDDT in one Series, plus the chain mapping used."""
        reference = PDBFrame.ensure(reference)
        mine, theirs = self.frame.match(reference, **match_kw)
        mapping = dict(zip(mine.chains, theirs.chains))
        return pd.Series({
            'rmsd': self.rmsd(reference, **match_kw),
            'tm_score': self.tm_score(reference, **match_kw),
            'lddt': self.lddt(reference, **match_kw),
            'plddt': self.plddt(),
            'chain_map': mapping,
        }, name=self.frame.name)


# ------------------------------------------------------------- coordinate maths

def _d0(length: int) -> float:
    """TM-score distance scale for a protein of `length` residues."""
    return max(1.24 * np.cbrt(length - 15) - 1.8, 0.5) if length > 15 else 0.5


def tm_score_coords(model: np.ndarray, ref: np.ndarray, length: int | None = None,
                    superpose: bool = True, max_iter: int = 20) -> float:
    """TM-score of paired Cα coordinates (one row per residue).

    Args:
        length: residues in the reference, which the score is normalised by.
            Defaults to the number of pairs; pass the full reference length so
            unpaired reference residues count against the model.
        superpose: search for the superposition that maximises the score - start
            from the best fit of all pairs, then refit on the pairs closer than d0
            until that set stops changing. Without it, score the current pose.

    The search is the simple heuristic, so the score is a lower bound on what
    TM-score / TM-align report (usually equal for models of the same protein).
    """
    length = length or len(ref)
    d0 = _d0(length)

    def score(xyz):
        d = np.linalg.norm(xyz - ref, axis=1)
        return float((1 / (1 + (d / d0) ** 2)).sum() / length), d

    if not superpose:
        return score(model)[0]
    best, subset = 0.0, np.ones(len(model), bool)
    for _ in range(max_iter):
        value, d = score(transform(model, *kabsch(model[subset], ref[subset])))
        best = max(best, value)
        new = d < d0
        if new.sum() < 3 or (new == subset).all():
            break
        subset = new
    return best


def lddt_coords(model: np.ndarray, ref: np.ndarray, residue: np.ndarray,
                cutoff: float = 15.0, thresholds=(0.5, 1.0, 2.0, 4.0),
                per_residue: bool = False):
    """Local Distance Difference Test - superposition-free local accuracy, 0 to 1.

    Takes every pair of atoms in different residues that are within `cutoff` Å in
    the reference, and scores how many of those distances the model keeps within
    each threshold.

    Args:
        residue: residue label per atom (any hashable), so same-residue pairs are
            skipped.
        per_residue: return one score per residue (mean over the pairs touching
            it, in order of first appearance) instead of the global score.
    """
    i, j = cKDTree(ref).query_pairs(cutoff, output_type='ndarray').T
    codes, uniques = _codes(residue)
    keep = codes[i] != codes[j]
    i, j = i[keep], j[keep]
    d_ref = np.linalg.norm(ref[i] - ref[j], axis=1)
    d_mod = np.linalg.norm(model[i] - model[j], axis=1)
    pair_score = (np.abs(d_mod - d_ref)[:, None] < np.asarray(thresholds)).mean(axis=1)
    if not per_residue:
        return float(pair_score.mean()) if len(pair_score) else 0.0
    n = len(uniques)
    total = np.bincount(codes[i], pair_score, n) + np.bincount(codes[j], pair_score, n)
    count = np.bincount(codes[i], minlength=n) + np.bincount(codes[j], minlength=n)
    with np.errstate(invalid='ignore'):
        return total / count


def _codes(labels: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Integer code per label, numbered in order of first appearance."""
    uniques, first, codes = np.unique(labels, return_index=True, return_inverse=True)
    order = np.argsort(first)
    rank = np.empty_like(order)
    rank[order] = np.arange(len(order))
    return rank[codes], uniques[order]


def dockq_coords(model: np.ndarray, ref: np.ndarray, residue: np.ndarray,
                 ligand: np.ndarray, backbone: np.ndarray, contact_cutoff: float = 5.0,
                 interface_cutoff: float = 10.0) -> dict:
    """DockQ for a two-body interface (Basu & Wallner 2016).

    Args:
        model, ref: paired heavy-atom coordinates.
        residue: residue label per atom.
        ligand: True for atoms of the ligand side, False for the receptor side.
        backbone: True for backbone atoms.
        contact_cutoff: residues touch when any heavy atoms are this close (Å).
        interface_cutoff: interface residues are this close to the partner (Å).

    Returns:
        ``dockq``, its three parts - ``fnat`` (share of reference contacts the model
        keeps), ``irmsd`` (backbone RMSD of interface residues after superposing
        them) and ``lrmsd`` (ligand backbone RMSD after superposing the receptor) -
        plus ``fnonnat`` and the contact counts.
    """
    ligand = np.asarray(ligand, bool)
    backbone = np.asarray(backbone, bool)
    rec_idx, lig_idx = np.flatnonzero(~ligand), np.flatnonzero(ligand)

    def contacts(xyz, cutoff):
        i, j, _ = atom_pairs_within(xyz[rec_idx], xyz[lig_idx], cutoff)
        return set(zip(residue[rec_idx[i]], residue[lig_idx[j]]))

    native, predicted = contacts(ref, contact_cutoff), contacts(model, contact_cutoff)
    fnat = len(native & predicted) / len(native) if native else 0.0
    fnonnat = len(predicted - native) / len(predicted) if predicted else 0.0

    interface = {r for pair in contacts(ref, interface_cutoff) for r in pair}
    sel = backbone & np.isin(residue, list(interface))
    irmsd = rmsd(model[sel], ref[sel], superpose=True) if sel.any() else np.nan

    rec_bb, lig_bb = backbone & ~ligand, backbone & ligand
    moved = transform(model, *kabsch(model[rec_bb], ref[rec_bb]))
    lrmsd = rmsd(moved[lig_bb], ref[lig_bb])

    score = (fnat + 1 / (1 + (irmsd / 1.5) ** 2) + 1 / (1 + (lrmsd / 8.5) ** 2)) / 3
    return {'dockq': float(score), 'fnat': fnat, 'irmsd': float(irmsd),
            'lrmsd': float(lrmsd), 'fnonnat': fnonnat,
            'native_contacts': len(native), 'model_contacts': len(predicted)}
