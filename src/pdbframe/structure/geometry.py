"""Coordinate maths on plain ``(N, 3)`` numpy arrays: superposition and distances."""

from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

__all__ = ['kabsch', 'transform', 'rmsd', 'atom_pairs_within']


def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotation and translation that best superpose `mobile` onto `target`.

    Both are ``(N, 3)`` with row i of one paired to row i of the other. Returns
    ``(R, t)`` such that ``mobile @ R.T + t`` is the least-squares fit to `target`.
    """
    if mobile.shape != target.shape or len(mobile) == 0:
        raise ValueError(f'need two equal, non-empty (N, 3) arrays, got '
                         f'{mobile.shape} and {target.shape}')
    mc, tc = mobile.mean(axis=0), target.mean(axis=0)
    H = (mobile - mc).T @ (target - tc)
    U, _, Vt = np.linalg.svd(H)
    # flip the last axis if needed so R is a rotation, not a reflection
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, tc - mc @ R.T


def transform(xyz: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Apply a rotation and translation from :func:`kabsch`."""
    return xyz @ R.T + t


def rmsd(a: np.ndarray, b: np.ndarray, superpose: bool = False) -> float:
    """Root-mean-square deviation of paired points, optionally after superposition."""
    if a.shape != b.shape or len(a) == 0:
        raise ValueError(f'need two equal, non-empty (N, 3) arrays, got {a.shape} and {b.shape}')
    if superpose:
        a = transform(a, *kabsch(a, b))
    return float(np.sqrt(((a - b) ** 2).sum(axis=1).mean()))


def atom_pairs_within(a: np.ndarray, b: np.ndarray, cutoff: float):
    """Every (i, j, distance) with ``|a[i] - b[j]| <= cutoff``, via KD-trees."""
    pairs = cKDTree(a).sparse_distance_matrix(cKDTree(b), cutoff, output_type='ndarray')
    return pairs['i'], pairs['j'], pairs['v']
