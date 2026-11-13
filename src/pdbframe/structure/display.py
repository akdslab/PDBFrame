"""3D viewing in Jupyter with py3Dmol, behind ``pf.py3Dmol()`` (``pip install 'pdbframe[view]'``)."""

from __future__ import annotations

from . import io

__all__ = ['view', 'COLORS', 'STYLES']

#: Default palette: structure i (or chain i) gets COLORS[i].
COLORS = ['#4E79A7', '#F28E2B', '#59A14F', '#E15759', '#76B7B2', '#EDC948',
          '#B07AA1', '#FF9DA7', '#9C755F', '#BAB0AC']

STYLES = ('cartoon', 'stick', 'line', 'sphere')


def _per_structure(value, n: int, what: str) -> list:
    values = list(value) if isinstance(value, (list, tuple)) else [value] * n
    if len(values) != n:
        raise ValueError(f'{len(values)} {what} for {n} structures')
    return values


def view(frames, color_by: str | None = None, colors=None, style='cartoon',
         width: int = 800, height: int = 600):
    """Show one or more structures in an interactive viewer.

    Args:
        frames: a PDBFrame or a list of them, drawn in the same coordinate frame.
        color_by: ``'chain'`` (one color per chain ID), ``'structure'`` (one color
            per structure) or ``'bfactor'`` (red-white-blue over 50-90, the pLDDT
            scale of predicted models: blue is confident). Default: ``'chain'``
            for one structure, ``'structure'`` for several.
        colors: override the palette - a list with one color per structure (for
            ``'structure'``), or a ``{chain: color}`` dict (for ``'chain'``). Any
            CSS color name or hex code.
        style: ``'cartoon'``, ``'stick'``, ``'line'`` or ``'sphere'``, or a list
            with one style per structure (e.g. a cartoon plus a selection in sticks).

    Returns:
        The py3Dmol view; it renders when it is the last value in a cell, and can be
        customised further with the py3Dmol API.
    """
    try:
        import py3Dmol
    except ImportError:
        raise ImportError("view needs py3Dmol: pip install 'pdbframe[view]'") from None
    frames = list(frames) if isinstance(frames, (list, tuple)) else [frames]
    styles = _per_structure(style, len(frames), 'styles')
    bad = set(styles) - set(STYLES)
    if bad:
        raise ValueError(f'style must be one of {STYLES}, not {sorted(bad)}')
    color_by = color_by or ('chain' if len(frames) == 1 else 'structure')

    all_chains = sorted({c for f in frames for c in f.chains})
    if color_by == 'chain':
        chain_colors = {c: COLORS[i % len(COLORS)] for i, c in enumerate(all_chains)}
        chain_colors.update(colors or {})
    elif color_by == 'structure':
        structure_colors = _per_structure(colors, len(frames), 'colors') if colors \
            else [COLORS[i % len(COLORS)] for i in range(len(frames))]
    elif color_by != 'bfactor':
        raise ValueError("color_by must be 'chain', 'structure' or 'bfactor'")

    v = py3Dmol.view(width=width, height=height)
    for model, (frame, st) in enumerate(zip(frames, styles)):
        v.addModel(io.to_string(frame.df, 'pdb', name=frame.name), 'pdb')
        if color_by == 'chain':
            for chain in frame.chains:
                v.setStyle({'model': model, 'chain': chain}, {st: {'color': chain_colors[chain]}})
        elif color_by == 'structure':
            v.setStyle({'model': model}, {st: {'color': structure_colors[model]}})
        else:
            scheme = {'prop': 'b', 'gradient': 'rwb', 'min': 50, 'max': 90}
            v.setStyle({'model': model}, {st: {'colorscheme': scheme}})
    v.zoomTo()
    return v
