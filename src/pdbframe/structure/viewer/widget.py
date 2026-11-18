"""Python side of the viewer: PDBFrames and options in, a self-contained HTML widget out.

The page itself lives in ``viewer.js`` and ``viewer.css`` next to this file; this
module only prepares the data they read.

    >>> pocket = ref.near(ref.select(chains='C'), cutoff=4)
    >>> ref.view(regions={'pocket': pocket, 'helix': ('A', 50, 85)}, sequence=True)
"""

from __future__ import annotations

import json
import string
import uuid
from importlib.resources import files
from pathlib import Path

import numpy as np
import pandas as pd

from .. import io
from ..display import COLORS
from ..matching import match

__all__ = ['view', 'view_many', 'Viewer', 'REGION_COLORS', 'HIGHLIGHT_COLOR', 'PICK_COLOR',
           'COMPARE_COLOR', 'STRUCTURE_COLORS', 'ARROW_COLORS']

#: Default region colors - none of them is in the chain palette.
REGION_COLORS = ['#FFD700', '#FF1493', '#00CED1', '#7CFC00', '#8A2BE2', '#FF4500']
#: Color of hovered and selected residues.
HIGHLIGHT_COLOR = '#00FF00'
#: Color of picked residues / chains and the distance lines between them.
PICK_COLOR = '#FF8C00'
#: Default color of a single structure passed as ``compare``.
COMPARE_COLOR = '#A9A9A9'
#: Default colors when several structures are compared - one each, none in the chain palette.
STRUCTURE_COLORS = ['#A9A9A9', '#C7A27C', '#9AB0C9', '#B5C99A', '#C99AB5', '#9AC9C0', '#D4B483',
                    '#A89AC9']
#: Default arrow colors, one per arrow in turn.
ARROW_COLORS = ['#000000', '#D62728', '#1F77B4', '#2CA02C', '#9467BD', '#8C564B']
#: Atoms kept by ``compact=True``: enough for cartoons and backbone measurements.
COMPACT_ATOMS = ('N', 'CA', 'C', 'O')
# value colors: low -> high
_VALUE_SCALE = [(0.0, (214, 39, 40)), (0.5, (255, 255, 255)), (1.0, (31, 119, 180))]
# value ranges that have a natural scale
_KNOWN_RANGES = {'plddt': (50.0, 90.0), 'bfactor': (50.0, 90.0), 'lddt': (0.0, 1.0)}

_JS_URL = 'https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js'
_STYLES = ('cartoon', 'stick', 'line', 'sphere')
# one-character chain IDs for the PDB text handed to 3Dmol.js, so chain IDs of
# any length work; the real IDs are kept for the labels
_CHAIN_IDS = string.ascii_uppercase + string.ascii_lowercase + string.digits
_RES_KEYS = ['chain', 'residue_number', 'insertion_code']

_PAGE = """<div id="pfv-{uid}" class="pfv">
<style>{css}</style>
<div class="bar"></div>
<div class="panel"><div class="tracks"></div><div class="legend"></div></div>
<div class="info"></div>
<div class="picks"></div>
<div class="view"><div class="axes"></div></div>
</div>
<script>
(function() {{
var D = {data};
{js}
}})();
</script>
"""


class Viewer:
    """The widget: shows itself in a notebook; ``save()`` writes a standalone page."""

    def __init__(self, html: str):
        self.html = html

    def _repr_html_(self) -> str:
        return self.html

    def save(self, path: str | Path) -> Path:
        """Write a standalone HTML page (opening it needs internet access for 3Dmol.js)."""
        path = Path(path)
        path.write_text('<!DOCTYPE html>\n<html><head><meta charset="utf-8"></head>'
                        f'<body>{self.html}</body></html>\n')
        return path


# --------------------------------------------------------------------- inputs

def _residues_of(frame, spec) -> list[tuple[str, int]]:
    """(chain, observed_index) of the residues a region spec covers.

    A spec is a chain ID, a ``(chain, start, end)`` residue-number range
    (inclusive), a PDBFrame (e.g. a selection of `frame`), or a list of these.
    """
    res = frame.residues
    if isinstance(spec, str):
        hit = res[res['chain'] == spec]
    elif isinstance(spec, tuple) and len(spec) == 3 and isinstance(spec[0], str):
        chain, start, end = spec
        hit = res[(res['chain'] == chain) & res['residue_number'].between(start, end)]
    elif hasattr(spec, 'residues') and hasattr(spec, 'df'):
        # observed_index is recomputed per frame, so pair by the file's numbering
        hit = res.merge(spec.residues[_RES_KEYS], on=_RES_KEYS)
    elif isinstance(spec, (list, tuple)):
        return list(dict.fromkeys(k for s in spec for k in _residues_of(frame, s)))
    else:
        raise TypeError(f'cannot read a region from {type(spec).__name__}: use a chain '
                        'ID, a (chain, start, end) tuple, a PDBFrame or a list of these')
    return list(zip(hit['chain'], hit['observed_index'].astype(int)))


def _interfaces(frame, interfaces, cutoff: float) -> dict[str, list[tuple[str, int]]]:
    """``{name: residues}`` for chain-chain interfaces: residues of either chain
    within `cutoff` Å of the other, as (chain, observed_index)."""
    pairs = frame.contacts(cutoff=cutoff)
    if interfaces is True:
        wanted = list(dict.fromkeys(zip(pairs['chain_1'], pairs['chain_2'])))
    else:
        wanted = [tuple(p) for p in interfaces]
    out = {}
    for a, b in wanted:
        hit = pairs[((pairs['chain_1'] == a) & (pairs['chain_2'] == b)) |
                    ((pairs['chain_1'] == b) & (pairs['chain_2'] == a))]
        if hit.empty:
            raise ValueError(f'chains {a!r} and {b!r} have no residues within {cutoff} Å')
        side1 = hit[['chain_1', 'residue_number_1', 'insertion_code_1']].set_axis(_RES_KEYS, axis=1)
        side2 = hit[['chain_2', 'residue_number_2', 'insertion_code_2']].set_axis(_RES_KEYS, axis=1)
        sides = pd.concat([side1, side2]).drop_duplicates()
        keys = frame.residues.merge(sides, on=_RES_KEYS)
        out[f'{a}–{b} interface'] = list(zip(keys['chain'], keys['observed_index'].astype(int)))
    return out


def _value_colors(frame, values, column, value_range):
    """Per-residue values -> (colors and values aligned with frame.residues, legend)."""
    res = frame.residues
    if isinstance(values, str):
        if values not in ('bfactor', 'plddt'):
            raise ValueError("values must be a table, 'bfactor' or 'plddt'")
        table = frame.metrics.plddt(per_residue=True)
        column = 'plddt'
        name = values
    else:
        table = pd.DataFrame(values)
        missing = [c for c in _RES_KEYS if c not in table]
        if missing:
            raise KeyError(f'values table needs columns {_RES_KEYS}; missing {missing}')
        if column is None:
            known = [c for c in ('lddt', 'plddt', 'value') if c in table]
            numeric = [c for c in table.select_dtypes('number').columns
                       if c not in ('residue_number', 'observed_index', 'atoms')]
            if not (known or numeric):
                raise ValueError('values table has no numeric column to color by')
            column = (known or numeric)[-1]
        name = column
    merged = res.merge(table[_RES_KEYS + [column]], on=_RES_KEYS, how='left')
    v = merged[column].to_numpy(float)
    if np.isnan(v).all():
        raise ValueError(f'values table matches no residues of {frame.name!r}: its chain, '
                         'residue_number and insertion_code must be that structure\'s')
    lo, hi = value_range or _KNOWN_RANGES.get(name, (np.nanmin(v), np.nanmax(v)))
    t = np.clip((v - lo) / ((hi - lo) or 1.0), 0, 1)
    colors = []
    for x in t:
        if np.isnan(x):
            colors.append(None)
            continue
        for (t0, c0), (t1, c1) in zip(_VALUE_SCALE[:-1], _VALUE_SCALE[1:]):
            if x <= t1:
                f = (x - t0) / (t1 - t0)
                colors.append('#%02x%02x%02x' % tuple(round(a + f * (b - a)) for a, b in zip(c0, c1)))
                break
    legend = {'name': name, 'lo': float(lo), 'hi': float(hi),
              'colors': ['#%02x%02x%02x' % c for _, c in _VALUE_SCALE]}
    return colors, [None if np.isnan(x) else round(float(x), 3) for x in v], legend


def _shown(frame, sequence) -> list[str]:
    if sequence is True:
        return frame.chains
    if sequence is False or sequence is None:
        return []
    shown = [sequence] if isinstance(sequence, str) else list(sequence)
    missing = set(shown) - set(frame.chains)
    if missing:
        raise KeyError(f'sequence names chains that are not in {frame.name!r}: {sorted(missing)}')
    return shown


def _model(frame, palette, shown, values=None, name=None, style=None, compact=False):
    """One structure's data: PDB text and per-chain residue rows."""
    inner = dict(zip(frame.chains, _CHAIN_IDS))
    # 3Dmol.js gets residues numbered by observed_index, which is unique within a
    # chain even when the file uses insertion codes; labels show the real numbers
    df = frame.df[frame.df['atom_name'].isin(COMPACT_ATOMS)] if compact else frame.df
    df = df.copy()
    df['chain'] = df['chain'].map(inner)
    df['residue_number'] = df['observed_index']
    df['insertion_code'] = ''
    res = frame.residues
    vcol, vval = values if values else ([None] * len(res), [None] * len(res))
    res = res.assign(value_color=vcol, value=vval)
    chains = [{'id': inner[c], 'label': c, 'color': palette[c], 'shown': c in shown,
               'residues': [[int(r.observed_index), f'{r.residue_number}{r.insertion_code}',
                             r.residue_name, r.aa, r.value_color, r.value]
                            for r in g.itertuples(index=False)]}
              for c, g in res.groupby('chain', sort=False)]
    return {'name': name or frame.name, 'pdb': io.to_string(df, 'pdb', name=frame.name), 'chains': chains,
            'style': style}, inner


def _residue_pairs(frame, other) -> list[tuple[str, int, str, int]]:
    """(chain, observed_index) pairs of `frame` and `other`, paired by sequence."""
    r0, r1 = match(frame.df, other.df)
    a = frame.df.iloc[r0][['chain', 'observed_index']].to_numpy()
    b = other.df.iloc[r1][['chain', 'observed_index']].to_numpy()
    pairs = {(c0, int(i0), c1, int(i1)) for (c0, i0), (c1, i1) in zip(a, b)}
    return sorted(pairs, key=lambda p: (p[2], p[3]))


def _links(pairs_by_model, inners) -> dict:
    """Residue links between model 0 and every other model, both directions."""
    links = {}
    for m, pairs in pairs_by_model.items():
        for c0, i0, c1, i1 in pairs:
            k0, k1 = f'0:{inners[0][c0]}:{i0}', f'{m}:{inners[m][c1]}:{i1}'
            links.setdefault(k0, []).append(k1)
            links.setdefault(k1, []).append(k0)
    return links


def _layout(frames, inners, pairs_by_model) -> list[dict]:
    """Sequence panel rows: each lead chain, then its partners in the other structures,
    with gap cells inserted so paired residues share a column.

    Returns groups ``{'rows': [{'chain': 'm:id', 'cells': [index or None, ...]}]}``;
    chains of other structures paired with no lead chain get a group of their own.
    """
    groups, placed = [], set()
    for c0, g0 in frames[0].residues.groupby('chain', sort=False):
        idx0 = g0['observed_index'].astype(int).tolist()
        col = {i: n for n, i in enumerate(idx0)}               # observed_index -> lead position
        partners = []                                          # (row key, at, before)
        for m, pairs in pairs_by_model.items():
            lead_pos = {(c1, i1): col[i0] for cc0, i0, c1, i1 in pairs if cc0 == c0}
            res_m = frames[m].residues
            for c1 in dict.fromkeys(c for c, _ in lead_pos):
                at, before, last = {}, {}, -1                  # before[p]: inserted ahead of lead position p
                for i in res_m.loc[res_m['chain'] == c1, 'observed_index'].astype(int):
                    p = lead_pos.get((c1, i))
                    if p is None:
                        before.setdefault(last + 1, []).append(i)
                    else:
                        at[p], last = i, p
                partners.append((f'{m}:{inners[m][c1]}', at, before))
                placed.add((m, c1))
        width = {p: max([len(b.get(p, [])) for _, _, b in partners] or [0])
                 for p in range(len(idx0) + 1)}

        def cells(at, before):
            out = []
            for p in range(len(idx0) + 1):
                ins = before.get(p, [])
                out += ins + [None] * (width[p] - len(ins))
                if p < len(idx0):
                    out.append(at.get(p))
            return out

        rows = [{'chain': f'0:{inners[0][c0]}', 'cells': cells(dict(enumerate(idx0)), {})}]
        rows += [{'chain': key, 'cells': cells(at, before)} for key, at, before in partners]
        groups.append({'rows': rows})
    for m in range(1, len(frames)):
        for c, g in frames[m].residues.groupby('chain', sort=False):
            if (m, c) not in placed:
                groups.append({'rows': [{'chain': f'{m}:{inners[m][c]}',
                                         'cells': g['observed_index'].astype(int).tolist()}]})
    return groups


def _fit(lead, other, superpose):
    """`other` superposed on `lead`: on all chains, or on those paired with the given lead chains."""
    if superpose is True or superpose == 'all':
        return other.superpose(lead)
    fit = [superpose] if isinstance(superpose, str) else list(superpose)
    missing = set(fit) - set(lead.chains)
    if missing:
        raise KeyError(f'superpose names chains that are not in {lead.name!r}: {sorted(missing)}')
    mine = [c for c, target in other.map_chains(lead).items() if target in fit]
    if not mine:
        raise ValueError(f'{other.name!r} has no chains matching {fit} to superpose on')
    return other.superpose(lead, chains=mine)


def _per_structure(value, n: int, what: str) -> list:
    """`value` for each of `n` structures: a list is taken as one entry per structure,
    anything else applies to the first structure only."""
    if isinstance(value, list):
        if len(value) != n:
            raise ValueError(f'{what}: {len(value)} entries for {n} structures')
        return value
    return [value] + [None] * (n - 1)


def _rotation(spec) -> tuple[list, list]:
    """An orientation - a 3x3 rotation, or (rotation, translation) - as JSON lists."""
    if isinstance(spec, (tuple, list)) and len(spec) == 2 and np.ndim(spec[0]) == 2:
        R, t = np.asarray(spec[0], float), np.asarray(spec[1], float)
    else:
        R, t = np.asarray(spec, float), np.zeros(3)
    if R.shape != (3, 3) or t.shape != (3,):
        raise ValueError('an orientation is a 3x3 rotation matrix or a (rotation, translation) pair')
    if not np.allclose(R @ R.T, np.eye(3), atol=1e-4) or np.linalg.det(R) < 0:
        raise ValueError('an orientation needs a proper rotation matrix (orthonormal, det +1)')
    return R.tolist(), t.tolist()


def _opacity(opacity) -> float:
    """How solid the structure is drawn, 0 to 1."""
    value = float(opacity)
    if not 0.0 <= value <= 1.0:
        raise ValueError(f'opacity must be between 0 and 1, not {opacity}')
    return value


def _arrows(arrows) -> list[dict]:
    """Arrows as JSON: each a dict or a ``(start, end[, color[, label]])`` tuple."""
    if arrows is None:
        return []
    if isinstance(arrows, dict) and 'start' in arrows:
        arrows = [arrows]
    out = []
    for i, a in enumerate(arrows):
        if not isinstance(a, dict):
            a = dict(zip(('start', 'end', 'color', 'label'), a))
        unknown = set(a) - {'start', 'end', 'color', 'label', 'radius', 'group', 'head', 'info'}
        if unknown:
            raise KeyError(f'arrow {i}: unknown keys {sorted(unknown)}; use start, end, '
                           'color, label, radius, group, head, info')
        try:
            start, end = (np.asarray(a[k], float).reshape(3) for k in ('start', 'end'))
        except (KeyError, ValueError):
            raise ValueError(f'arrow {i}: needs a start and an end point, each 3 numbers') from None
        if not (np.isfinite(start).all() and np.isfinite(end).all()) or np.allclose(start, end):
            raise ValueError(f'arrow {i}: start and end must be finite, different points')
        out.append({'start': start.tolist(), 'end': end.tolist(),
                    'color': a.get('color') or ARROW_COLORS[i % len(ARROW_COLORS)],
                    'label': None if a.get('label') is None else str(a['label']),
                    'radius': float(a.get('radius', 0.4)),
                    'group': str(a.get('group') or 'Arrows'),
                    'head': bool(a.get('head', True)),
                    'info': None if a.get('info') is None else str(a['info'])})
    return out


def _spheres(spheres) -> list[dict]:
    """Spheres as JSON: each a dict or a ``(center[, radius[, color[, label]]])`` tuple."""
    if spheres is None:
        return []
    if isinstance(spheres, dict) and 'center' in spheres:
        spheres = [spheres]
    out = []
    for i, sp in enumerate(spheres):
        if not isinstance(sp, dict):
            sp = dict(zip(('center', 'radius', 'color', 'label'), sp))
        unknown = set(sp) - {'center', 'radius', 'color', 'label', 'group', 'info'}
        if unknown:
            raise KeyError(f'sphere {i}: unknown keys {sorted(unknown)}; use center, '
                           'radius, color, label, group, info')
        try:
            center = np.asarray(sp['center'], float).reshape(3)
        except (KeyError, ValueError):
            raise ValueError(f'sphere {i}: needs a center of 3 numbers') from None
        if not np.isfinite(center).all():
            raise ValueError(f'sphere {i}: center must be finite')
        out.append({'center': center.tolist(),
                    'color': sp.get('color') or ARROW_COLORS[i % len(ARROW_COLORS)],
                    'label': None if sp.get('label') is None else str(sp['label']),
                    'radius': float(sp.get('radius', 0.6)),
                    'group': str(sp.get('group') or 'Spheres'),
                    'info': None if sp.get('info') is None else str(sp['info'])})
    return out


# ---------------------------------------------------------------------- widget

def view(frame, regions=None, sequence=False, compare=None, values=None,
         value_column: str | None = None, value_range: tuple | None = None,
         interfaces=None, interface_cutoff: float = 5.0,
         region_colors: dict | None = None, region_styles: dict | None = None,
         colors=None, compare_color: str | None = None,
         style='cartoon', surface: bool = False, ligands: bool = True,
         axes: bool = True, width: int = 900, height: int = 500,
         names: list | None = None, hidden_regions=None, region_groups: dict | None = None,
         orientations: dict | None = None, orientation_center: str = 'structure',
         compact: bool = False, superpose=False, arrows=None, spheres=None,
         opacity: float = 1.0, orientation=None, sequence_only: bool = False) -> Viewer:
    """3D view of a structure, with each chain's sequence linked to it.

    Args:
        frame: the PDBFrame to show; it leads when there are several structures.
        regions: ``{name: spec}`` of residues to mark. A spec is a chain ID, a
            ``(chain, start, end)`` residue-number range (inclusive), a PDBFrame
            such as ``frame.select(...)`` or ``frame.near(...)``, or a list of these.
            With several structures, a list of such dicts gives one per structure,
            resolved against that structure; regions with the same name are one region.
        sequence: which chains' sequences show at first - ``False`` for none (the
            default), ``True`` for all, or a list of chain IDs of `frame` (the
            compared chains paired with them show too). The *Sequences* menu
            changes it.
        compare: other structures to draw in the same space, e.g. models superposed
            on `frame`: one PDBFrame, a list, or ``{name: PDBFrame}``. Each one's
            residues are paired with `frame`'s by sequence, so hovering or selecting
            a residue lights up its partners in all of them, and their sequences
            line up under `frame`'s.
        values: per-residue values to color by - a table with ``chain``,
            ``residue_number``, ``insertion_code`` and a value column (such as
            ``model.metrics.lddt(ref, per_residue=True)``), or ``'plddt'`` /
            ``'bfactor'`` for the B-factor column. Red is low, blue is high. A list
            gives one per structure (``None`` for none).
        value_column: which column of `values` to use (default: ``lddt``,
            ``plddt`` or ``value`` if present, else the last numeric column).
        value_range: ``(low, high)`` of the color scale (default: 0-1 for lDDT,
            50-90 for pLDDT, else the data's range).
        interfaces: mark chain-chain interfaces of `frame` as regions - ``True``
            for every pair of chains in contact, or a list of pairs like ``[('A', 'C')]``.
        interface_cutoff: residues of two chains within this distance (Å) form
            their interface.
        region_colors: ``{name: color}`` to override region colors.
        region_styles: ``{name: style}`` to also draw a region in another style on
            top of the main one, e.g. ``{'pocket': 'stick'}``.
        colors: ``{chain: color}`` to override `frame`'s chain colors. A list gives
            one entry per structure: a ``{chain: color}`` dict, one color for the
            whole structure, or ``None``.
        compare_color: one color for every compared structure. Default: grey for a
            single one, a different color each when there are several.
        style: 3D style - ``'cartoon'``, ``'stick'``, ``'line'`` or ``'sphere'``, or
            a list with one per structure.
        surface: start with the molecular surface shown (a checkbox toggles it).
        ligands: draw ligands / ions as sticks (they are only in the frame when
            it was read with ``hetero=True``).
        axes: show x, y, z arrows in the bottom-left corner.
        width, height: size of the 3D panel in pixels.
        names: names for all structures, `frame` first (default: their own names,
            or the keys when `compare` is a dict).
        hidden_regions: names of regions (or region groups) to start switched off;
            click a region in the legend to switch it on.
        region_groups: ``{group: [region names]}`` - a legend entry that switches
            several regions on and off together.
        orientations: ``{label: rotation}`` - a button each that turns the view so
            the given frame is seen down -z with y up. A rotation is a 3x3 matrix or
            ``(R, t)`` as returned by superposition; the coordinates are not changed.
        orientation_center: ``'structure'`` keeps the structure centered (only the
            rotation is used); ``'frame'`` also centers on the point the
            transform maps to the origin.
        orientation: the orientation the view opens in, and *Reset view* returns to:
            a rotation as for `orientations`, or the label of one of them. Only the
            camera turns; the coordinates are not changed. Default: as the
            coordinates are, looking down -z.
        compact: send only backbone atoms (N, CA, C, O) to the page - smaller pages
            for many structures; side chains and ligands are then not shown.
        superpose: fit each compared structure onto `frame` (by C-alpha atoms)
            before drawing - ``True`` or ``'all'`` to fit on all chains, or a chain
            ID / list of `frame`'s chain IDs to fit on only those (each structure's
            matching chains are found by sequence). Only the drawn copies move; the frames passed in
            are not changed. Default ``False``: draw them as they are.
        arrows: arrows to draw in `frame`'s coordinates, e.g. vectors or axes - a
            list of ``(start, end)``, ``(start, end, color, label)`` or dicts with
            ``start``, ``end`` and optionally ``color``, ``label`` (text at the
            tip), ``radius`` (Å, default 0.4) and ``group``. Each group (default
            ``'Arrows'``) gets a checkbox that shows and hides its arrows.
            ``head=False`` draws a plain cylinder, for something that has no
            direction; ``info`` is text shown when the arrow is clicked.
        spheres: balls to draw in `frame`'s coordinates, e.g. the nodes of a graph
            whose edges are `arrows` - a list of ``(center,)``,
            ``(center, radius, color, label)`` or dicts with ``center`` and
            optionally ``color``, ``label``, ``radius`` (A, default 0.6),
            ``group`` and ``info``. Groups share the arrows' checkboxes, so an
            edge group and its nodes can carry the same name and switch together.
        opacity: how solid the structure is drawn, 0 to 1 (default 1). Below 1 it
            keeps its shape and colors but lets drawn shapes show through.
        sequence_only: open the sequences only (every chain, unless `sequence`
            names some), with their regions and values; a *Show 3D* button opens
            the 3D view under them.
    """
    others = compare if isinstance(compare, (list, tuple)) else \
        list(compare.values()) if isinstance(compare, dict) else \
        [] if compare is None else [compare]
    frames = [frame] + list(others)
    n = len(frames)
    if superpose is not False and superpose is not None and n > 1:
        frames = [frame] + [_fit(frame, f, superpose) for f in frames[1:]]
    if names is None:
        names = [frame.name] + (list(compare.keys()) if isinstance(compare, dict)
                                else [f.name for f in others])
    if len(names) != n:
        raise ValueError(f'{len(names)} names for {n} structures')
    given = [str(x) for x in names]
    # keep names distinct: repeated ones get their position, e.g. 'model (2)'
    names = [f'{x} ({i + 1})' if given.count(x) > 1 else x for i, x in enumerate(given)]

    styles = _per_structure(style, n, 'style') if isinstance(style, list) else [style] * n
    region_styles = region_styles or {}
    for st in [*styles, *region_styles.values()]:
        if st not in _STYLES:
            raise ValueError(f'style must be one of {_STYLES}, not {st!r}')
    region_dicts = [r or {} for r in _per_structure(regions, n, 'regions')]
    region_names = list(dict.fromkeys(k for r in region_dicts for k in r))
    unknown = set(region_styles) - set(region_names)
    if unknown:
        raise KeyError(f'region_styles names regions that are not given: {sorted(unknown)}')
    if orientation_center not in ('structure', 'frame'):
        raise ValueError("orientation_center must be 'structure' or 'frame'")
    if isinstance(orientation, str):
        if orientation not in (orientations or {}):
            raise ValueError(f'orientation {orientation!r} is not one of the orientations: '
                             f'{list(orientations or {})}')
        orientation = orientations[orientation]
    for f in frames:
        if len(f.chains) > len(_CHAIN_IDS):
            raise ValueError(f'view shows at most {len(_CHAIN_IDS)} chains per structure')

    # colors: the lead's chain palette; compared structures one color each
    color_specs = _per_structure(colors, n, 'colors') if isinstance(colors, list) else [colors] + [None] * (n - 1)
    palettes = []
    for m, (f, spec) in enumerate(zip(frames, color_specs)):
        if m == 0:
            base = {c: COLORS[i % len(COLORS)] for i, c in enumerate(sorted(f.chains))}
        else:
            whole = compare_color or (COMPARE_COLOR if n == 2 else
                                      STRUCTURE_COLORS[(m - 1) % len(STRUCTURE_COLORS)])
            base = {c: whole for c in f.chains}
        if isinstance(spec, str):
            base = {c: spec for c in f.chains}
        elif isinstance(spec, dict):
            base.update(spec)
        palettes.append(base)

    value_specs = _per_structure(values, n, 'values')
    value_data = [_value_colors(f, v, value_column, value_range) if v is not None else None
                  for f, v in zip(frames, value_specs)]

    if sequence_only and (sequence is False or sequence is None):
        sequence = True
    shown = _shown(frame, sequence)
    pairs_by_model = {m: _residue_pairs(frame, frames[m]) for m in range(1, n)}
    models, inners = [], []
    for m, f in enumerate(frames):
        if m == 0:
            show = shown
        else:          # `sequence` names the lead's chains; show their partners too
            partner_of = {c1: c0 for c0, _, c1, _ in pairs_by_model[m]}
            show = [c for c in f.chains if partner_of.get(c) in shown]
        model, inner = _model(f, palettes[m], show,
                              values=value_data[m][:2] if value_data[m] else None,
                              name=names[m], style=styles[m], compact=compact)
        models.append(model)
        inners.append(inner)
    links = _links(pairs_by_model, inners)
    layout = _layout(frames, inners, pairs_by_model)

    # regions: same name across structures -> one region
    resolved = {}
    for m, (f, rd) in enumerate(zip(frames, region_dicts)):
        for name, spec in rd.items():
            keys = _residues_of(f, spec)
            if not keys:
                raise ValueError(f'region {name!r} covers no residues of {f.name!r}')
            resolved.setdefault(name, []).extend(f'{m}:{inners[m][c]}:{k}' for c, k in keys)
    if interfaces:
        for name, keys in _interfaces(frame, interfaces, interface_cutoff).items():
            resolved.setdefault(name, []).extend(f'0:{inners[0][c]}:{k}' for c, k in keys)
    marked = [{'name': str(name), 'extra': region_styles.get(name), 'keys': keys,
               'color': (region_colors or {}).get(name, REGION_COLORS[i % len(REGION_COLORS)])}
              for i, (name, keys) in enumerate(resolved.items())]

    groups = {str(g): [str(x) for x in members] for g, members in (region_groups or {}).items()}
    missing = {x for members in groups.values() for x in members} - set(resolved)
    if missing:
        raise KeyError(f'region_groups names regions that are not given: {sorted(missing)}')
    hidden = [str(x) for x in ([hidden_regions] if isinstance(hidden_regions, str)
                               else hidden_regions or [])]
    missing = set(hidden) - set(resolved) - set(groups)
    if missing:
        raise KeyError(f'hidden_regions names regions that are not given: {sorted(missing)}')
    hidden = list(dict.fromkeys(x for h in hidden for x in groups.get(h, [h])))

    legend_values = []
    for vd in value_data:
        if vd and vd[2] not in legend_values:
            legend_values.append(vd[2])
    uid = uuid.uuid4().hex[:12]
    data = {'uid': uid, 'js_url': _JS_URL, 'models': models, 'links': links, 'layout': layout,
            'regions': marked, 'region_groups': groups, 'hidden_regions': hidden,
            'values': legend_values[0] if legend_values else None, 'value_scales': legend_values,
            'orientations': {str(k): _rotation(v) for k, v in (orientations or {}).items()},
            'orientation_center': orientation_center,
            'orientation': None if orientation is None else _rotation(orientation),
            'style': styles[0], 'highlight': HIGHLIGHT_COLOR, 'pick': PICK_COLOR,
            'surface': bool(surface), 'ligands': bool(ligands) and not compact,
            'axes': bool(axes), 'width': width, 'height': height,
            'sequence_only': bool(sequence_only),
            'arrows': _arrows(arrows), 'spheres': _spheres(spheres),
            'opacity': _opacity(opacity)}
    here = files(__package__)
    html = _PAGE.format(uid=uid,
                        css=here.joinpath('viewer.css').read_text().replace('#ROOT', f'#pfv-{uid}'),
                        js=here.joinpath('viewer.js').read_text(),
                        data=json.dumps(data).replace('</', '<\\/'))
    return Viewer(html)


def view_many(frames, names: list | None = None, **kwargs) -> Viewer:
    """Several structures in one view: ``pdbframe.view([ref, m1, m2], names=[...])``.

    The first structure leads (the others are paired with it); options are those
    of :func:`view`, with per-structure lists in the same order as `frames`.
    """
    if not isinstance(frames, (list, tuple)):
        frames = [frames]
    if not frames:
        raise ValueError('view_many needs at least one structure')
    return view(frames[0], compare=list(frames[1:]) or None, names=names, **kwargs)
