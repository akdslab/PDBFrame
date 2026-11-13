"""The PDBFrame object: a protein structure as a pandas atom table.

    >>> pf = PDBFrame.read('model.pdb')
    >>> pf                                   # the atom table
    >>> pf.sequences                         # {chain: sequence}
    >>> pf.select(chains='A', atoms='ca')    # a new PDBFrame
    >>> pf.metrics.rmsd(ref), pf.metrics.lddt(ref)   # chains matched by sequence
    >>> pf.metrics.summary(ref)              # all scores at once
    >>> pf.superpose(ref).write('model_on_ref.pdb')

Every method returns a new PDBFrame (or a number / table); nothing is modified in
place. ``pf.df`` is the plain DataFrame for anything not covered here.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..sequence import SeqFrame
from . import deposited, display, io, selection
from .geometry import kabsch, transform
from .matching import chain_map, match
from .selection import ATOM_SETS

__all__ = ['PDBFrame', 'ATOM_SETS']


class PDBFrame:
    """A protein structure: one row per atom.

    Args:
        df: an atom table with at least ``chain``, ``residue_number``,
            ``residue_name``, ``atom_name``, ``x``, ``y``, ``z``. Missing optional
            columns are filled in (see :data:`pdbframe.structure.io.COLUMNS`).
        name: a label, used in printouts and written files.
        metadata: facts about the structure as a whole, kept in :attr:`metadata`
            (e.g. the RCSB entry's release date and resolution; see
            ``from_rcsb(..., metadata=True)``).
    ``observed_index`` is always the 1-based position of the residue within its
    chain in *this* frame, so it is recomputed after a selection. The numbering
    from the file is in ``residue_number`` / ``insertion_code``.

    Every method returns a new, independent PDBFrame. The new one shares column
    data with its parent until either of them changes that column (pandas
    copy-on-write), so deriving frames is cheap and changing one never affects
    the other. :meth:`copy` makes a fully separate copy up front.
    """

    #: Facts about the whole structure, ``{key: value}``; empty unless given. Derived
    #: frames (selections, transforms, ...) keep it; ``+`` starts empty.
    metadata: dict

    #: The deposited sequence behind :attr:`residues` (see :mod:`.deposited`), or None
    #: (formats without one). Internal: read it through :attr:`residues`.
    _seqres: pd.DataFrame | None = None

    def __init__(self, df: pd.DataFrame, name: str | None = None,
                 metadata: dict | None = None):
        self.df = io.normalize(df)
        self.name = name or 'structure'
        self.metadata = dict(metadata or {})
        self._seqres = None

    @classmethod
    def _wrap(cls, df: pd.DataFrame, name: str | None, metadata: dict | None = None,
              seqres: pd.DataFrame | None = None) -> 'PDBFrame':
        """A PDBFrame around a table that is already normalized - no checks, no copy.
        `seqres` is the deposited sequence (:mod:`.deposited`) of exactly these atoms."""
        frame = cls.__new__(cls)
        frame.df, frame.name = df, name or 'structure'
        frame.metadata = dict(metadata or {})
        frame._seqres = seqres
        return frame

    @classmethod
    def read(cls, path: str | Path, model: int = 0, hetero: bool = False,
             name: str | None = None) -> 'PDBFrame':
        """Read a PDB / mmCIF / CSV / FASTA file (``.gz`` too). See :func:`.io.read`.

        A FASTA file gives a frame of sequences without coordinates; see
        :meth:`from_sequences`.
        """
        path = Path(path)
        df, seqres = io.read_seqres(path, model=model, hetero=hetero)
        return cls._wrap(df, name=name or path.name.split('.')[0], seqres=seqres)

    @classmethod
    def from_rcsb(cls, pdb_id: str, assembly: int | None = None, fmt: str = 'cif',
                  model: int = 0, hetero: bool = False, save: str | Path | None = None,
                  name: str | None = None, metadata: bool = False) -> 'PDBFrame':
        """Download an entry from the RCSB PDB, e.g. ``PDBFrame.from_rcsb('3tfk')``.

        Args:
            pdb_id: the 4-character PDB ID.
            assembly: a biological assembly number (``1`` is the usual one);
                default ``None``, the deposited structure.
            fmt: ``'cif'`` (default; every entry has one), ``'pdb'``, or ``'fasta'``
                for the entry's full sequences without coordinates (chains named as
                in the structure; see :meth:`from_sequences`).
            model, hetero: as in :meth:`read`.
            save: also write the downloaded file - a folder, or a file path.
            name: defaults to the PDB ID, upper case.
            metadata: also fetch the entry's metadata (:meth:`rcsb_info`: title,
                release date, method, resolution, ...) into :attr:`metadata`; one
                more request.

        See :mod:`pdbframe.structure.rcsb`.
        """
        from . import rcsb
        name = name or pdb_id.strip().upper()
        if fmt == 'fasta':
            text = rcsb.fetch(pdb_id, assembly=assembly, fmt='fasta', save=save)
            frame = cls._wrap(io.from_sequences(io.parse_fasta(text)), name)
        else:
            structure, parser = rcsb._read(pdb_id, assembly=assembly, fmt=fmt, save=save)
            frame = cls.from_biopython(structure, model=model, hetero=hetero, name=name)
            if fmt == 'cif':
                frame._seqres = io._parsed_seqres(parser, frame.df)
        if metadata:
            frame.metadata = rcsb.info(pdb_id).iloc[0].to_dict()
        return frame

    @staticmethod
    def rcsb_info(pdb_ids, timeout: float = 60) -> pd.DataFrame:
        """RCSB metadata of entries, e.g. ``PDBFrame.rcsb_info(['3tfk', '1ao7'])``:
        ``pdb_id``, ``title``, ``release_date``, ``method``, ``resolution``,
        ``n_entities``, one row per entry. See :func:`pdbframe.structure.rcsb.info`."""
        from . import rcsb
        return rcsb.info(pdb_ids, timeout=timeout)

    @classmethod
    def from_sequences(cls, sequences, name: str | None = None) -> 'PDBFrame':
        """A PDBFrame of sequences only - no coordinates.

        One row per residue (a CA row with empty ``x``, ``y``, ``z``), numbered from 1,
        so selecting, sequences, residues, renaming and matching to a structure work
        as usual; anything that needs coordinates raises an error.

        Args:
            sequences: ``{chain: sequence}``, a list of sequences (chains A, B, ...),
                or one sequence (chain A). See :func:`.io.from_sequences`.
        """
        return cls._wrap(io.from_sequences(sequences), name or 'sequences')

    @classmethod
    def from_biopython(cls, structure, model: int = 0, hetero: bool = False,
                       name: str | None = None) -> 'PDBFrame':
        """A PDBFrame from a Bio.PDB Structure (or one Model of it).

        `model` and `hetero` work as in :meth:`read`; the name defaults to the
        structure's ID. See :func:`.io.from_biopython`.
        """
        default = structure.get_parent().id if getattr(structure, 'level', None) == 'M' \
            and structure.get_parent() is not None else structure.id
        return cls._wrap(io.from_biopython(structure, model=model, hetero=hetero),
                         name=name or str(default))

    @classmethod
    def ensure(cls, data) -> 'PDBFrame':
        """`data` as a PDBFrame: a PDBFrame, a file path, an atom table or a Bio.PDB structure."""
        if isinstance(data, PDBFrame):
            return data
        if isinstance(data, (str, Path)):
            return cls.read(data)
        if isinstance(data, pd.DataFrame):
            return cls(data)
        if getattr(data, 'level', None) in ('S', 'M'):
            return cls.from_biopython(data)
        raise TypeError(f'cannot make a PDBFrame from {type(data).__name__}')

    def _same(self, df: pd.DataFrame, seqres=False) -> 'PDBFrame':
        """A new frame for `df`, this one's atoms changed row by row (same residues, so
        nothing to recompute). The deposited sequence follows renamed chains and
        renumbered residues (:func:`.deposited.remap`), unless `seqres` gives it.

        Like every derived frame, it has this frame's class, so subclasses survive
        selecting, transforming, renaming, ...
        """
        if seqres is False:
            seqres = deposited.remap(self._seqres, self.df, df)
        return type(self)._wrap(df, self.name, getattr(self, 'metadata', None), seqres)

    def _subset(self, rows, unobserved: np.ndarray | None = None) -> 'PDBFrame':
        """A new frame of some of this one's atoms (a boolean mask or row positions).

        Only ``observed_index`` changes, since it is the position within the chain. The
        unobserved residues kept are those `unobserved` marks (one bool per deposited
        residue), or by default those inside what is kept (:func:`.deposited.subset`).
        """
        rows = np.asarray(rows)
        sub = (self.df[rows] if rows.dtype == bool else self.df.iloc[rows]).reset_index(drop=True)
        sub['observed_index'] = io.observed_index(sub)
        return type(self)._wrap(sub, self.name, getattr(self, 'metadata', None),
                                deposited.subset(self._seqres, sub, unobserved))

    # ---------------------------------------------------------------- content

    @property
    def coords(self) -> np.ndarray:
        """``(N, 3)`` array of atom coordinates."""
        return self.df[['x', 'y', 'z']].to_numpy(float)

    @property
    def has_coordinates(self) -> bool:
        """False for a frame of sequences only (from FASTA or :meth:`from_sequences`)."""
        return bool(self.df[['x', 'y', 'z']].notna().any(axis=None))

    @property
    def chains(self) -> list[str]:
        return self.df['chain'].unique().tolist()

    @property
    def residues(self) -> pd.DataFrame:
        """One row per residue: numbering, name, one-letter code, atom count, and
        whether it was ``observed``.

        Read from mmCIF, it is every residue of the deposited sequence: the residues an
        experiment did not resolve are there too, in sequence order, with ``observed``
        False, ``atoms`` 0 and no ``observed_index`` (they are not in :attr:`df`). So is
        a modified residue (a HETATM record in the chain) read without ``hetero=True``,
        as it has no atoms here either.
        ``full_index`` is each residue's position in its chain's deposited sequence
        (empty without one: other formats, ligands). A selection keeps the unobserved
        residues that belong to it (see :meth:`select`).
        """
        return deposited.residue_table(self.df, self._seqres)

    @property
    def unobserved(self) -> pd.DataFrame:
        """The rows of :attr:`residues` without atoms (none unless read from mmCIF)."""
        res = self.residues
        return res[~res['observed']].reset_index(drop=True)

    @property
    def full_sequences(self) -> dict[str, str]:
        """``{chain: sequence}`` with the unobserved residues (:attr:`residues`): the
        deposited sequences of an mmCIF file. :attr:`sequences` has the observed ones
        only; the two are the same for other formats."""
        res = self.residues
        return {c: ''.join(g['aa']) for c, g in res.groupby('chain', sort=False)}

    @property
    def sequences(self) -> dict[str, str]:
        """``{chain: one-letter sequence}`` of the residues with atoms; non-standard
        residues are ``X``. :attr:`full_sequences` adds the unobserved ones."""
        res = deposited.observed_residues(self.df)
        return {c: ''.join(g['aa']) for c, g in res.groupby('chain', sort=False)}

    def seqframe(self, chain: str, full: bool = False) -> SeqFrame:
        """One chain's sequence as a :class:`SeqFrame`, named after the chain.

        Args:
            full: the full sequence, unobserved residues included (:attr:`full_sequences`),
                not only the residues with atoms (:attr:`sequences`).
        """
        if chain not in self.chains:
            raise KeyError(f'no chain {chain!r}; chains are {self.chains}')
        return SeqFrame((self.full_sequences if full else self.sequences)[chain], name=chain)

    def seqframes(self, chains: list[str] | str | None = None,
                  full: bool = False) -> dict[str, SeqFrame]:
        """``{chain: SeqFrame}`` for the given chains (all by default), in chain order;
        `full` as in :meth:`seqframe`.

        Unpack them into :meth:`SeqFrame.align` to align the chains with each other:
        ``SeqFrame.align(*pf.seqframes(['A', 'B']).values())``.
        """
        chains = self.chains if chains is None else [chains] if isinstance(chains, str) else list(chains)
        return {c: self.seqframe(c, full) for c in chains}

    def __len__(self) -> int:
        return len(self.df)

    def _summary(self) -> str:
        chains = self.chains
        shown = ', '.join(chains[:10]) + (', ...' if len(chains) > 10 else '')
        n_res = len(self.df.drop_duplicates(['chain', 'observed_index']))
        atoms = f'{len(self)} atoms' if self.has_coordinates else 'sequences only'
        return f'PDBFrame {self.name!r}: {len(chains)} chains ({shown}), {n_res} residues, {atoms}'

    def __repr__(self) -> str:
        return f'{self._summary()}\n{self.df!r}'

    def _repr_html_(self) -> str:
        return f'<p><b>{self._summary()}</b></p>{self.df._repr_html_()}'

    # --------------------------------------------------------------- selecting

    def select(self, chains=None, residues=None, atoms=None, mask=None,
               invert: bool = False) -> 'PDBFrame':
        """Keep the atoms matching every filter given.

        Unobserved residues (:attr:`residues`) have no atoms, so `chains` and
        `residues` pick them as they pick the others; with `mask` (or `atoms` and
        `invert` together) a run of them is kept when the selection keeps the observed
        residues on both sides of it.

        Args:
            chains: a chain ID or list of them.
            residues: residue numbers (from the file) - a ``(start, end)`` tuple,
                inclusive, or a list / range of numbers. Combine with `chains` to
                pick a stretch of one chain.
            atoms: ``'all'``, ``'heavy'`` (no hydrogens), ``'backbone'`` (N, CA,
                C, O), ``'ca'``, or a list of atom names.
            mask: a boolean array, one value per atom.
            invert: keep everything *except* the selection.
        """
        keep = selection.select_mask(self.df, chains=chains, residues=residues, atoms=atoms,
                                     mask=mask)
        unobserved = None
        if self._seqres is not None and mask is None and not (invert and atoms is not None):
            picked = selection.select_mask(self._seqres, chains=chains, residues=residues)
            unobserved = ~picked if invert else picked
        return self._subset(~keep if invert else keep, unobserved)

    def near(self, other: 'PDBFrame', cutoff: float = 5.0) -> 'PDBFrame':
        """Whole residues of this structure with any atom within `cutoff` Å of `other`."""
        other = PDBFrame.ensure(other)
        for f in (self, other):
            io.require_coordinates(f.df, f.name, 'measure distances')
        return self.select(mask=selection.near_mask(self.df, other.df, cutoff))

    def contacts(self, other: 'PDBFrame | None' = None, cutoff: float = 5.0) -> pd.DataFrame:
        """Residue pairs in contact, with their closest atom-atom distance.

        With `other`, pairs between this structure and `other`; without it, pairs
        between different chains of this structure (its interfaces).
        """
        other = None if other is None else PDBFrame.ensure(other)
        for f in (self, other or self):
            io.require_coordinates(f.df, f.name, 'measure distances')
        return selection.contacts(self.df, None if other is None else other.df, cutoff)

    # ----------------------------------------------------------------- editing

    def rename_chains(self, mapping: 'dict[str, str] | PDBFrame') -> 'PDBFrame':
        """Rename chains.

        Args:
            mapping: ``{'A': 'H', 'B': 'L'}`` (unlisted chains keep their ID), or a
                reference PDBFrame - each chain then takes the ID of its counterpart
                in the reference, found by sequence (:meth:`map_chains`).
        """
        if isinstance(mapping, PDBFrame):
            mapping = self.map_chains(mapping)
        new = [mapping.get(c, c) for c in self.chains]
        if len(set(new)) < len(new):
            raise ValueError(f'renaming {self.chains} -> {new} would merge chains')
        # chains are never merged, so observed_index stays valid
        return self._same(self.df.assign(chain=self.df['chain'].map(dict(zip(self.chains, new)))))

    def renumber(self, start: int = 1) -> 'PDBFrame':
        """Number each chain's residues consecutively from `start`; drops insertion codes.

        Unobserved residues (:attr:`residues`) count too, so the numbers skip over a gap.
        """
        res = self.residues
        res['new_number'] = res.groupby('chain', sort=False).cumcount() + start
        res['new_insertion_code'] = ''
        seen = res[res['observed']].astype({'observed_index': int})
        new = pd.MultiIndex.from_frame(self.df[['chain', 'observed_index']])
        number = seen.set_index(['chain', 'observed_index'])['new_number'].reindex(new)
        seqres = (None if self._seqres is None else deposited.renumbered(
            self._seqres, res[deposited.KEYS + ['new_number', 'new_insertion_code']]))
        return self._same(self.df.assign(residue_number=number.to_numpy(), insertion_code=''),
                          seqres)

    def transform(self, rotation: np.ndarray, translation: np.ndarray) -> 'PDBFrame':
        """Apply ``xyz @ rotation.T + translation`` to every atom."""
        xyz = transform(self.coords, np.asarray(rotation), np.asarray(translation))
        return self._same(self.df.assign(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2]))

    def __add__(self, other: 'PDBFrame') -> 'PDBFrame':
        """Combine two structures into one; their chain IDs must not overlap."""
        other = PDBFrame.ensure(other)
        shared = set(self.chains) & set(other.chains)
        if shared:
            raise ValueError(f'both structures have chains {sorted(shared)}; '
                             'use rename_chains() on one first')
        frame = type(self)(pd.concat([self.df, other.df]), name=f'{self.name}+{other.name}')
        frame._seqres = deposited.concat(self._seqres, other._seqres)
        return frame

    def copy(self) -> 'PDBFrame':
        """A fully separate copy, with its own data from the start.

        Rarely needed: every PDBFrame method already returns an independent frame.
        """
        return self._same(self.df.copy(deep=True),
                          None if self._seqres is None else self._seqres.copy(deep=True))

    # --------------------------------------------------------------- comparing

    def map_chains(self, reference: 'PDBFrame', min_identity: float = 0.7) -> dict[str, str]:
        """Which chain of `reference` each of this structure's chains corresponds to.

        Found by sequence, so chain IDs need not agree; see :func:`.matching.chain_map`.
        Returns ``{my_chain: reference_chain}``; chains without a counterpart are absent.
        """
        return chain_map(self.df, PDBFrame.ensure(reference).df, min_identity)

    def match(self, other: 'PDBFrame', by: str = 'sequence', chains='auto',
              min_identity: float = 0.7) -> 'tuple[PDBFrame, PDBFrame]':
        """The atoms this structure and `other` have in common, in the same order.

        Chains are paired by sequence unless told otherwise (`chains`: ``'auto'``,
        a list of my chains, ``'same'`` for equal IDs, or ``{my_chain: other_chain}``).
        Residues are paired by sequence alignment (``by='sequence'``) or residue
        number (``by='number'``), atoms by name. See :func:`.matching.match`.

        Returns:
            ``(mine, theirs)`` - two PDBFrames, row i of one paired with row i of
            the other.
        """
        other = PDBFrame.ensure(other)
        r1, r2 = match(self.df, other.df, by=by, chains=chains, min_identity=min_identity)
        if len(r1) == 0:
            raise ValueError(f'{self.name} and {other.name} have no atoms in common')
        return self._subset(r1), other._subset(r2)

    def _paired(self, other, atoms, match_kw) -> 'tuple[pd.DataFrame, pd.DataFrame]':
        """Matched atom tables of self and `other`, restricted to `atoms`."""
        mine, theirs = self.match(other, **match_kw)
        for f in (mine, theirs):
            io.require_coordinates(f.df, f.name, 'be superposed or scored')
        keep = selection.atom_mask(mine.df, atoms) & selection.atom_mask(theirs.df, atoms)
        if not keep.any():
            raise ValueError(f'no {atoms!r} atoms in common')
        return mine.df[keep], theirs.df[keep]

    def superpose(self, other: 'PDBFrame', atoms='ca', **match_kw) -> 'PDBFrame':
        """This structure moved onto `other`, by the best fit of their common `atoms`.

        `match_kw` (``by``, ``chains``, ``min_identity``) go to :meth:`match`, so
        the fit can use only some chains while the whole structure moves.
        """
        a, b = self._paired(other, atoms, match_kw)
        R, t = kabsch(a[['x', 'y', 'z']].to_numpy(float), b[['x', 'y', 'z']].to_numpy(float))
        return self.transform(R, t)

    @property
    def metrics(self):
        """Scores of this structure against a reference: ``pf.metrics.rmsd(ref)``, ...

        See :class:`.metrics.Metrics`.
        """
        from .metrics import Metrics
        return Metrics(self)

    # ------------------------------------------------------------------ output

    def write(self, path: str | Path) -> Path:
        """Save as PDB / mmCIF / CSV / FASTA, picked from the extension.

        mmCIF keeps the unobserved residues (:attr:`residues`) in its sequence tables,
        so reading the file back gives them again; PDB, CSV and FASTA leave them out,
        with a warning. FASTA headers look like ``>name|Chain A``, which :meth:`read`
        turns back into the same chains.
        """
        return io.write(self.df, path, name=self.name, seqres=self._seqres)

    def to_string(self, fmt: str = 'pdb') -> str:
        """The structure as PDB (``'pdb'``) or mmCIF (``'cif'``) text (see :meth:`write`)."""
        return io.to_string(self.df, fmt, name=self.name, seqres=self._seqres)

    def to_fasta(self, path: str | Path | None = None) -> str:
        """Chain sequences as FASTA text (``>name_chain``), saved to `path` if given."""
        text = ''.join(f'>{self.name}_{c}\n{s}\n' for c, s in self.sequences.items())
        if path is not None:
            Path(path).write_text(text)
        return text

    def to_biopython(self):
        """The structure as a Bio.PDB Structure."""
        return io.to_biopython(self.df, self.name)

    def view(self, regions=None, sequence=False, **kwargs):
        """3D view with each chain's sequence linked to it; see :func:`.viewer.view`.

        `regions` marks residues (``{name: residues}``). `sequence` picks the chains
        whose sequences show at first (``True`` for all; none by default). Other
        options include ``compare=`` (a second structure, residues linked by
        sequence), ``values=`` (color by per-residue values such as lDDT or
        pLDDT), ``interfaces=``, ``surface=`` and ``axes=``. The **?** button in
        the viewer lists the mouse controls.
        """
        from .viewer import view
        return view(self, regions=regions, sequence=sequence, **kwargs)

    def py3Dmol(self, *others: 'PDBFrame', **kwargs):
        """py3Dmol view of this and any `others` overlaid; see :func:`.display.view`.

        Use it to compare several structures, color by structure or B-factor, or
        keep customizing the returned py3Dmol view.
        """
        return display.view([self, *others], **kwargs)
