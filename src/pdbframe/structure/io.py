"""Reading and writing structure files as atom tables.

Every structure is one pandas DataFrame, one row per atom, with the columns in
:data:`COLUMNS`. Coordinates live only in ``x``, ``y``, ``z``.

Supported formats, picked from the file extension (``.gz`` allowed for reading):
PDB (``.pdb``, ``.ent``), mmCIF (``.cif``), CSV (``.csv``, a saved atom table) and
FASTA (``.fasta``, ``.fa``, ``.faa``, ``.fas``).

A FASTA file gives a table without coordinates: one row per residue (a CA row
with ``x``, ``y``, ``z`` empty), so sequences can be selected and compared like
structures. See :func:`from_sequences`.
"""

from __future__ import annotations

import gzip
import io
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import MMCIFIO, PDBIO, MMCIFParser, PDBParser
from Bio.Data.IUPACData import protein_letters_1to3
from Bio.PDB.Polypeptide import protein_letters_3to1_extended
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.StructureBuilder import StructureBuilder

from . import deposited

__all__ = ['COLUMNS', 'read', 'read_seqres', 'seqres_table',
           'from_biopython', 'from_sequences', 'parse_fasta', 'fasta_text', 'write',
           'normalize', 'observed_index', 'require_coordinates', 'to_biopython', 'to_string']

#: The atom table, in column order.
#: ``residue_number`` + ``insertion_code`` are the numbering from the file;
#: ``observed_index`` is the 1-based position of the residue within its chain.
COLUMNS = ['chain', 'residue_number', 'insertion_code', 'observed_index',
           'residue_name', 'aa', 'atom_name', 'element', 'altloc',
           'x', 'y', 'z', 'occupancy', 'bfactor', 'hetero']

# Defaults for columns a hand-built table may leave out.
_DEFAULTS = {'insertion_code': '', 'altloc': '', 'occupancy': 1.0, 'bfactor': 0.0,
             'hetero': False}

_FORMATS = {'.pdb': 'pdb', '.ent': 'pdb', '.pdb1': 'pdb', '.cif': 'cif',
            '.mmcif': 'cif', '.csv': 'csv', '.fasta': 'fasta', '.fa': 'fasta',
            '.faa': 'fasta', '.fas': 'fasta'}


def _format(path: Path) -> str:
    suffixes = [s.lower() for s in path.suffixes]
    if suffixes and suffixes[-1] == '.gz':
        suffixes = suffixes[:-1]
    if not suffixes or suffixes[-1] not in _FORMATS:
        raise ValueError(f'unknown structure format: {path.name} '
                         f'(expected one of {sorted(_FORMATS)})')
    return _FORMATS[suffixes[-1]]


def read(path: str | Path, model: int = 0, hetero: bool = False) -> pd.DataFrame:
    """Read a PDB / mmCIF / CSV / FASTA file into an atom table.

    Args:
        path: the file. A trailing ``.gz`` is decompressed on the fly.
        model: which model to read (0-based), for files holding several.
        hetero: keep HETATM records (ligands, ions). Water is always dropped.
    """
    return read_seqres(path, model=model, hetero=hetero)[0]


def read_seqres(path: str | Path, model: int = 0, hetero: bool = False
                ) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    """:func:`read`, and the file's deposited sequences (:func:`seqres_table`).

    Only mmCIF files have them (their ``_pdbx_poly_seq_scheme``); for other formats,
    and mmCIF files without that table, the second value is None.
    """
    path = Path(path)
    fmt = _format(path)
    if fmt == 'csv':
        return normalize(pd.read_csv(path, keep_default_na=False,
                                     na_values={'x': [''], 'y': [''], 'z': ['']})), None

    opener = gzip.open if path.suffix.lower() == '.gz' else open
    if fmt == 'fasta':
        with opener(path, 'rt') as handle:
            return from_sequences(parse_fasta(handle.read())), None
    parser = PDBParser(QUIET=True) if fmt == 'pdb' else MMCIFParser(QUIET=True)
    with opener(path, 'rt') as handle:
        structure = parser.get_structure(path.stem, handle)
    df = from_biopython(structure, model=model, hetero=hetero)
    return df, (_parsed_seqres(parser, df) if fmt == 'cif' else None)


def _parsed_seqres(parser, df: pd.DataFrame) -> pd.DataFrame | None:
    """The :func:`seqres_table` of the mmCIF file `parser` (an MMCIFParser) just read."""
    mmcif = getattr(parser, '_mmcif_dict', None)
    return seqres_table(mmcif, df) if mmcif else None


def _text_list(value) -> list[str]:
    return [value] if isinstance(value, str) else list(value)


def seqres_table(mmcif: dict, df: pd.DataFrame | None = None) -> pd.DataFrame | None:
    """Every chain's deposited sequence, from an mmCIF dictionary's ``_pdbx_poly_seq_scheme``.

    A structure file has atoms only for the residues that were observed; this table
    has every residue of the sequence that was deposited, with its place in the file's
    numbering, so the two can be matched residue by residue (see :mod:`.deposited`).

    Args:
        mmcif: the file as a ``Bio.PDB.MMCIF2Dict`` dictionary.
        df: the file's atom table. ``observed`` is then True for the residues the file
            lists as observed that `df` has (a modified residue read without
            ``hetero=True`` has no atoms there), and chains `df` lacks are left out.
            Without it, ``observed`` is what the file lists.

    Returns:
        One row per residue (columns :data:`.deposited.COLUMNS`), chains by author chain
        ID as in the atom table, in chain and sequence order, or None if the file has no
        ``_pdbx_poly_seq_scheme``. Where two residues share a position
        (micro-heterogeneity) the first is kept.
    """
    key = '_pdbx_poly_seq_scheme.'
    fields = ['pdb_strand_id', 'seq_id', 'mon_id', 'pdb_seq_num']
    if any(key + f not in mmcif for f in fields):
        return None
    raw = pd.DataFrame({f: _text_list(mmcif[key + f]) for f in fields})
    n = len(raw)
    for optional in ['pdb_ins_code', 'pdb_mon_id', 'auth_seq_num']:
        raw[optional] = _text_list(mmcif[key + optional]) if key + optional in mmcif else ['.'] * n
    raw = raw[pd.to_numeric(raw['pdb_seq_num'], errors='coerce').notna()]
    table = pd.DataFrame({
        'chain': raw['pdb_strand_id'].astype(str),
        'full_index': raw['seq_id'].astype(int),
        'residue_name': raw['mon_id'].astype(str),
        'residue_number': raw['pdb_seq_num'].astype(int),
        'insertion_code': raw['pdb_ins_code'].replace({'.': '', '?': ''}).astype(str),
    }).drop_duplicates(['chain', 'full_index']).reset_index(drop=True)
    table['aa'] = table['residue_name'].map(protein_letters_3to1_extended).fillna('X')
    first = raw.drop_duplicates(['pdb_strand_id', 'seq_id'])
    table['observed'] = ((first['pdb_mon_id'] != '?')
                         & (first['auth_seq_num'] != '?')).to_numpy()
    if df is not None:
        present = pd.MultiIndex.from_frame(df[deposited.KEYS].drop_duplicates())
        table['observed'] &= pd.MultiIndex.from_frame(table[deposited.KEYS]).isin(present)
        order = {c: i for i, c in enumerate(dict.fromkeys(df['chain']))}
        table = table[table['chain'].isin(order)]
        table = table.assign(_order=table['chain'].map(order)).sort_values(
            ['_order', 'full_index'], kind='stable').drop(columns='_order')
    return table[deposited.COLUMNS].reset_index(drop=True)


def from_biopython(structure, model: int = 0, hetero: bool = False) -> pd.DataFrame:
    """Atom table of a Bio.PDB Structure (or a single Model of one).

    Args:
        structure: a ``Bio.PDB.Structure.Structure``, or a ``Model`` taken from one.
        model: which model to use (0-based) when given a whole structure.
        hetero: keep HETATM records (ligands, ions). Water is always dropped.
    """
    level = getattr(structure, 'level', None)
    if level == 'M':
        models = [structure]
        model = 0
    elif level == 'S':
        models = list(structure)
    else:
        raise TypeError(f'expected a Bio.PDB Structure or Model, not {type(structure).__name__}')
    if not models:
        raise ValueError(f'structure {structure.id!r} holds no models')
    if model >= len(models):
        raise IndexError(f'structure has {len(models)} models, asked for model {model}')

    rows = []
    for chain in models[model]:
        for residue in chain:
            flag, number, icode = residue.id
            if flag == 'W' or (flag != ' ' and not hetero):
                continue
            name = residue.get_resname()
            for atom in residue:
                x, y, z = atom.coord
                rows.append((chain.id, number, icode.strip(), name, atom.get_name(),
                             atom.element, atom.get_altloc().strip(), x, y, z,
                             atom.get_occupancy(), atom.get_bfactor(), flag != ' '))
    df = pd.DataFrame(rows, columns=['chain', 'residue_number', 'insertion_code',
                                     'residue_name', 'atom_name', 'element', 'altloc',
                                     'x', 'y', 'z', 'occupancy', 'bfactor', 'hetero'])
    return normalize(df)


# '>3TFK_1|Chains A, C[auth B]|name|organism' (RCSB) or '>name|Chain A' (fasta_text)
_CHAINS_FIELD = re.compile(r'^([^|]*)\|Chains? ([^|]+)(?:\||$)')
_ENTITY = re.compile(r'_(\d+)$')


def parse_fasta(text: str) -> dict[str, str]:
    """``{chain: sequence}`` from FASTA text.

    Headers with a chains field - RCSB's ``>3TFK_1|Chains A, B|...`` or this
    package's ``>name|Chain A`` - give every listed chain its own copy of the
    sequence, under its author chain ID (``C[auth B]`` -> ``B``), so the chains
    match the structure's. RCSB entities come in entity order. Any other header
    names one chain by its first word.
    """
    records = []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith('>'):
            records.append([line[1:].strip(), []])
        elif line and records:
            records[-1][1].append(line)
    if not records:
        raise ValueError('no FASTA records found')
    entries = []                                   # (sort key, chains, sequence)
    for n, (head, lines) in enumerate(records):
        seq = ''.join(lines).replace('-', '').replace('*', '').upper()
        m = _CHAINS_FIELD.match(head)
        if m:
            chains = [re.sub(r'.*\[auth (\w+)\]', r'\1', c.strip()) for c in m[2].split(',')]
            entity = _ENTITY.search(m[1])
            key = int(entity[1]) if entity else n
        else:
            chains, key = [head.split()[0] if head else f'seq{n + 1}'], n
        entries.append((key, chains, seq))
    out = {}
    for _, chains, seq in sorted(entries, key=lambda e: e[0]):
        for c in chains:
            if c in out:
                raise ValueError(f'FASTA names chain {c!r} twice')
            out[c] = seq
    return out


def from_sequences(sequences) -> pd.DataFrame:
    """An atom table without coordinates: one CA row per residue, ``x/y/z`` empty.

    Args:
        sequences: ``{chain: sequence}``, a list of sequences (chains A, B, ...),
            or one sequence (chain A). One-letter codes; ``X`` and unknown letters
            become residue ``UNK``.
    """
    if isinstance(sequences, str):
        sequences = {'A': sequences}
    elif not isinstance(sequences, dict):
        sequences = dict(zip('ABCDEFGHIJKLMNOPQRSTUVWXYZ', sequences))
    rows = []
    for chain, seq in sequences.items():
        for i, letter in enumerate(str(seq).strip().upper(), 1):
            rows.append((str(chain), i, protein_letters_1to3.get(letter, 'UNK').upper()))
    if not rows:
        raise ValueError('no residues in the given sequences')
    df = pd.DataFrame(rows, columns=['chain', 'residue_number', 'residue_name'])
    df = df.assign(atom_name='CA', element='C', x=np.nan, y=np.nan, z=np.nan)
    return normalize(df)


def fasta_text(df: pd.DataFrame, name: str = 'structure') -> str:
    """Chain sequences as FASTA, with RCSB-style headers that :func:`parse_fasta` reads back."""
    res = df.drop_duplicates(['chain', 'observed_index'])
    return ''.join(f'>{name}|Chain {c}\n{"".join(g["aa"])}\n'
                   for c, g in res.groupby('chain', sort=False))


def normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Fill in derived / optional columns and put them in :data:`COLUMNS` order.

    Needs at least ``chain``, ``residue_number``, ``residue_name``, ``atom_name``
    and ``x``, ``y``, ``z``. Extra columns are kept at the end.
    """
    required = ['chain', 'residue_number', 'residue_name', 'atom_name', 'x', 'y', 'z']
    missing = [c for c in required if c not in df]
    if missing:
        raise KeyError(f'atom table is missing columns: {missing}')
    # a new table object; its columns are shared with `df` until changed (copy-on-write)
    df = df.copy(deep=False)
    for col, value in _DEFAULTS.items():
        if col not in df:
            df[col] = value
    df['chain'] = df['chain'].astype(str)
    df['insertion_code'] = df['insertion_code'].fillna('').astype(str)
    df['altloc'] = df['altloc'].fillna('').astype(str)
    df['residue_number'] = df['residue_number'].astype(int)
    if 'element' not in df:
        df['element'] = df['atom_name'].str.lstrip('0123456789').str[0]
    df['aa'] = df['residue_name'].map(protein_letters_3to1_extended).fillna('X')
    df['observed_index'] = observed_index(df)
    return df[COLUMNS + [c for c in df.columns if c not in COLUMNS]].reset_index(drop=True)


def observed_index(df: pd.DataFrame) -> np.ndarray:
    """1-based position of each atom's residue within its chain, in file order.

    Vectorized: number every distinct (chain, number, insertion code) in order of
    first appearance, then rank those numbers within each chain.
    """
    group = df.groupby(['chain', 'residue_number', 'insertion_code'], sort=False).ngroup()
    return group.groupby(df['chain'].to_numpy(), sort=False).rank(method='dense').astype(int).to_numpy()


def require_coordinates(df: pd.DataFrame, name: str, what: str) -> None:
    """Raise a clear error when a table has no coordinates (e.g. read from FASTA)."""
    if df[['x', 'y', 'z']].isna().all(axis=None):
        raise ValueError(f'{name!r} has no coordinates (it holds sequences only, e.g. from '
                         f'FASTA), so it cannot {what}')


def to_biopython(df: pd.DataFrame, name: str = 'structure'):
    """Build a Bio.PDB Structure (one model) from an atom table."""
    require_coordinates(df, name, 'be written as a structure or shown in 3D')
    builder = StructureBuilder()
    builder.init_structure(name)
    builder.init_model(0)
    chain = residue = None
    for i, row in enumerate(df.itertuples(index=False), 1):
        if row.chain != chain:
            builder.init_seg(' ')
            builder.init_chain(row.chain)
            chain, residue = row.chain, None
        res_id = (row.residue_number, row.insertion_code or ' ')
        if res_id != residue:
            flag = f'H_{row.residue_name}' if row.hetero else ' '
            builder.init_residue(row.residue_name, flag, *res_id)
            residue = res_id
        builder.init_atom(row.atom_name, np.array([row.x, row.y, row.z], dtype='f'),
                          float(row.bfactor), float(row.occupancy), row.altloc or ' ',
                          row.atom_name, i, row.element)
    return builder.get_structure()


def _writer(fmt: str, df: pd.DataFrame, name: str, seqres: pd.DataFrame | None = None):
    w = PDBIO() if fmt == 'pdb' else MMCIFIO()
    w.set_structure(to_biopython(df, name))
    if fmt == 'cif' and seqres is not None:
        buffer = io.StringIO()
        w.save(buffer)
        buffer.seek(0)
        w = MMCIFIO()
        w.set_dict(deposited.cif_tables(MMCIF2Dict(buffer), seqres))
    return w


def _warn_unobserved(fmt: str, name: str, seqres: pd.DataFrame | None) -> None:
    """Warn that a PDB / CSV / FASTA file cannot hold the unobserved residues."""
    if seqres is None or fmt == 'cif' or seqres['observed'].all():
        return
    unseen = seqres[~seqres['observed']].groupby('chain', sort=False).size()
    listed = ', '.join(f'{c} ({n})' for c, n in unseen.items())
    warnings.warn(f'{name}: {fmt.upper()} output leaves out the unobserved residues (in '
                  f'the deposited sequence, without atoms) of chains {listed}; write .cif '
                  'to keep them', stacklevel=3)


def to_string(df: pd.DataFrame, fmt: str = 'pdb', name: str = 'structure',
              seqres: pd.DataFrame | None = None) -> str:
    """The atom table as PDB or mmCIF text; mmCIF also lists the deposited sequence
    `seqres` (:mod:`.deposited`), unobserved residues included."""
    _warn_unobserved(fmt, name, seqres)
    buffer = io.StringIO()
    _writer(fmt, df, name, seqres).save(buffer)
    return buffer.getvalue()


def write(df: pd.DataFrame, path: str | Path, name: str = 'structure',
          seqres: pd.DataFrame | None = None) -> Path:
    """Write an atom table to PDB / mmCIF / CSV / FASTA, picked from the extension.

    mmCIF also gets the deposited sequence `seqres` (:mod:`.deposited`), unobserved
    residues included; the other formats cannot hold them (with a warning).
    """
    path = Path(path)
    fmt = _format(path)
    _warn_unobserved(fmt, name, seqres)
    if fmt == 'csv':
        df.to_csv(path, index=False)
    elif fmt == 'fasta':
        path.write_text(fasta_text(df, name))
    else:
        _writer(fmt, df, name, seqres).save(str(path))
    return path
