"""Loading structures straight from the RCSB PDB.

    >>> PDBFrame.from_rcsb('3tfk')                 # the deposited structure, mmCIF
    >>> PDBFrame.from_rcsb('3tfk', assembly=1)     # biological assembly 1
    >>> PDBFrame.from_rcsb('3tfk', save='data/')   # also keep the downloaded file
    >>> PDBFrame.from_rcsb('3tfk', fmt='fasta')    # the full sequences, no coordinates
    >>> PDBFrame.rcsb_info(['3tfk', '1ao7'])       # title, release date, resolution, ...

Structures come from https://files.rcsb.org/download, sequences from
https://www.rcsb.org/fasta/entry, entry metadata from https://data.rcsb.org. mmCIF is the default: the PDB format is not
available for large entries. The FASTA sequences are the full deposited
constructs, so they include residues that have no coordinates in the structure.
"""

from __future__ import annotations

import io as _io
import json
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser

__all__ = ['url', 'fetch', 'read', 'info']

BASE_URL = 'https://files.rcsb.org/download'
FASTA_URL = 'https://www.rcsb.org/fasta/entry'
DATA_URL = 'https://data.rcsb.org/rest/v1/core/entry'


def url(pdb_id: str, assembly: int | None = None, fmt: str = 'cif') -> str:
    """Download URL of an entry (``assembly=None``) or one of its biological assemblies."""
    pdb_id = pdb_id.strip().upper()
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        raise ValueError(f'not a PDB ID: {pdb_id!r}')
    if fmt == 'fasta':
        if assembly:
            raise ValueError('FASTA is per entry; assembly does not apply')
        return f'{FASTA_URL}/{pdb_id}'
    if fmt == 'cif':
        name = f'{pdb_id}-assembly{assembly}.cif' if assembly else f'{pdb_id}.cif'
    elif fmt == 'pdb':
        name = f'{pdb_id}.pdb{assembly}' if assembly else f'{pdb_id}.pdb'
    else:
        raise ValueError(f"fmt must be 'cif', 'pdb' or 'fasta', not {fmt!r}")
    return f'{BASE_URL}/{name}'


def fetch(pdb_id: str, assembly: int | None = None, fmt: str = 'cif',
          save: str | Path | None = None, timeout: float = 60) -> str:
    """The text of an RCSB file; `save` (a folder or a file path) also writes it to disk."""
    address = url(pdb_id, assembly, fmt)
    try:
        with urllib.request.urlopen(address, timeout=timeout) as response:
            text = response.read().decode('utf-8')
    except urllib.error.HTTPError as e:
        what = f'assembly {assembly} of ' if assembly else ''
        raise ValueError(f'RCSB has no {fmt} file for {what}{pdb_id.upper()} '
                         f'(HTTP {e.code}: {address})') from None
    if save is not None:
        save = Path(save)
        filename = address.rsplit('/', 1)[1] + ('.fasta' if fmt == 'fasta' else '')
        target = save / filename if save.is_dir() or not save.suffix else save
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(text)
    return text


def read(pdb_id: str, assembly: int | None = None, fmt: str = 'cif',
         save: str | Path | None = None):
    """An RCSB entry as a Bio.PDB Structure (see :func:`fetch` for the arguments)."""
    return _read(pdb_id, assembly, fmt, save)[0]


def _read(pdb_id: str, assembly: int | None = None, fmt: str = 'cif',
          save: str | Path | None = None):
    """:func:`read`, and the parser it used (an mmCIF parser keeps the file's dictionary)."""
    text = fetch(pdb_id, assembly, fmt, save)
    parser = MMCIFParser(QUIET=True) if fmt == 'cif' else PDBParser(QUIET=True)
    return parser.get_structure(pdb_id.upper(), _io.StringIO(text)), parser


# ------------------------------------------------------------------- metadata

def _entry(pdb_id: str, timeout: float) -> dict:
    """One entry's metadata, as RCSB's Data API gives it."""
    address = f'{DATA_URL}/{pdb_id}'
    try:
        with urllib.request.urlopen(address, timeout=timeout) as response:
            return json.loads(response.read().decode('utf-8'))
    except urllib.error.HTTPError as e:
        raise ValueError(f'RCSB has no entry {pdb_id} (HTTP {e.code}: {address})') from None


def info(pdb_ids, timeout: float = 60) -> pd.DataFrame:
    """Metadata of RCSB entries, one request per entry.

    Args:
        pdb_ids: one PDB ID or a list of them, any case.

    Returns:
        One row per entry: ``pdb_id``, ``title``, ``release_date`` (a datetime),
        ``method`` (e.g. ``'X-ray'``), ``resolution`` (Å; NaN when there is none,
        as for NMR) and ``n_entities`` (number of polymer entities).
    """
    rows = []
    for pid in [pdb_ids] if isinstance(pdb_ids, str) else list(pdb_ids):
        pid = url(pid).rsplit('/', 1)[1].removesuffix('.cif')     # checks the ID
        e = _entry(pid, timeout)
        entry = e.get('rcsb_entry_info', {})
        resolution = entry.get('resolution_combined') or [None]
        rows.append({'pdb_id': pid, 'title': e.get('struct', {}).get('title'),
                     'release_date': e.get('rcsb_accession_info', {})
                                      .get('initial_release_date'),
                     'method': entry.get('experimental_method'),
                     'resolution': resolution[0],
                     'n_entities': entry.get('polymer_entity_count')})
    out = pd.DataFrame(rows, columns=['pdb_id', 'title', 'release_date', 'method',
                                      'resolution', 'n_entities'])
    out['release_date'] = pd.to_datetime(out['release_date'], errors='coerce')
    out['resolution'] = pd.to_numeric(out['resolution'])
    return out
