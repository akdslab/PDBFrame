---
name: analyze-structure
description: Answer questions about a protein structure file (PDB / mmCIF) with the pdbframe library - chains, sequences, residues, selections, contacts, interfaces, distances. Use when asked what a structure contains, which residues touch something, or how far apart parts are.
---

# Analyze a structure with pdbframe

`pip install git+https://github.com/akdslab/PDBFrame` (in its repo: `uv run python ...`).

```python
from pdbframe import PDBFrame
pf = PDBFrame.read('x.cif')                     # .pdb .cif .csv, .gz ok; hetero=True keeps ligands
pf = PDBFrame.from_biopython(structure)         # or from a Bio.PDB Structure / Model
pf = PDBFrame.from_rcsb('3tfk')                 # download; assembly=1, fmt='pdb', save='dir/'
pf = PDBFrame.from_rcsb('3tfk', metadata=True)  # pf.metadata: title, release_date, method, resolution
PDBFrame.rcsb_info(['3tfk', '1ao7'])            # that metadata as a table, one row per entry
seqs = PDBFrame.from_rcsb('3tfk', fmt='fasta')  # or .read('x.fasta'): sequences, no coordinates
pf.chains, pf.sequences                         # ['A', ...], {chain: one-letter sequence}
pf.residues                                     # one row per residue
pf.seqframes()                                  # {chain: SeqFrame} for aligning (align-sequences skill)
pf.df                                           # one row per atom (plain pandas)
```

**Columns that matter.** `residue_number` + `insertion_code` = the file's numbering (use for
reporting and `select(residues=...)`). `observed_index` = 1-based position in the chain *within
this frame*, recomputed after every selection - never report it as a residue number.
`bfactor` holds pLDDT in predicted models.

**Selecting** - every call returns a new, independent PDBFrame:
```python
pf.select(chains='A', residues=(50, 85), atoms='backbone')   # atoms: all|heavy|backbone|ca|[names]
pf.select(chains='C', invert=True)
pocket = pf.near(pf.select(chains='C'), cutoff=4)             # whole residues within 4 Å
pf.contacts(cutoff=4)            # residue pairs between chains: chain_1, residue_number_1, ..., distance
pf.contacts(other_frame, cutoff=5)
```

**Distances**: `pf.select(...).coords` gives an (N, 3) numpy array; for residue-residue minimum
distance use `contacts` between two selections, or `scipy.spatial.distance.cdist` on coords.

**Pitfalls**
- A FASTA frame (`has_coordinates` False) selects and matches like a structure, but anything
  needing coordinates raises. Matching it to the structure shows unmodeled residues.
- A selection's `sequences` concatenates only the selected residues; a scattered region is not
  a contiguous stretch. Report residues from `.residues` (name + `residue_number`).
- `pf + other` needs distinct chain IDs (`rename_chains` first).
- PDB output allows only 1-character chain IDs; write `.cif` otherwise.

Every method and option: `pdbframe docs pdbframe` (prints `docs/pdbframe.md`).
