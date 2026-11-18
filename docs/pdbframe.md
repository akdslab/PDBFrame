# PDBFrame reference

A `PDBFrame` is a protein structure as a pandas table, one row per atom. Every method
returns a new frame (or a number or table); nothing is changed in place.

```python
from pdbframe import PDBFrame

pf = PDBFrame.read('notebooks/data/3tfk_af2.pdb')    # or PDBFrame.from_rcsb('3tfk')
pf                       # one-line summary, then the atom table
```

## Building one

| | |
|---|---|
| `PDBFrame.read(path, model=0, hetero=False, name=None)` | PDB, mmCIF, CSV or FASTA, `.gz` too. `hetero=True` also reads ligands, ions and water |
| `PDBFrame.from_rcsb(pdb_id, assembly=None, fmt='cif', save=None, metadata=False)` | download from the RCSB. `assembly=1` gives the biological assembly; `fmt='fasta'` gives the full deposited sequences; `save=` keeps the file; `metadata=True` also fetches the entry's metadata into `pf.metadata` |
| `PDBFrame.rcsb_info(pdb_ids)` | RCSB metadata, one row per entry: `title`, `release_date`, `method`, `resolution`, `n_entities` |
| `PDBFrame.from_sequences(seqs)` | sequences only, no coordinates; `seqs` is `{chain: seq}`, a list, or one string |
| `PDBFrame.from_biopython(structure)` | a Bio.PDB Structure or Model |
| `PDBFrame(df, name=None)` | an atom table you built yourself |
| `PDBFrame.ensure(x)` | any of the above inputs, as a PDBFrame |

## The atom table

`pf.df` has one row per atom:

| column | |
|---|---|
| `chain` | chain ID |
| `residue_number`, `insertion_code` | residue numbering as in the file (`105`, `'A'`) |
| `observed_index` | position of the residue in its chain in *this* frame: 1, 2, 3, ... It is recomputed after every selection |
| `residue_name`, `aa` | `'GLU'`, `'E'` (`'X'` for non-standard residues) |
| `atom_name`, `element`, `altloc` | atom |
| `x`, `y`, `z` | coordinates (empty for sequences only) |
| `occupancy`, `bfactor` | `bfactor` holds the pLDDT in predicted models |
| `hetero` | True for HETATM records |

Extra columns you add are kept through selections and other operations.

## Contents

| | gives |
|---|---|
| `pf.name` | label used in printouts and written files |
| `pf.metadata` | facts about the whole structure, a dict (e.g. from `from_rcsb(..., metadata=True)`: title, release date, method, resolution). Kept by selections and other derived frames; `+` starts empty |
| `pf.chains` | chain IDs, in order |
| `pf.sequences` | `{chain: sequence}` |
| `pf.full_sequences` | `{chain: sequence}` with the unobserved residues: the deposited sequences of an mmCIF file (the same as `pf.sequences` for other formats) |
| `pf.residues` | one row per residue: `chain`, `observed_index`, `residue_number`, `insertion_code`, `residue_name`, `aa`, `atoms`, `observed`, `full_index`. From mmCIF it has every residue of the deposited sequence, the unobserved ones too (below) |
| `pf.unobserved` | the rows of `pf.residues` without atoms |
| `pf.coords` | `(N, 3)` numpy array |
| `pf.has_coordinates` | False for a sequences-only frame |
| `pf.seqframe('A', full=False)` | one chain's sequence as a SeqFrame (below); `full=True` with the unobserved residues |
| `pf.seqframes(chains=None, full=False)` | `{chain: SeqFrame}` for those chains (all by default) |
| `len(pf)` | number of atoms |

### Unobserved residues

An experimental structure has atoms only for the residues that were resolved; the
others (a disordered loop, a chain's ends) are missing from the file's atoms. An mmCIF
file still lists them in its deposited sequence (`_pdbx_poly_seq_scheme`), and
`pf.residues` has them, in sequence order, with `observed=False`, `atoms=0` and no
`observed_index`. `pf.df` never has them, so nothing that uses coordinates sees them.

```python
pf = PDBFrame.read('9rup.cif')
pf.sequences['b']           # ...YFSETQFSNS...  residues with atoms
pf.full_sequences['b']      # ...YFSETQRNKGNFPGRFSGRQFSNS...  the deposited sequence
pf.unobserved               # the residues without atoms
```

They follow the frame. `rename_chains`, `renumber` (which counts them, so numbers skip
the gap) and `+` carry them. A selection keeps the ones that belong to it: `chains=`
and `residues=` pick them like any residue; a `mask` (also `near`, `match`) keeps a run
of them when it keeps the observed residues on both sides. `pf.write('out.cif')` writes
them into the file's sequence tables, so reading it back gives the same residues; PDB,
CSV and FASTA cannot hold them and warn. A modified residue (a HETATM record in the
chain) read without `hetero=True` is listed the same way, as it has no atoms either.

## Selecting

| | gives |
|---|---|
| `pf.select(chains=, residues=, atoms=, mask=, invert=False)` | atoms matching every filter given |
| `pf.near(other, cutoff=5.0)` | whole residues of `pf` with any atom within `cutoff` Å of `other` |
| `pf.contacts(other=None, cutoff=5.0)` | residue pairs in contact, with their closest atom-atom distance. Without `other`: pairs between different chains of `pf` |

`select` filters:
- `chains='A'` or `['A', 'B']`.
- `residues=(10, 20)`: a range of the file's residue numbers, inclusive. A list or range
  of numbers also works.
- `atoms='all'`, `'heavy'`, `'backbone'` (N, CA, C, O), `'ca'`, or a list of atom names.
- `mask=`: a boolean array, one value per atom.
- `invert=True`: everything *except* the selection.

```python
pf.select(chains='A', residues=(50, 86), atoms='ca')
pf.select(chains=['D', 'E']).near(pf.select(chains='C'), 4.0)
```

## Editing

| | gives |
|---|---|
| `pf.rename_chains({'A': 'H'})` | chains renamed; or pass another frame to take its IDs, matched by sequence |
| `pf.renumber(start=1)` | each chain numbered consecutively, unobserved residues counted; insertion codes dropped |
| `pf.transform(R, t)` | every atom moved: `xyz @ R.T + t` |
| `pf + other` | both structures in one frame; chain IDs must not overlap |
| `pf.copy()` | a fully separate copy (rarely needed: frames are already independent) |

## Comparing two structures

| | gives |
|---|---|
| `pf.map_chains(ref, min_identity=0.7)` | `{my_chain: ref_chain}`, found by sequence, so chain IDs need not agree |
| `pf.match(ref, by='sequence', chains='auto', min_identity=0.7)` | `(mine, theirs)`: the atoms both have, row i paired with row i |
| `pf.superpose(ref, atoms='ca')` | `pf` moved onto `ref` by the best fit of their common atoms |

`match` and `superpose` options:
- `by='sequence'` aligns each chain's sequence, so numbering may differ and residues may
  be missing. `by='number'` pairs equal residue numbers.
- `chains='auto'` finds counterparts by sequence. Also: a list of my chains, `'same'`
  (equal IDs), or `{my_chain: ref_chain}`.
- Mismatched residues (mutations) pair backbone atoms only.
- Passing `chains=` to `superpose` fits on those chains only, while the whole structure
  moves.

## Scores: `pf.metrics`

All take the reference first. Chains are matched as in `match`, and the `match` options
can be passed too.

| | gives |
|---|---|
| `pf.metrics.rmsd(ref, atoms='ca', superpose=True)` | RMSD in Å over the common atoms, after superposing unless `superpose=False` |
| `pf.metrics.tm_score(ref, superpose=True)` | TM-score (0-1), normalised by the reference's length |
| `pf.metrics.lddt(ref, atoms='heavy', cutoff=15.0, per_residue=False)` | lDDT (0-1); no superposition needed. `per_residue=True` gives a table |
| `pf.metrics.dockq(ref, side1, side2=None)` | DockQ of the interface between two chains or chain groups: `fnat`, `irmsd`, `lrmsd`, `dockq` |
| `pf.metrics.plddt(per_residue=False)` | mean B-factor, i.e. the pLDDT of predicted models |
| `pf.metrics.summary(ref)` | RMSD (Cα), TM-score, lDDT, pLDDT and the chain mapping in one Series |

## Viewing

`pf.view(...)` opens an interactive 3D view in a notebook, with each chain's sequence
linked to it. Hover or click a residue to see it in both. The **?** button lists the
mouse controls, and `.save('view.html')` writes a standalone page.

| option | |
|---|---|
| `regions={name: spec}` | residues to mark. A spec is a chain ID, a `(chain, start, end)` range, a frame such as `pf.select(...)`, or a list of these. Each region can be switched on and off in the legend |
| `hidden_regions=[...]`, `region_groups={group: [names]}`, `region_colors=`, `region_styles=` | regions that start off; legend entries that switch several at once; colors; an extra style per region |
| `sequence=True` or `['A', 'B']` | chains whose sequences show at first |
| `sequence_only=True` | the sequences only (every chain, unless `sequence=` names some), with their regions and values; *Show 3D* opens the structure under them |
| `compare=other` (a list, or `{name: frame}`) | more structures in the same view, residues linked by sequence |
| `superpose=True` or chain IDs | fit the compared structures onto `pf` first (only the drawn copies move) |
| `values=` | color by per-residue values: a table such as `model.metrics.lddt(ref, per_residue=True)`, or `'plddt'` / `'bfactor'`. Also `value_column=`, `value_range=` |
| `interfaces=True` or `[('A', 'C')]` | mark chain-chain interfaces as regions (`interface_cutoff=`) |
| `style=` | `'cartoon'`, `'stick'`, `'line'` or `'sphere'` |
| `colors=`, `compare_color=` | chain colors; one color for the compared structures |
| `surface=True`, `ligands=True`, `axes=True` | molecular surface, ligands as sticks, x/y/z arrows |
| `arrows=[(start, end, color, label), ...]` | vectors drawn in `pf`'s coordinates, one checkbox per `group`. `head=False` for a plain cylinder, `info=` for text shown when it is clicked |
| `spheres=[(center, radius, color, label), ...]` | balls in `pf`'s coordinates, sharing the arrows' group checkboxes - the nodes to the arrows' edges |
| `opacity=0.4` | how solid the structure is drawn, so shapes inside it show through |
| `orientations={label: R}` | buttons that turn the camera to saved orientations |
| `orientation=R` or a label | the orientation the view opens in (and *Reset view* returns to). Only the camera turns; coordinates are unchanged |
| `names=`, `width=`, `height=`, `compact=True` | structure names; panel size; backbone only (smaller pages) |

`pf.py3Dmol(*others, ...)` gives a py3Dmol view instead, for further customizing.

## Writing

| | gives |
|---|---|
| `pf.write('out.pdb')` | PDB, mmCIF, CSV or FASTA, picked from the extension. Only mmCIF keeps the unobserved residues; the others warn |
| `pf.to_string('pdb')` / `'cif'` | the file text |
| `pf.to_fasta(path=None)` | FASTA text (`>name_chain`), saved if `path` is given |
| `pf.to_biopython()` | a Bio.PDB Structure |

To convert a structure file to FASTA: `PDBFrame.read('in.pdb').write('out.fasta')`. Its
headers (`>name|Chain A`) read back into the same chains; `to_fasta()`'s `>name_A` headers
read back as chains named `name_A`.

## SeqFrame

`pf.seqframe('A')` gives one chain's sequence as a SeqFrame, for aligning and comparing
sequences; `pf.seqframes()` gives every chain's. To align chains with each other:
`SeqFrame.align(*pf.seqframes(['C', 'D']).values())`. See [seqframe.md](seqframe.md).

## Subclassing

Frames derived from a subclass (by `select`, `transform`, `copy`, `+`, ...) are of that
subclass. That lets a library extend PDBFrame with its own methods; tcrstruct's
`TCRFrame` does this.
