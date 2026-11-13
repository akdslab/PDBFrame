---
name: compare-structures
description: Score predicted protein models against a reference structure (RMSD, TM-score, lDDT, pLDDT, DockQ per interface) with the pdbframe library, for one model or a folder of them. Use when evaluating AlphaFold/Boltz/other predictions or comparing two structures of the same protein.
---

# Compare models with a reference

**Many models** - run the bundled script (needs pdbframe installed):
```bash
python scripts/score.py reference.cif models/ [more_models ...] -o scores.csv
```
One row per model: `rmsd`, `tm_score`, `lddt`, `plddt`, `chain_map`, `dockq_<A>-<B>` for every
chain pair in contact in the reference, and `error` if a model failed.

**One model, in Python**
```python
from pdbframe import PDBFrame
ref, model = PDBFrame.read('ref.cif'), PDBFrame.read('model.pdb')
model.metrics.summary(ref)                  # rmsd (CA, superposed), tm_score, lddt, plddt, chain_map
model.metrics.dockq(ref, 'A', 'C')          # sides use the REFERENCE's chain IDs
model.metrics.lddt(ref, per_residue=True)   # where the model is wrong
model.metrics.rmsd(ref, chains=['H'])       # restrict to model chains
moved = model.superpose(ref, chains=['H'])  # fit on some chains, move everything
```
Pairing is automatic: chains by sequence (IDs need not match; check `model.map_chains(ref)`),
residues by alignment (numbering may differ), atoms by name. Force it with
`chains={'H': 'A'}` or `chains='same'`.

**Reading the numbers**
- TM-score: > 0.5 same fold; < 0.17 random. Normalised by the reference length.
- lDDT: 0-1, superposition-free, local accuracy.
- DockQ: < 0.23 incorrect, 0.23-0.49 acceptable, 0.49-0.80 medium, >= 0.80 high.
- pLDDT (model's own confidence): > 90 very high, 70-90 confident, 50-70 low, < 50 very low.
- A large global RMSD with good per-chain RMSD means domains/chains are placed differently -
  check `{c: model.metrics.rmsd(ref, chains=[c]) for c in model.chains}`.

TM-score here keeps the sequence pairing, so it can sit slightly below TM-align's value.

Every method and option: `pdbframe docs pdbframe` (prints `docs/pdbframe.md`).
