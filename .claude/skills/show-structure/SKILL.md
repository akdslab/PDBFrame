---
name: show-structure
description: Make an interactive 3D view of a protein structure with its sequence linked (pdbframe view()), marking regions, interfaces, a compared model or per-residue scores, and save it as a shareable HTML page. Use when asked to show, visualize or highlight parts of a structure.
---

# Show a structure

Pick what answers the question, then save a page (it opens in any browser; needs internet
for 3Dmol.js). In a notebook the returned object displays itself.

**From the command line** (bundled script):
```bash
python scripts/make_view.py ref.cif -o page.html \
    --region helix=A:50-85 --region peptide=C --near pocket=C:4 --sticks pocket --sequence
python scripts/make_view.py ref.cif -o page.html --compare model.pdb --lddt --interfaces
```

**In Python**
```python
pf.view(regions={'pocket': pf.near(pf.select(chains='C'), 4), 'helix': ('A', 50, 85),
                 'peptide': 'C'},            # spec: chain | (chain, start, end) | PDBFrame | list
        region_styles={'pocket': 'stick'},   # also draw a region as sticks
        sequence=True,                       # or a list of chains; none shown by default
        interfaces=True,                     # or [('A', 'C')]
        compare=model.superpose(pf),         # residues linked by sequence; grey
        values=model.metrics.lddt(pf, per_residue=True),   # or 'plddt'; red low -> blue high
        surface=False, style='cartoon').save('page.html')
pf.py3Dmol(other, colors=['lightgrey', 'orange'])   # plain overlay of several structures
```

**Several structures and switchable regions**
```python
import pdbframe
pdbframe.view([ref, m1, m2], names=['crystal', 'AF3', 'Boltz'],   # first leads; all paired with it
              style=['cartoon', 'cartoon', 'line'],                # per-structure options: lists
              superpose=['A'],                                     # fit on lead chain A ('all': all)
              values=[m1.metrics.lddt(ref, per_residue=True), None, None])
# a values table must be keyed by that structure's residues; lDDT tables use the reference's
pf.view(regions={...}, region_groups={'CDRs': ['CDR1', 'CDR2']}, hidden_regions=['CDRs'],
        orientations={'standard': R}, orientation_center='structure',  # or 'frame'
        orientation='standard',              # open in that orientation (or a 3x3 rotation)
        arrows=[(start, end, 'red', 'axis'),  # vectors in structure coordinates, or dicts with
                {'start': s, 'end': e, 'group': 'normals'}])  # color, label, radius, group
```
Legend chips switch regions (Shift-click: only this one). *Structures* / *One at a time*
(◀ ▶, arrow keys) flip between models. `compact=True` keeps big ensembles small.

Tell the user how to use the page: hover links sequence and 3D; drag across the sequence
selects a range; Ctrl/⌘-click picks residues or chains and lists distances; Ctrl+drag pans,
scroll zooms; the **?** button lists everything. Say which color is which.

Every method and option: `pdbframe docs pdbframe` (prints `docs/pdbframe.md`).
