# Changelog

All notable changes to this project are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Minor edits bump the last
number. Changes that would break existing code are announced under *Deprecated* first.

## [0.2.1] - 2026-11-18

### Added
- `view(spheres=...)`: balls drawn in the structure's coordinates, with `center`,
  `radius`, `color`, `label`, `group` and `info`. They share the arrows' group
  checkboxes, so the edges and nodes of a graph can switch together.
- `view(opacity=...)`: how solid the structure is drawn, 0 to 1. Below 1 it keeps its
  shape and colors and lets drawn shapes show through.
- `view(sequence_only=True)`: the sequences only, with their regions and values; a *Show 3D*
  button opens the 3D view under them.
- Arrows take `head=False`, drawing a plain cylinder for something with no direction,
  and `info=`, text shown in the view when that arrow or sphere is clicked.

## [0.2.0] - 2026-11-13

First public release.

- `PDBFrame`: a structure as a table of atoms. Reads PDB, mmCIF, CSV and FASTA (gzipped
  too) or fetches from the RCSB PDB; selections, contacts, interfaces, renaming,
  renumbering and joining; writes PDB, mmCIF, CSV and FASTA.
- Unobserved residues: from an mmCIF file, `pf.residues` has every residue of the
  deposited sequence (`observed`, `observed_index`, `full_index`), and mmCIF output keeps
  them.
- Scores against a reference: RMSD, TM-score, lDDT, pLDDT, DockQ per interface.
- `SeqFrame`: a sequence or an alignment as a per-position table, with identity,
  similarity and gap statistics.
- Interactive 3D views linked to the sequence, with regions, interfaces and compared
  models.
- The docs and four Claude Code skills come with the package: `pdbframe docs`,
  `pdbframe skills install`.
