# Changelog

All notable changes to this project are listed here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). Minor edits bump the last
number. Changes that would break existing code are announced under *Deprecated* first.

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
