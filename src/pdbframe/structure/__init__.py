"""Structure tools: the PDBFrame object.

    pdbframe.py   the PDBFrame class
    io.py         reading / writing PDB, mmCIF and CSV as atom tables
    deposited.py  the deposited sequence: unobserved residues beside the atoms
    rcsb.py       downloading entries and their metadata from the RCSB PDB
    selection.py  picking atoms and residues, and finding contacts
    matching.py   pairing the chains, residues and atoms of two structures
    geometry.py   superposition and distances on coordinate arrays
    metrics.py    all scores: pf.metrics (RMSD, TM-score, lDDT, DockQ, pLDDT) and the maths
    viewer/       pf.view(): 3D view with each chain's sequence linked to it
    display.py    pf.py3Dmol(): py3Dmol view of one or more structures
"""

from .pdbframe import ATOM_SETS, PDBFrame

__all__ = ['PDBFrame', 'ATOM_SETS']
