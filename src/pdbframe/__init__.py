"""pdbframe - pandas-based protein sequence and structure objects.

    pdbframe.sequence   SeqFrame: a sequence, or an alignment, as a per-position DataFrame
    pdbframe.structure  PDBFrame: a structure as an atom DataFrame, with comparison metrics
    pdbframe.resources  the docs and Claude Code skills: pdbframe.docs(), install_skills()
"""

from .sequence import SeqFrame, detect_type
from .structure import PDBFrame
from .structure.viewer import view_many as view
from .resources import docs, install_skills

__all__ = ['SeqFrame', 'PDBFrame', 'detect_type', 'view', 'docs', 'install_skills']
__version__ = '0.2.0'
