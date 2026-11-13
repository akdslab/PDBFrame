"""The structure viewer behind ``pf.view()``: 3D, with each chain's sequence linked to it.

    widget.py    Python side: turns a PDBFrame and options into the widget's data
    viewer.js    browser side: sequence panel, 3D view and their interactions
    viewer.css   browser side: layout

What it does:

- Hover a residue in the sequence and it lights up in 3D, and the reverse.
- Drag across the sequence to select a range; click a chain label to select the chain.
- Ctrl/Cmd-click residues or chain labels (in either panel) to pick them: their details
  and the distances between them are listed, and drawn in 3D.
- Marked regions are shaded in their own color in both.
- The sequence panel can be hidden, and a chain menu picks which sequences it shows.
- Arrows in the corner show the x, y, z axes as the view turns.

Plain HTML + JavaScript on 3Dmol.js (loaded from a CDN, as py3Dmol does), so it
needs no extra Python packages. It renders in Jupyter / VS Code notebooks, and
``.save('view.html')`` writes a standalone page.
"""

from .widget import (COMPARE_COLOR, HIGHLIGHT_COLOR, PICK_COLOR, REGION_COLORS, STRUCTURE_COLORS,
                     Viewer, view, view_many)

__all__ = ['view', 'view_many', 'Viewer', 'REGION_COLORS', 'HIGHLIGHT_COLOR', 'PICK_COLOR',
           'COMPARE_COLOR', 'STRUCTURE_COLORS']
