"""Score predicted models against a reference structure with pdbframe.

    python score.py REFERENCE MODEL_OR_DIR [MODEL_OR_DIR ...] [-o scores.csv]

One row per model: RMSD (CA, after superposition), TM-score, lDDT, pLDDT, the
chain pairing used, and DockQ for every chain-chain interface of the reference.
Chains are paired by sequence, so model chain IDs need not match the reference.
"""

import argparse
import itertools
import sys
from pathlib import Path

import pandas as pd

from pdbframe import PDBFrame

STRUCTURE_SUFFIXES = {'.pdb', '.ent', '.cif', '.mmcif', '.gz'}


def model_files(paths):
    for p in map(Path, paths):
        if p.is_dir():
            yield from sorted(f for f in p.iterdir() if f.suffix.lower() in STRUCTURE_SUFFIXES)
        else:
            yield p


def interfaces(ref, cutoff=5.0):
    """Reference chain pairs with residues within `cutoff` Å of each other."""
    pairs = ref.contacts(cutoff=cutoff)
    return sorted({tuple(sorted(p)) for p in zip(pairs['chain_1'], pairs['chain_2'])})


def score(ref, path, sides):
    model = PDBFrame.read(path)
    row = {'model': path.name}
    try:
        row.update(model.metrics.summary(ref).to_dict())
        for a, b in sides:
            row[f'dockq_{a}-{b}'] = model.metrics.dockq(ref, a, b)['dockq']
    except Exception as e:                      # keep going; report the failure
        row['error'] = f'{type(e).__name__}: {e}'
    return row


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('reference')
    ap.add_argument('models', nargs='+', help='model files or folders of them')
    ap.add_argument('-o', '--out', help='write a CSV here (default: print the table)')
    args = ap.parse_args(argv)

    ref = PDBFrame.read(args.reference)
    sides = interfaces(ref)
    rows = [score(ref, p, sides) for p in model_files(args.models)]
    table = pd.DataFrame(rows)
    if args.out:
        table.to_csv(args.out, index=False)
        print(f'{len(table)} models scored -> {args.out}', file=sys.stderr)
    else:
        print(table.to_string(index=False))


if __name__ == '__main__':
    main()
