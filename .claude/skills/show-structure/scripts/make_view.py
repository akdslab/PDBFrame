"""Save an interactive pdbframe view() page for a structure.

    python make_view.py STRUCTURE -o page.html
        [--region NAME=CHAIN:START-END ...]   residue-number range, inclusive
        [--region NAME=CHAIN ...]             a whole chain
        [--near NAME=CHAIN:CUTOFF ...]        residues within CUTOFF Å of CHAIN
        [--interfaces] [--compare MODEL [--lddt]] [--plddt]
        [--sequence] [--surface] [--sticks NAME ...]

--compare draws MODEL superposed on STRUCTURE with residues linked by sequence;
--lddt then colors STRUCTURE by the model's per-residue lDDT. --plddt colors by
the B-factor column instead. --sticks also draws the named regions as sticks.
"""

import argparse
import re

from pdbframe import PDBFrame


def region_spec(frame, text):
    """'A:50-85' -> ('A', 50, 85); 'A' -> 'A'."""
    m = re.fullmatch(r'([^:]+):(-?\d+)-(-?\d+)', text)
    return (m[1], int(m[2]), int(m[3])) if m else text


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('structure')
    ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--region', action='append', default=[], metavar='NAME=SPEC')
    ap.add_argument('--near', action='append', default=[], metavar='NAME=CHAIN:CUTOFF')
    ap.add_argument('--interfaces', action='store_true')
    ap.add_argument('--compare')
    ap.add_argument('--lddt', action='store_true')
    ap.add_argument('--plddt', action='store_true')
    ap.add_argument('--sequence', action='store_true', help='show every sequence at first')
    ap.add_argument('--surface', action='store_true')
    ap.add_argument('--sticks', nargs='*', default=[], metavar='NAME')
    args = ap.parse_args(argv)

    pf = PDBFrame.read(args.structure, hetero=True)
    regions = {}
    for item in args.region:
        name, spec = item.split('=', 1)
        regions[name] = region_spec(pf, spec)
    for item in args.near:
        name, spec = item.split('=', 1)
        chain, cutoff = spec.rsplit(':', 1)
        target = pf.select(chains=chain)
        regions[name] = pf.near(target, float(cutoff)).select(chains=chain, invert=True)

    kwargs = {}
    if args.compare:
        model = PDBFrame.read(args.compare).superpose(pf)
        kwargs['compare'] = model
        if args.lddt:
            kwargs['values'] = model.metrics.lddt(pf, per_residue=True)
    if args.plddt:
        kwargs['values'] = 'plddt'
    view = pf.view(regions=regions or None, sequence=args.sequence, interfaces=args.interfaces,
                   surface=args.surface, region_styles={n: 'stick' for n in args.sticks}, **kwargs)
    print(view.save(args.out))


if __name__ == '__main__':
    main()
