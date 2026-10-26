"""Command line: ``clustercert EDGES [LATTICE] [--shell K FILE ...] [--all] [--symmetry [--n-up N]] [--json OUT]``."""
import argparse
import json
import sys

import numpy as np

from .certify import certify
from .lattice import LATTICES
from .report import summary
from .symmetry import symmetry


def read_edges(path):
    """One bond per line, two site labels separated by spaces or a comma. ``#`` starts a comment."""
    edges = []
    with open(path) as f:
        for number, line in enumerate(f, 1):
            fields = line.split('#')[0].replace(',', ' ').split()
            if not fields:
                continue
            if len(fields) < 2:
                raise SystemExit(f'{path}:{number}: expected two site labels, got {line.strip()!r}')
            edges.append(tuple(int(x) if x.lstrip('-').isdigit() else x for x in fields[:2]))
    return edges


def main(argv=None):
    parser = argparse.ArgumentParser(prog='clustercert', description=__doc__.split(':', 1)[1].strip())
    parser.add_argument('edges', help='edge-list file')
    parser.add_argument('lattice', nargs='?', choices=sorted(LATTICES),
                        help='lattice to test (default: try every built-in lattice)')
    parser.add_argument('--shell', nargs=2, action='append', metavar=('K', 'FILE'), default=[],
                        help='further-neighbour bonds: FILE lists the bonds of neighbour shell K (2, 3, ...)')
    parser.add_argument('--all', action='store_true', help='list every matching supercell')
    parser.add_argument('--symmetry', action='store_true',
                        help='also list the symmetry sectors: momenta, little groups and irreps')
    parser.add_argument('--n-up', type=int, metavar='N',
                        help='with --symmetry: block size of every sector with N up spins (spin 1/2)')
    parser.add_argument('--json', metavar='OUT', help='write the full result, including the site map, to OUT')
    args = parser.parse_args(argv)
    edges = read_edges(args.edges)
    if args.shell:
        if not args.lattice:
            parser.error('--shell needs a lattice')
        edges = {1: edges, **{int(k): read_edges(f) for k, f in args.shell}}
    lattices = [LATTICES[args.lattice]] if args.lattice else list(LATTICES.values())
    results = []
    for lattice in lattices:
        result = certify(edges, lattice, find_all=args.all)
        if args.lattice or result['certified']:
            print(summary(result, lattice))
        if args.symmetry and result['certified']:
            print(symmetry(result, lattice).summary(n_up=args.n_up))
        results.append(result)
    if not args.lattice and not any(r['certified'] for r in results):
        print('not a periodic cluster of any built-in lattice:')
        for r in results:
            print(f'  {r["lattice"]}: {r["reason"]}')
    if args.json:
        with open(args.json, 'w') as f:
            json.dump([_jsonable(r) for r in results], f, indent=1)
    return 0 if any(r['certified'] for r in results) else 1


def _jsonable(result):
    out = {k: v for k, v in result.items() if k != 'sites'}
    if 'sites' in result:
        out['sites'] = [dict(site=s, unit_cell=list(c), sublattice=a) for s, (c, a) in result['sites'].items()]
    return json.loads(json.dumps(out, default=lambda x: x.tolist() if isinstance(x, np.ndarray) else str(x)))


if __name__ == '__main__':
    sys.exit(main())
