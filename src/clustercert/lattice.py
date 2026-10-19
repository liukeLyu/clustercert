"""Periodic lattices, each described by the contents of one unit cell.

Every periodic lattice repeats one unit cell. A :class:`Lattice` records

- ``dim``: the dimension d;
- ``n_basis``: the number n of sites in each unit cell, called sublattices
  and numbered 0..n-1;
- ``bonds``: one entry ``(a, b, t)`` per bond of a unit cell. It means that
  sublattice a in unit cell c is bonded to sublattice b in unit cell c + t,
  for every c. Here c and t are integer vectors, counted in primitive vectors;
- ``vectors`` and ``positions`` (optional geometry): the primitive vectors
  a_1..a_d (one per row) and the Cartesian position of each sublattice in
  unit cell 0.

A site of the infinite lattice is a pair (c, a), at position
c_1 a_1 + ... + c_d a_d + positions[a]. The certificate uses only the
integers (n_basis, bonds); the geometry is used to report period vectors,
positions and momenta, and to build further-neighbour shells.
"""
from itertools import product
from math import ceil, sqrt
from typing import NamedTuple

import numpy as np

__all__ = ['Lattice', 'SQUARE', 'TRIANGULAR', 'HONEYCOMB', 'KAGOME', 'CHECKERBOARD', 'CUBIC',
           'BCC', 'FCC', 'DIAMOND', 'PYROCHLORE', 'LATTICES', 'line_graph_lattice', 'lattice_from_geometry',
           'neighbour_distances', 'neighbour_shell', 'bond_families']


class Lattice(NamedTuple):
    """An infinite periodic lattice, given by one unit cell (see the module docstring)."""
    name: str
    dim: int
    n_basis: int                # sublattices (sites) per unit cell
    bonds: tuple                # (a, b, t): sublattice a in unit cell c -- sublattice b in unit cell c + t
    vectors: tuple = None       # primitive vectors, one per row (Cartesian)
    positions: tuple = None     # Cartesian position of each sublattice in unit cell 0

    def describe(self):
        """The lattice in words: primitive vectors, sublattices and bonds."""
        num = lambda x: f'{x:.4g}' if abs(x) > 1e-12 else '0'  # noqa: E731
        vec = lambda v: '(' + ', '.join(num(x) for x in v) + ')'  # noqa: E731
        lines = [f'{self.name}: {self.dim}D, {self.n_basis} sublattice(s) per unit cell, '
                 f'{len(self.bonds)} bonds per unit cell']
        if self.vectors is not None:
            lines += [f'  primitive vector a{i+1} = {vec(v)}' for i, v in enumerate(self.vectors)]
            lines += [f'  sublattice {a} at {vec(p)}' for a, p in enumerate(self.positions)]
        for a, b, t in self.bonds:
            line = (f'  bond ({a}, {b}, {tuple(t)}): sublattice {a} in unit cell c -- '
                    f'sublattice {b} in unit cell c + {tuple(t)}')
            if self.vectors is not None:
                r = np.array(self.positions[b]) + np.array(t) @ np.array(self.vectors) - np.array(self.positions[a])
                line += f', length {num(float(np.linalg.norm(r)))}'
            lines.append(line)
        return '\n'.join(lines)

    def degrees(self):
        """Coordination number of each basis site in the infinite lattice."""
        deg = [0]*self.n_basis
        for a, b, _ in self.bonds:
            deg[a] += 1
            deg[b] += 1
        return deg


def _canonical(a, b, t):
    """One representative of the bond (a, b, t) = (b, a, -t)."""
    t = tuple(int(x) for x in t)
    reverse = (b, a, tuple(-x for x in t))
    return min((a, b, t), reverse)


def _make(name, bonds, vectors, positions):
    canonical = sorted({_canonical(*bond) for bond in bonds})
    if len(canonical) != len(bonds):
        raise ValueError(f'{name}: a bond is listed twice')
    if any(a == b and not any(t) for a, b, t in canonical):
        raise ValueError(f'{name}: a site is bonded to itself')
    return Lattice(name, len(vectors), len(positions), tuple(canonical),
                 tuple(map(tuple, vectors)), tuple(map(tuple, positions)))


def _e(dim, i):
    return tuple(int(k == i) for k in range(dim))


def _corner_sharing(name, vectors):
    """Kagome (dim 2) or pyrochlore (dim 3): corner-sharing triangles or tetrahedra.

    Basis site 0 at the cell origin and site i at a_i/2. The up simplex of cell
    c is (c, 0..dim); the down simplex is (c, 0) with (c - e_i, i).
    """
    dim = len(vectors)
    zero = (0,)*dim
    up = [(a, b, zero) for a in range(dim+1) for b in range(a+1, dim+1)]
    down = [(0, i, tuple(-x for x in _e(dim, i-1))) for i in range(1, dim+1)]
    down += [(i, j, tuple(x-y for x, y in zip(_e(dim, i-1), _e(dim, j-1))))
             for i in range(1, dim+1) for j in range(i+1, dim+1)]
    positions = [(0.0,)*dim] + [tuple(x/2 for x in v) for v in vectors]
    return _make(name, up+down, vectors, positions)


def line_graph_lattice(lattice, name=None):
    """Lattice of the line graph: one site per bond, bonded when two bonds share a site.

    Kagome is the line graph of honeycomb, pyrochlore of diamond, and
    checkerboard of square. Site k sits at the midpoint of bond k.
    """
    bonds = set()
    for (k, (ak, bk, tk)), (m, (am, bm, tm)) in product(enumerate(lattice.bonds), repeat=2):
        tk, tm = np.array(tk), np.array(tm)
        # bond k at cell 0 has ends (0, ak), (tk, bk); bond m at cell s has ends (s, am), (s + tm, bm)
        for end_k, site_k in ((0*tk, ak), (tk, bk)):
            for end_m, site_m in ((0*tm, am), (tm, bm)):
                if site_k == site_m:
                    s = tuple(int(x) for x in end_k - end_m)
                    if k != m or any(s):
                        bonds.add(_canonical(k, m, s))
    positions = None
    if lattice.vectors is not None:
        vec = np.array(lattice.vectors)
        pos = np.array(lattice.positions)
        positions = [(pos[a] + pos[b] + np.array(t) @ vec)/2 for a, b, t in lattice.bonds]
        # move every site into unit cell 0: shifting site k by s_k cells changes bond (a, b, t) to t + s_b - s_a
        shift = [np.floor(np.linalg.solve(vec.T, p) + 1e-9).astype(int) for p in positions]
        positions = [p - s @ vec for p, s in zip(positions, shift)]
        bonds = {_canonical(a, b, tuple(np.array(t) + shift[b] - shift[a])) for a, b, t in bonds}
    return _make(name or f'line graph of {lattice.name}', sorted(bonds), lattice.vectors, positions)


def lattice_from_geometry(name, vectors, positions, bond_length, tol=1e-8):
    """Bond every pair of sites at Cartesian distance ``bond_length``.

    ``vectors``: primitive vectors (rows). ``positions``: basis sites in
    Cartesian coordinates. Every periodic image is searched, so the result is
    complete. For lattices whose bonds are not fixed by one distance (for
    example checkerboard), list the bonds yourself in a :class:`Lattice`.
    """
    vec = np.array(vectors, float)
    pos = np.array(positions, float)
    # |t @ vec| <= bond_length + 2 max|pos|, and |t|_inf <= |t @ vec| / smallest singular value
    reach = ceil((bond_length + 2*np.abs(pos).sum(1).max()) / np.linalg.svd(vec, compute_uv=False).min())
    bonds = set()
    for t in product(range(-reach, reach+1), repeat=len(vec)):
        shift = np.array(t) @ vec
        for a, b in product(range(len(pos)), repeat=2):
            if abs(np.linalg.norm(pos[b] + shift - pos[a]) - bond_length) < tol:
                bonds.add(_canonical(a, b, t))
    if not bonds:
        raise ValueError(f'{name}: no pair of sites at distance {bond_length}')
    return _make(name, sorted(bonds), vec.tolist(), pos.tolist())


def neighbour_distances(lattice, count, tol=1e-6):
    """The ``count`` smallest distances between distinct sites of the lattice."""
    vec, pos = np.array(lattice.vectors, float), np.array(lattice.positions, float)
    sigma = np.linalg.svd(vec, compute_uv=False).min()
    spread = 2*np.abs(pos).sum(1).max()
    reach = 1
    while True:
        # every pair whose cells differ by more than ``reach`` is at least ``guaranteed`` apart
        guaranteed = sigma*(reach+1) - spread
        found = []
        for t in product(range(-reach, reach+1), repeat=len(vec)):
            r = pos[None, :, :] + np.array(t) @ vec - pos[:, None, :]
            found.extend(np.linalg.norm(r, axis=-1).ravel().tolist())
        shells = []
        for d in sorted(found):
            if d > tol and (not shells or d - shells[-1] > tol):
                shells.append(d)
        if len(shells) >= count and shells[count-1] < guaranteed - tol:
            return shells[:count]
        reach += 1


def neighbour_shell(lattice, k, tol=1e-6):
    """Lattice whose bonds join every pair of sites at the k-th neighbour distance.

    Shell 1 is the nearest neighbours. A shell can contain inequivalent bonds
    at the same distance (on kagome, third neighbours across a hexagon and
    along a line); :func:`bond_families` separates them.
    """
    if lattice.vectors is None:
        raise ValueError(f'{lattice.name} has no geometry, so its neighbour shells are undefined')
    distance = neighbour_distances(lattice, k, tol)[-1]
    shell = lattice_from_geometry(f'{lattice.name} shell {k}', lattice.vectors, lattice.positions, distance, tol)
    return shell._replace(name=f'{lattice.name} shell {k} (distance {distance:.4g})')


def _steps(lattice, a, b, t):
    """J1 steps from sublattice a in cell 0 to sublattice b in cell t, and the number of shortest J1 paths."""
    neighbours = {}
    for x, y, u in lattice.bonds:
        neighbours.setdefault(x, []).append((y, np.array(u)))
        neighbours.setdefault(y, []).append((x, -np.array(u)))
    start, goal = ((0,)*lattice.dim, a), (tuple(t), b)
    paths, layer, seen, depth = {start: 1}, {start: 1}, {start}, 0
    while goal not in paths:
        depth += 1
        nxt = {}
        for (c, x), count in layer.items():
            for y, u in neighbours.get(x, ()):
                site = (tuple(np.add(c, u).tolist()), y)
                if site not in seen:
                    nxt[site] = nxt.get(site, 0) + count
        if not nxt:
            return None, 0
        seen.update(nxt)
        paths.update(nxt)
        layer = nxt
    return depth, paths[goal]


def bond_families(lattice, shells=3, tol=1e-6):
    """The bond families J1, J2, ... up to the given neighbour shell, each as a Lattice.

    A distance shell can hold inequivalent bonds. Kagome's third shell, for
    example, joins sites along a straight line (2 J1 steps) and across a
    hexagon (3 J1 steps). Such a shell is split by two numbers that every
    lattice symmetry preserves: how many J1 steps separate the two sites, and
    along how many shortest J1 paths. Equivalent bonds always agree on both,
    so a bond type is never broken up. Keys are 1, 2, 3 for an unsplit shell
    and '3a', '3b', ... (fewest J1 steps first) for a split one; each name
    says what distinguishes it.
    """
    families = {1: lattice}
    for k in range(2, shells + 1):
        shell = neighbour_shell(lattice, k, tol)
        distance = neighbour_distances(lattice, k, tol)[-1]
        groups = {}
        for bond in shell.bonds:
            groups.setdefault(_steps(lattice, *bond), []).append(bond)
        if len(groups) == 1:
            families[k] = shell
            continue
        for letter, key in zip('abcdefgh', sorted(groups, key=lambda g: (g[0] is None, g[0] or 0, -g[1]))):
            steps, paths = key
            name = (f'{lattice.name} J{k}{letter} (distance {distance:.4g}, {steps} J1 steps, '
                    f'{paths} shortest path{"s" if paths != 1 else ""})')
            families[f'{k}{letter}'] = _make(name, groups[key], lattice.vectors, lattice.positions)
    return families


_HALF_SQRT3 = sqrt(3)/2
_TRI = ((1.0, 0.0), (0.5, _HALF_SQRT3))
_FCC = ((0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0))

SQUARE = _make('square', [(0, 0, (1, 0)), (0, 0, (0, 1))], ((1.0, 0.0), (0.0, 1.0)), [(0.0, 0.0)])
TRIANGULAR = _make('triangular', [(0, 0, (1, 0)), (0, 0, (0, 1)), (0, 0, (-1, 1))], _TRI, [(0.0, 0.0)])
HONEYCOMB = _make('honeycomb', [(0, 1, (0, 0)), (0, 1, (-1, 0)), (0, 1, (0, -1))],
                  _TRI, [(0.0, 0.0), (0.5, _HALF_SQRT3/3)])
KAGOME = _corner_sharing('kagome', _TRI)
CHECKERBOARD = line_graph_lattice(SQUARE, 'checkerboard')
CUBIC = _make('cubic', [(0, 0, (1, 0, 0)), (0, 0, (0, 1, 0)), (0, 0, (0, 0, 1))],
              ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)), [(0.0, 0.0, 0.0)])
BCC = _make('bcc', [(0, 0, (1, 0, 0)), (0, 0, (0, 1, 0)), (0, 0, (0, 0, 1)), (0, 0, (1, 1, 1))],
            ((-0.5, 0.5, 0.5), (0.5, -0.5, 0.5), (0.5, 0.5, -0.5)), [(0.0, 0.0, 0.0)])
FCC = _make('fcc', [(0, 0, (1, 0, 0)), (0, 0, (0, 1, 0)), (0, 0, (0, 0, 1)),
                    (0, 0, (1, -1, 0)), (0, 0, (1, 0, -1)), (0, 0, (0, 1, -1))], _FCC, [(0.0, 0.0, 0.0)])
DIAMOND = _make('diamond', [(0, 1, (0, 0, 0)), (0, 1, (1, 0, 0)), (0, 1, (0, 1, 0)), (0, 1, (0, 0, 1))],
                _FCC, [(0.0, 0.0, 0.0), (-0.25, -0.25, -0.25)])
PYROCHLORE = _corner_sharing('pyrochlore', _FCC)

LATTICES = {m.name: m for m in (SQUARE, TRIANGULAR, HONEYCOMB, KAGOME, CHECKERBOARD,
                                CUBIC, BCC, FCC, DIAMOND, PYROCHLORE)}
