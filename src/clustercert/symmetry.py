"""Symmetry sectors of a certified cluster: its space group, momenta, stars, little groups and irreps.

Everything is derived from the lattice definition (primitive vectors, sublattice
positions, bonds) and the certificate (periods, site map); nothing is tabulated
per lattice. The steps:

1. **Space group of the infinite lattice.** A symmetry is x -> W x + t in
   fractional coordinates. Its linear part W is an integer matrix (it maps
   lattice translations to lattice translations) that preserves lengths,
   W^T G W = G with G the Gram matrix of the primitive vectors; these are
   finitely many. For each W, the translation part t must move sublattice 0
   onto some sublattice, which leaves n_basis candidates; a candidate is kept
   if it maps every site to a site and every bond to a bond. Glides and screws
   (diamond, pyrochlore) come out of the same search.
2. **Space group of the cluster.** An operation survives the folding if
   W Lambda = Lambda (exact integer test). Combined with the translations of
   the supercell, each one becomes a permutation of your site labels, checked
   to be an automorphism of your bonds.
3. **Momenta, stars and little groups.** The point operations move the allowed
   momenta among themselves; the orbits are the stars. The little group G_k of
   a momentum k is every symmetry whose point operation fixes k up to a
   reciprocal lattice vector.
4. **Irreps.** The irreps of G_k in which translations act as e^{-ik.T} are
   the components of the representation induced from that character of the
   translations. It has dimension |little co-group|, and each irrep of
   dimension d appears in it d times. It is split numerically, by
   diagonalising a random matrix that commutes with it. Non-symmorphic and
   projective cases need no special treatment.
5. **Sector sizes.** The number of times an irrep occurs in the spin-1/2
   Hilbert space follows from its characters and the cycle structure of each
   permutation: a permutation fixes 2^(number of cycles) basis states.

Conventions: a translation by the cell vector c acts on a state of momentum
k as U(T_c) psi = e^{-ik.c} psi, so the one-magnon state sum_x e^{ik.x}|x>
has momentum k. Momenta are fractional coordinates y in the reciprocal
primitive vectors, k.c = 2 pi y.c. The names of point groups (Schoenflies),
irreps (Mulliken) and momenta (Γ, M, K, X, ...) are a cosmetic layer computed
from the operations; see docs/TERMINOLOGY.md.
"""
from dataclasses import dataclass, field
from fractions import Fraction
from itertools import product

import networkx as nx
import numpy as np

from .certify import _family_lattices, _lookup, _reduce, cluster, hnf
from .supercell import momenta

__all__ = ['space_group', 'symmetry', 'Symmetry', 'Sector', 'SpaceGroupOperation']

TOL = 1e-6


# ---------------------------------------------------------------- the infinite lattice

@dataclass(frozen=True)
class SpaceGroupOperation:
    """One symmetry x -> W x + t of the infinite lattice, in fractional coordinates.

    It sends site (c, a) to site (W c + shift[a], sublattice[a]).
    """
    linear: np.ndarray          # W, integer d x d, acting on cell coordinates (columns)
    translation: tuple          # t, fractional (Fractions, approximated), in [0, 1)^d
    sublattice: tuple           # a -> sublattice of the image
    shift: tuple                # a -> integer cell shift of the image
    name: str                   # E, C2, C3, C4, C6, i, σ, S4, S6, S3; with the translation if any

    def __repr__(self):
        return f'SpaceGroupOperation({self.name})'


def _fractional(lattice):
    if lattice.vectors is None:
        raise ValueError(f'{lattice.name} has no geometry (vectors, positions), so its rotations are undefined')
    a = np.array(lattice.vectors, float)
    return a, np.linalg.solve(a.T, np.array(lattice.positions, float).T).T     # rows f with f @ a = position


def _isometries(a):
    """Integer W with W^T G W = G: the linear maps of the lattice of translations that keep lengths (holohedry)."""
    g = a @ a.T
    dim = len(a)
    reach = int(np.ceil(np.sqrt(g.diagonal().max()) / np.linalg.svd(a, compute_uv=False).min() + TOL))
    vecs = np.array(list(product(range(-reach, reach + 1), repeat=dim)))
    norms = np.einsum('ij,jk,ik->i', vecs, g, vecs)
    columns = [vecs[np.abs(norms - g[j, j]) < TOL] for j in range(dim)]
    out = []
    for cols in product(*columns):
        w = np.array(cols).T
        if np.allclose(w.T @ g @ w, g, atol=TOL):
            out.append(w)
    return out


def _cartesian(w, a):
    """The Cartesian matrix R of the linear part W (x = a^T f), embedded in 3D for a 2D lattice."""
    r = a.T @ w @ np.linalg.inv(a.T)
    if len(a) == 2:
        r = np.block([[r, np.zeros((2, 1))], [np.zeros((1, 2)), np.ones((1, 1))]])
    return r


def _classify(r):
    """(kind, order, axis) of a 3x3 orthogonal matrix: kind in E, C, i, σ, S; axis a unit vector or None."""
    det = np.linalg.det(r)
    proper = r if det > 0 else -r
    angle = np.arccos(np.clip((np.trace(proper) - 1) / 2, -1, 1))
    axis = None
    if angle > TOL:
        vals, vecs = np.linalg.eig(proper)
        axis = np.real(vecs[:, np.argmin(np.abs(vals - 1))])
        axis = axis / np.linalg.norm(axis)
        axis = axis * np.sign(axis[np.flatnonzero(np.abs(axis) > TOL)[0]])
    if det > 0:
        return ('E', 1, None) if angle < TOL else ('C', int(round(2 * np.pi / angle)), axis)
    if angle < TOL:
        return 'i', 2, None
    if abs(angle - np.pi) < TOL:
        return 'σ', 2, axis                     # axis = normal of the mirror plane
    theta = (angle + np.pi) % (2 * np.pi)
    return 'S', int(round(2 * np.pi / min(theta, 2 * np.pi - theta))), axis


def _op_name(r):
    kind, order, _ = _classify(r)
    return {'E': 'E', 'i': 'i', 'σ': 'σ'}.get(kind, f'{kind}{order}')


def _frac(x, limit=48):
    return Fraction(float(x)).limit_denominator(limit)


def _fmt_vec(v):
    return '(' + ', '.join(str(x) for x in v) + ')'


def space_group(lattice):
    """Every symmetry of the infinite lattice, up to lattice translations, as :class:`SpaceGroupOperation`.

    Derived from the lattice definition alone: see the module docstring.
    The identity comes first. ``len({op.linear.tobytes() ...})`` is the
    order of the point group (48 for cubic lattices, 12 for hexagonal ones).
    """
    lattice = _lookup(lattice)
    a, f = _fractional(lattice)
    n = lattice.n_basis
    bonds = set(lattice.bonds)
    from .lattice import _canonical
    ops = []
    for w in _isometries(a):
        wf = f @ w.T
        for b in range(n):
            t = f[b] - wf[0]
            t = t - np.floor(t + TOL)
            sub, shift = [], []
            for x in wf + t:
                hits = [(c, np.round(x - f[c])) for c in range(n) if np.allclose(x - f[c], np.round(x - f[c]),
                                                                                 atol=TOL)]
                if len(hits) != 1:
                    break
                sub.append(hits[0][0])
                shift.append(tuple(int(v) for v in hits[0][1]))
            if len(sub) != n or len(set(sub)) != n:
                continue
            s = np.array(shift)
            if all(_canonical(sub[i], sub[j], w @ np.array(tau) + s[j] - s[i]) in bonds for i, j, tau in bonds):
                r = _cartesian(w, a)
                frac_t = tuple(_frac(x) for x in t)
                name = _op_name(r) + ('' if not any(frac_t) else f' + {_fmt_vec(frac_t)}')
                ops.append(SpaceGroupOperation(w.astype(int), frac_t, tuple(sub), tuple(shift), name))
    ops.sort(key=lambda op: (op.name != 'E', _order_key(_cartesian(op.linear, a)), op.translation))
    return ops


def _order_key(r):
    kind, order, _ = _classify(r)
    return 'ECiσS'.index(kind), -order if kind == 'C' else order


def _preserves(op, lattice):
    from .lattice import _canonical
    bonds = set(lattice.bonds)
    s = np.array(op.shift)
    sub = op.sublattice
    return all(_canonical(sub[i], sub[j], op.linear @ np.array(tau) + s[j] - s[i]) in bonds for i, j, tau in bonds)


# ---------------------------------------------------------------- names (cosmetic)

def _axes(info, kind, order=None):
    """Distinct axes (up to sign) of the elements of a kind (and order)."""
    out = []
    for k, n, ax in info:
        if k == kind and (order is None or n == order) and not any(abs(abs(ax @ b) - 1) < TOL for b in out):
            out.append(ax)
    return out


def _principal(info):
    """(order, axis) of the principal axis: highest-order rotation; S4 axis for D2d/S4; nearest z for D2."""
    orders = [n for k, n, _ in info if k == 'C']
    nmax = max(orders, default=1)
    if nmax == 1:
        return 1, None
    s4 = _axes(info, 'S', 4)
    if nmax == 2 and s4:
        return 2, s4[0]
    cands = _axes(info, 'C', nmax)
    z = np.array([0.0, 0, 1])
    return nmax, max(cands, key=lambda ax: (round(abs(ax @ z), 6), round(abs(ax[1]), 6)))


def _schoenflies(rs):
    """Schoenflies name of a crystallographic point group given as Cartesian 3x3 matrices."""
    info = [_classify(r) for r in rs]
    kinds = {k for k, _, _ in info}
    has_i = 'i' in kinds
    if len(_axes(info, 'C', 3)) >= 2:
        if _axes(info, 'C', 4):
            return 'Oh' if has_i else 'O'
        if _axes(info, 'S', 4):
            return 'Td'
        return 'Th' if has_i else 'T'
    n, z = _principal(info)
    mirrors = _axes(info, 'σ')
    if n == 1:
        return 'Cs' if mirrors else 'Ci' if has_i else 'C1'
    c2perp = [ax for ax in _axes(info, 'C', 2) if abs(ax @ z) < TOL]
    sigma_h = any(abs(abs(m @ z) - 1) < TOL for m in mirrors)
    sigma_v = any(abs(m @ z) < TOL for m in mirrors)
    if c2perp:
        return f'D{n}' + ('h' if sigma_h else 'd' if sigma_v else '')
    if sigma_h:
        return f'C{n}h'
    if sigma_v:
        return f'C{n}v'
    if _axes(info, 'S', 2 * n):
        return f'S{2 * n}'
    return f'C{n}'


def _mulliken(rs, chars):
    """Mulliken name of the irrep with characters ``chars`` of the point group ``rs``, or None if not standard.

    A/B: 1-dimensional, symmetric/antisymmetric under the principal rotation;
    E, T: 2- and 3-dimensional; 1/2: symmetric/antisymmetric under the
    reference C2 perpendicular to the principal axis, else the reference
    vertical mirror (the one closest to containing the x axis); g/u: under
    inversion; ' and ": under the horizontal mirror (no inversion).
    """
    info = [_classify(r) for r in rs]
    d = int(round(chars[0].real))
    if abs(chars[0] - d) > TOL or d not in (1, 2, 3):
        return None
    has_i = any(k == 'i' for k, _, _ in info)

    def char(pred):
        vals = [c for (k, n, ax), c in zip(info, chars) if pred(k, n, ax)]
        return vals[0] / d if vals else None

    def sign(x):
        return None if x is None or abs(x.imag) > TOL or abs(abs(x.real) - 1) > TOL else int(round(x.real))

    def suffix():
        out = ''
        if has_i:
            out += 'g' if sign(char(lambda k, n, ax: k == 'i')) == 1 else 'u'
        elif sigma_h is not None:
            out += "'" if sign(sigma_h) == 1 else '"'
        return out

    cubic = len(_axes(info, 'C', 3)) >= 2
    if cubic:
        sigma_h = None
        letter = {1: 'A', 2: 'E', 3: 'T'}[d]
        ref = char(lambda k, n, ax: k == 'C' and n == 4)
        if ref is None:
            ref = char(lambda k, n, ax: k == 'S' and n == 4)
        sub = ''
        if d != 2 and ref is not None:
            s = sign(ref) if d == 1 else (1 if ref.real > 0 else -1)
            sub = '1' if s == 1 else '2'
        if d == 1 and sign(char(lambda k, n, ax: k == 'C' and n == 3)) is None:
            return None                         # complex 1-dimensional irreps of T, Th
        return letter + sub + suffix()
    n, z = _principal(info)
    mirrors = [(c, ax) for (k, _, ax), c in zip(info, chars) if k == 'σ']
    horizontal = [c for c, ax in mirrors if z is not None and abs(abs(ax @ z) - 1) < TOL]
    sigma_h = horizontal[0] / d if horizontal and not has_i else None
    if z is None:                                   # C1, Ci, Cs
        if mirrors and not has_i:
            sigma_h = mirrors[0][0] / d
        return 'A' + suffix() if d == 1 else None
    s4 = n == 2 and bool(_axes(info, 'S', 4))      # D2d, S4: the principal axis is the S4 axis
    principal = char(lambda k, m, ax: (k == 'S' and m == 4 if s4 else k == 'C' and m == n)
                     and abs(abs(ax @ z) - 1) < TOL)
    x = np.array([1.0, 0, 0]) - (z @ [1.0, 0, 0]) * z
    x = x / np.linalg.norm(x) if np.linalg.norm(x) > TOL else np.array([0.0, 1, 0])
    c2perp = [(c, ax) for (k, m, ax), c in zip(info, chars) if k == 'C' and m == 2 and abs(ax @ z) < TOL]
    vertical = [(c, ax) for c, ax in mirrors if abs(ax @ z) < TOL]
    if d == 1:
        p = sign(principal)
        if p is None:
            return None                             # complex pair of a cyclic group
        letter = 'A' if p == 1 else 'B'
        if n == 2 and not s4 and len(c2perp) == 2 and letter == 'B':   # D2, D2h: B1, B2, B3 by the axis that is kept
            y = np.cross(z, x)
            kept = [ax for c, ax in c2perp if sign(c) == 1]
            if not kept:
                return None
            return 'B' + ('2' if abs(kept[0] @ y) > abs(kept[0] @ x) else '3') + suffix()
        if n == 2 and not s4 and len(c2perp) == 2 and letter == 'A':
            if any(sign(c) != 1 for c, _ in c2perp):
                return 'B1' + suffix()
            return 'A' + suffix()
        ref = c2perp or vertical
        if ref:
            key = (lambda item: abs(item[1] @ x)) if c2perp else (lambda item: abs(item[1] @ np.cross(z, x)))
            c = max(ref, key=lambda item: round(key(item), 6))[0]
            return letter + ('1' if sign(c) == 1 else '2') + suffix()
        return letter + suffix()
    if d == 2:
        sub = ''
        if n == 6 and principal is not None:
            sub = '1' if principal.real > 0 else '2'
        return 'E' + sub + suffix()
    return None


def _bravais_points(a, isometries):
    """(a, [(name, |k| a / 2 pi, holohedry order)]): standard high-symmetry momenta of the Bravais lattice type.

    ``a`` is the conventional lattice constant. The type is read off the
    holohedry and the number of shortest translations, so it does not depend
    on the choice of primitive vectors or on the orientation.
    """
    dim = len(a)
    order = len(isometries)
    g = a @ a.T
    lengths = []
    for v in product(range(-2, 3), repeat=dim):
        if any(v):
            lengths.append(float(np.sqrt(np.array(v) @ g @ np.array(v))))
    dmin = min(lengths)
    shortest = sum(abs(x - dmin) < TOL for x in lengths)
    if dim == 2 and order == 12:
        return dmin, [('K', 2 / 3, 6), ('M', 1 / np.sqrt(3), 4)]
    if dim == 2 and order == 8:
        return dmin, [('X', 1 / 2, 4), ('M', 1 / np.sqrt(2), 8)]
    if dim == 3 and order == 48:
        if shortest == 6:       # simple cubic, cube edge = dmin
            return dmin, [('X', 1 / 2, 16), ('M', 1 / np.sqrt(2), 16), ('R', np.sqrt(3) / 2, 48)]
        if shortest == 12:      # fcc, cube edge = sqrt(2) dmin
            return np.sqrt(2) * dmin, [('X', 1, 16), ('L', np.sqrt(3) / 2, 12), ('W', np.sqrt(5) / 2, 8),
                                       ('K', 3 * np.sqrt(2) / 4, 4)]
        if shortest == 8:       # bcc, cube edge = 2 dmin / sqrt 3
            return 2 / np.sqrt(3) * dmin, [('H', 1, 48), ('N', 1 / np.sqrt(2), 8), ('P', np.sqrt(3) / 2, 24)]
    return dmin, []


def _momentum_names(a, ks):
    """Γ, or a standard high-symmetry name, or the fractional coordinates, for each momentum."""
    isometries = _isometries(a)
    scale, points = _bravais_points(a, isometries)
    recip = np.linalg.inv(a).T          # rows b_i / 2 pi
    names = []
    for y in ks:
        yv = np.array([float(v) for v in y])
        if not yv.any():
            names.append('Γ')
            continue
        shortest = min(np.linalg.norm((yv + np.array(m)) @ recip) for m in product(range(-2, 3), repeat=len(a)))
        fixed = sum(_fixes(w, y) for w in isometries)
        name = next((p for p, length, order in points
                     if abs(shortest * scale - length) < 1e-6 and fixed == order), None)
        names.append(name or _fmt_vec(y))
    return names


def _fixes(w, y):
    """Does the point operation W fix momentum y modulo reciprocal lattice vectors? (y -> y W^-1)"""
    image = _act(w, y)
    return image == tuple(y)


def _act(w, y):
    winv = np.round(np.linalg.inv(w)).astype(int)
    image = [sum(y[i] * int(winv[i, j]) for i in range(len(y))) for j in range(len(y))]
    return tuple(x - (x.numerator // x.denominator) for x in image)


# ---------------------------------------------------------------- the cluster

@dataclass
class Sector:
    """One symmetry sector: a momentum (with its star) and an irrep of its little group.

    ``label`` reads ``momentum.irrep``, for example ``Γ.A1`` or ``M.B2``.
    ``dim`` is the dimension d of the irrep of the little group. Every
    eigenvalue of a lattice-symmetric Hamiltonian in this sector is
    ``len(star) * dim``-fold degenerate. ``characters[i]`` is the character
    of ``Symmetry.operations[elements[i]]``.
    """
    label: str
    momentum: tuple
    momentum_name: str
    star: list
    little_cogroup: str
    irrep: str
    dim: int
    elements: list = field(repr=False)
    characters: np.ndarray = field(repr=False)
    star_operations: list = field(repr=False, default=None)     # star[i] = operations[star_operations[i]] k


@dataclass
class Symmetry:
    """Symmetry sectors of a certified cluster. Made by :func:`symmetry`; see :meth:`summary`.

    ``sites``: your site labels, in the order used by ``permutations``.
    ``operations``: names of the symmetries, one per row of ``permutations``;
    row g sends site ``sites[i]`` to ``sites[permutations[g, i]]``.
    ``translations``: indices of the pure translations, with their cell vectors.
    ``point_operations``: one name per kept point operation (coset of the translations).
    """
    lattice: str
    periods: tuple
    sites: list
    operations: list
    permutations: np.ndarray
    translation_cells: dict
    point_operations: list
    point_group: str
    lattice_point_group: str
    lattice_point_group_order: int
    momenta: list
    stars: list
    sectors: list
    n_automorphisms: int = None
    better_supercell: tuple = None      # (periods, point operations) of a matching supercell with more symmetry

    def permutation(self, index):
        """Operation ``index`` as a dict, your label -> your label."""
        return {self.sites[i]: self.sites[j] for i, j in enumerate(self.permutations[index])}

    def generators(self):
        """A few operation indices that generate the whole group."""
        gens, group = [], {self.permutations[0].tobytes()}
        for g, p in enumerate(self.permutations):
            if p.tobytes() in group:
                continue
            gens.append(g)
            group, frontier = {self.permutations[0].tobytes()}, [self.permutations[0]]
            while frontier:
                nxt = []
                for e in frontier:
                    for h in gens:
                        q = self.permutations[h][e]
                        if q.tobytes() not in group:
                            group.add(q.tobytes())
                            nxt.append(q)
                frontier = nxt
        return gens

    def sector(self, label):
        return next(s for s in self.sectors if s.label == label)

    def block_sizes(self, n_up=None, spin_flip=False):
        """Size of the Hamiltonian block of each sector, ``{label: size}``, for spin 1/2.

        ``n_up``: number of up spins (fixed total S^z); None for the whole
        Hilbert space. The size is the number of times the sector's irrep
        occurs: the dimension of the block when you keep one row of each
        irrep, so each of its eigenvalues stands for ``len(star) * dim``
        degenerate states. ``spin_flip=True`` (whole space, or n_up = N/2)
        splits each sector by its parity under flipping every spin, with
        labels ``Γ.A1+`` and ``Γ.A1-``.
        """
        n = len(self.sites)
        if spin_flip and n_up is not None and 2 * n_up != n:
            raise ValueError('flipping every spin is a symmetry only of the whole space or of n_up = N/2')
        cycles = [_cycles(p) for p in self.permutations]
        trace = [_trace(c, n_up) for c in cycles]
        flipped = [(2 ** len(c) if all(x % 2 == 0 for x in c) else 0) for c in cycles] if spin_flip else None
        out = {}
        for s in self.sectors:
            plain = _multiplicity(s, trace)
            if not spin_flip:
                out[s.label] = plain
                continue
            odd = _multiplicity(s, flipped)
            out[s.label + '+'], out[s.label + '-'] = (plain + odd) // 2, (plain - odd) // 2
        return out

    def projector(self, label, n_up=None, momentum=None):
        """Dense projector onto a sector, for small clusters: ``(states, P)``.

        ``states``: the basis, integers whose bit i is the spin (1 = up) of
        ``sites[i]``; all 2^N, or those with ``n_up`` up spins. ``P`` is
        (d/|G_k|) sum_g conj(chi(g)) U(g) over the little group, the
        projector onto the sector at the star's representative momentum, or at
        ``momentum`` (another member of the star) if given. The spin-flip
        parity is not included.
        """
        s = self.sector(label)
        n = len(self.sites)
        states = np.arange(2 ** n) if n_up is None else np.array(
            [x for x in range(2 ** n) if bin(x).count('1') == n_up])
        p = np.zeros((len(states), len(states)), complex)
        columns = np.arange(len(states))
        weight = s.dim / len(s.elements)
        for g, chi in zip(s.elements, s.characters):     # U(g) sends basis state j to state rows[j]
            rows = np.searchsorted(states, _permute_states(states, self.permutations[g]))
            np.add.at(p, (rows, columns), weight * np.conj(chi))
        if momentum is not None and tuple(momentum) != tuple(s.momentum):
            g = s.star_operations[s.star.index(tuple(momentum))]
            rows = np.searchsorted(states, _permute_states(states, self.permutations[g]))
            moved = np.empty_like(p)
            moved[np.ix_(rows, rows)] = p                   # U(g) P U(g)^T
            p = moved
        return states, p

    def summary(self, n_up=None, spin_flip=False):
        """The symmetry analysis as text.

        With ``n_up`` (an integer, or ``'all'`` for the whole Hilbert space),
        also the block size of every sector; see :meth:`block_sizes`.
        """
        n_t, n_p = len(self.translation_cells), len(self.point_operations)
        lines = [f'symmetry of the periodic {self.lattice} cluster, periods {[list(r) for r in self.periods]}',
                 f'  {len(self.sites)} sites, {n_t * n_p} lattice symmetries = {n_t} translations x {n_p} '
                 f'point operations (point group {self.point_group}; the infinite lattice has '
                 f'{self.lattice_point_group}, {self.lattice_point_group_order} operations)']
        if self.n_automorphisms is not None:
            more = self.n_automorphisms - n_t * n_p
            if more == 0:
                lines.append('  every automorphism of the bond graph is a lattice symmetry')
            else:
                lines.append(f'  the bond graph has {self.n_automorphisms} automorphisms, {more} more than the lattice '
                             'symmetries: levels in different sectors can be accidentally degenerate')
        if self.better_supercell:
            other, count = self.better_supercell
            lines.append(f'  the supercell with periods {[list(r) for r in other]} gives the same graph and keeps '
                         f'{count} point operations: pass it to certify(..., periods=...) to use it')
        counts = {}
        for name in self.point_operations:
            counts[name] = counts.get(name, 0) + 1
        lines.append('  point operations: ' + ', '.join(f'{k} x{v}' if v > 1 else k for k, v in counts.items()))
        lines.append(f'  {len(self.momenta)} momenta in {len(self.stars)} stars, k in reciprocal primitive vectors:')
        rows = [('name', 'k', 'star', 'little co-group', 'irreps, (dimension) if not 1')]
        for star in self.stars:
            secs = [s for s in self.sectors if s.momentum == star[0]]
            name = secs[0].momentum_name
            rows.append(('' if name.startswith('(') else name, _fmt_vec(star[0]), str(len(star)),
                         secs[0].little_cogroup, ' '.join(s.irrep + (f'({s.dim})' if s.dim > 1 else '') for s in secs)))
        widths = [max(len(r[c]) for r in rows) for c in range(4)]
        for r in rows:
            lines.append('    ' + '  '.join(x.ljust(w) if c != 2 else x.rjust(w)
                                            for c, (x, w) in enumerate(zip(r[:4], widths))) + '  ' + r[4])
        if n_up is not None:
            sizes = self.block_sizes(None if n_up == 'all' else n_up, spin_flip)
            space = 'all 2^N states' if n_up == 'all' else f'{n_up} up spins'
            lines.append(f'  block sizes ({space}; each level is (star x dimension)-fold degenerate):')
            for star in self.stars:
                labels = {x.label for x in self.sectors if x.momentum == star[0]}
                name = next(x.momentum_name for x in self.sectors if x.momentum == star[0])
                parts = [f'{label.split(".", 1)[1]} {size}' for label, size in sizes.items()
                         if (label[:-1] if spin_flip else label) in labels]
                lines.append(f'    {name}: ' + ', '.join(parts))
        return '\n'.join(lines)


def _cycles(p):
    seen = np.zeros(len(p), bool)
    out = []
    for i in range(len(p)):
        if not seen[i]:
            length, j = 0, i
            while not seen[j]:
                seen[j] = True
                j = p[j]
                length += 1
            out.append(length)
    return out


def _trace(cycles, n_up):
    """Basis states fixed by a permutation with these cycle lengths (all, or with n_up up spins)."""
    if n_up is None:
        return 2 ** len(cycles)
    poly = [1]
    for c in cycles:            # each cycle is all up or all down: multiply by (1 + x^c)
        new = poly + [0] * c
        for i, v in enumerate(poly):
            new[i + c] += v
        poly = new
    return poly[n_up] if n_up < len(poly) else 0


def _multiplicity(sector, trace):
    total = sum(np.conj(chi) * trace[g] for g, chi in zip(sector.elements, sector.characters))
    value = total / len(sector.elements)
    m = int(round(value.real))
    if abs(value - m) > 1e-6 * max(1, abs(value)):
        raise RuntimeError(f'{sector.label}: non-integer multiplicity {value}')
    return m


def _permute_states(states, perm):
    out = np.zeros_like(states)
    for i, j in enumerate(perm):
        out |= ((states >> i) & 1) << int(j)
    return out


def _family_list(result, lattice):
    """[(key, Lattice)], first the connected family from which the space group is derived."""
    if 'families' in result:
        keys = list(result['families'])
        base = lattice if lattice is not None else result['lattice']
        lattices = _family_lattices(dict.fromkeys(keys), _lookup(base))
        return [(k, lattices[k]) for k in keys]
    return [(1, _lookup(lattice if lattice is not None else result['lattice']))]


def symmetry(result, lattice=None, *, automorphisms=True):
    """Symmetry sectors of the cluster certified by ``result``, as a :class:`Symmetry`.

    ``lattice``: needed only for a lattice that is not built in (a
    :class:`Lattice` with geometry), or a dict of family lattices for a
    result with several bond families. ``automorphisms``: also count the
    automorphisms of the bond graph, to flag accidental symmetries.
    """
    if not result.get('certified'):
        raise ValueError(f'only a certified cluster has a symmetry analysis: {result.get("reason")}')
    families = _family_list(result, lattice)
    base = families[0][1]
    a, _ = _fractional(base)
    h = hnf(result['periods'])
    dim = len(h)
    try:
        sites = sorted(result['sites'])
    except TypeError:
        sites = sorted(result['sites'], key=str)
    n = len(sites)
    node_of = [result['sites'][s] for s in sites]
    cells_of = np.array([c for c, _ in node_of], dtype=np.int64)
    sub_of = np.array([b for _, b in node_of])
    radix = [h[j][j] for j in range(dim)]

    def code(cells, subs):
        out = np.zeros(len(subs), dtype=np.int64)
        for j in range(dim):
            out = out * radix[j] + cells[:, j]
        return out * base.n_basis + subs

    lookup = np.full(int(np.prod(radix)) * base.n_basis, -1)
    lookup[code(cells_of, sub_of)] = np.arange(n)

    def reduce_rows(x):
        x = x.copy()
        for j in range(dim):
            q = x[:, j] // h[j][j]
            x -= q[:, None] * np.array(h[j])
        return x

    # kept operations: W Lambda = Lambda, and every bond family preserved
    full = space_group(base)
    lattice_ws = {op.linear.tobytes(): op.linear for op in full}
    kept = [op for op in full
            if all(not any(_reduce(list(op.linear @ np.array(r)), h)) for r in h)
            and all(_preserves(op, lat) for _, lat in families[1:])]
    box = np.array(list(product(*(range(d) for d in radix))), dtype=np.int64).reshape(-1, dim)

    perms, names, affine, seen = [], [], [], {}
    for k, op in enumerate(kept):
        image_cells = cells_of @ op.linear.T + np.array(op.shift)[sub_of]
        image_sub = np.array(op.sublattice)[sub_of]
        for c in box:
            p = lookup[code(reduce_rows(image_cells + c), image_sub)]
            key = p.tobytes()
            if key in seen:
                continue
            seen[key] = len(perms)
            perms.append(p)
            shift = tuple(Fraction(int(x)) for x in c)
            name = op.name if not any(c) else f'{op.name}, then translate by {_fmt_vec(tuple(int(x) for x in c))}'
            names.append(name)
            affine.append((k, tuple(Fraction(t) + s for t, s in zip(op.translation, shift))))
    perms = np.array(perms)

    # every permutation is an automorphism of every bond family, in your labels
    label_index = {s: i for i, s in enumerate(sites)}
    node_label = {node: label for label, node in result['sites'].items()}
    family_edges = []
    for _, lat in families:
        folded = cluster(lat, h)
        edges = {frozenset((label_index[node_label[u]], label_index[node_label[v]])) for u, v in folded.edges}
        family_edges.append(edges)
        for p in perms:
            if {frozenset((int(p[i]), int(p[j]))) for i, j in map(tuple, edges)} != edges:
                raise RuntimeError('a lattice symmetry is not an automorphism of the cluster')

    # translations, cosets of the translations
    translation_cells = {seen[p.tobytes()]: tuple(int(x) for x in c) for c in box
                         for p in [lookup[code(reduce_rows(cells_of + c), sub_of)]]}
    trans_perm = {g: perms[g] for g in translation_cells}
    where = {}                      # element -> (coset, cell c) with element = rep o T_c
    reps = []
    for g in range(len(perms)):
        if perms[g].tobytes() in where:
            continue
        j = len(reps)
        reps.append(g)
        for t, c in translation_cells.items():
            where[perms[g][trans_perm[t]].tobytes()] = (j, c)
    rep_w = [kept[affine[g][0]].linear for g in reps]
    rep_t = [affine[g][1] for g in reps]
    point_names = [_op_name(_cartesian(w, a)) for w in rep_w]
    rs = [_cartesian(w, a) for w in rep_w]
    distinct_w = len({w.tobytes() for w in rep_w}) == len(rep_w)
    point_group = _schoenflies(rs) if distinct_w else f'order {len(reps)}'
    lattice_rs = [_cartesian(w, a) for w in lattice_ws.values()]

    # momenta and stars
    ks = momenta(h)
    knames = dict(zip(ks, _momentum_names(a, ks)))
    remaining, stars = list(ks), []
    while remaining:
        k0 = remaining[0]
        orbit = []
        for w in rep_w:
            image = _act(w, k0)
            if image not in orbit:
                orbit.append(image)
        orbit.sort(key=lambda y: ks.index(y))
        orbit.remove(k0)
        orbit.insert(0, k0)
        stars.append(orbit)
        remaining = [y for y in remaining if y not in orbit]
    stars.sort(key=lambda s: ks.index(s[0]))
    by_name = {}
    for star in stars:
        by_name.setdefault(knames[star[0]], []).append(star[0])
    for name, ys in by_name.items():
        if len(ys) > 1:         # e.g. X points that the cluster's own symmetries do not relate
            for y in ys:
                knames[y] = f'{name}{_fmt_vec(y)}'

    sectors = []
    rng = np.random.default_rng(0)
    for star in stars:
        y = star[0]
        yv = np.array([float(v) for v in y])
        little = [j for j, w in enumerate(rep_w) if _fixes(w, y)]
        star_ops = []
        for target in star:
            j = next(j for j, w in enumerate(rep_w) if _act(w, y) == target)
            star_ops.append(reps[j])

        def lam(c, yv=yv):      # the translation character, e^{-ik.c}
            return np.exp(-2j * np.pi * (yv @ np.array(c, float)))
        m = len(little)
        pos = {j: i for i, j in enumerate(little)}
        d_mats = np.zeros((m, m, m), complex)
        for a_i, j in enumerate(little):
            for b_i, i in enumerate(little):
                dest, c = where[perms[reps[j]][perms[reps[i]]].tobytes()]
                d_mats[a_i, pos[dest], b_i] = lam(c)
        irreps = _decompose(d_mats, rng)
        cog_rs = [rs[j] for j in little]
        cogroup = _schoenflies(cog_rs) if distinct_w else f'order {m}'
        origin = _origin(rep_w, rep_t, little, dim) if distinct_w else None
        named = []
        for chars in irreps:
            name = None
            if distinct_w and (origin is not None or not yv.any()):
                lifted = chars.copy()
                if origin is not None:
                    for a_i, j in enumerate(little):
                        shift = [origin[r] - sum(int(rep_w[j][r, s]) * origin[s] for s in range(dim)) - rep_t[j][r]
                                 for r in range(dim)]
                        lifted[a_i] = chars[a_i] * lam(shift)
                name = _mulliken(cog_rs, lifted)
            named.append(name)
        # A1g A2g ... Eg ... T2g, then the u (or ") irreps; unnamed ones by dimension
        order = sorted(range(len(irreps)), key=lambda i, named=named, irreps=irreps: (
            named[i] is None, (named[i] or '').endswith(('u', '"')), int(round(irreps[i][0].real)), named[i] or ''))
        final = [named[i] or str(index) for index, i in enumerate(order, 1)]
        repeats = {x: 0 for x in final}
        for index, i in enumerate(order):
            chars = irreps[i]
            name = final[index]
            if final.count(name) > 1:
                repeats[name] += 1
                name = f'{name}({repeats[name]})'
            elements, values = [], []
            for a_i, j in enumerate(little):
                for t, c in translation_cells.items():
                    elements.append(seen[perms[reps[j]][trans_perm[t]].tobytes()])
                    values.append(chars[a_i] * lam(c))
            label = f'{knames[y]}.{name}'
            sectors.append(Sector(label, y, knames[y], star, cogroup, name, int(round(chars[0].real)),
                                  elements, np.array(values), star_ops))

    n_aut = _count_automorphisms(family_edges, perms) if automorphisms else None
    better, most = None, len(reps)
    for other in result.get('matching', []):        # another supercell for the same graph, with more symmetry?
        count = len({op.linear.tobytes() for op in full
                     if all(not any(_reduce(list(op.linear @ np.array(r)), other)) for r in other)
                     and all(_preserves(op, lat) for _, lat in families[1:])})
        if count > most:
            better, most = (other, count), count
    return Symmetry(base.name, h, sites, names, perms, translation_cells, point_names, point_group,
                    _schoenflies(lattice_rs), len(lattice_rs), ks, stars, sectors, n_aut, better)


def _decompose(d_mats, rng, attempts=5):
    """Characters of the irreps in the representation D (matrices over coset reps), each irrep once."""
    m = d_mats.shape[1]
    for _ in range(attempts):
        x = rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m))
        x = x + x.conj().T
        c = sum(d @ x @ d.conj().T for d in d_mats)
        vals, vecs = np.linalg.eigh(c)
        groups, start = [], 0
        scale = max(1.0, np.abs(vals).max())
        for i in range(1, m + 1):
            if i == m or vals[i] - vals[i - 1] > 1e-7 * scale:
                groups.append(vecs[:, start:i])
                start = i
        found = []
        for v in groups:
            chars = np.einsum('ia,gij,ja->g', v.conj(), d_mats, v)
            if not any(np.allclose(chars, f, atol=1e-6) for f in found):
                found.append(chars)
        dims = [int(round(f[0].real)) for f in found]
        norms = [np.vdot(f, f).real / m for f in found]
        if sum(x * x for x in dims) == m and np.allclose(norms, 1, atol=1e-6):
            return found
    raise RuntimeError('could not split the little-group representation into irreps')


def _origin(ws, ts, little, dim, denominator=24):
    """A point fixed (mod lattice) by every little co-group operation, the most symmetric one; None if none."""
    grid = np.array(list(product(range(denominator), repeat=dim))) / denominator

    def fixed_by(j):
        w = ws[j].astype(float)
        t = np.array([float(x) for x in ts[j]])
        r = grid @ w.T + t - grid
        return np.all(np.abs(r - np.round(r)) < 1e-9, axis=1)

    ok = np.all([fixed_by(j) for j in little], axis=0)
    if not ok.any():
        return None
    score = np.sum([fixed_by(j) for j in range(len(ws))], axis=0)
    best = np.flatnonzero(ok)[np.argmax(score[ok])]
    return [Fraction(int(round(x * denominator)), denominator) for x in grid[best]]


def _count_automorphisms(family_edges, perms):
    """Order of the automorphism group of the bond graph, bond families kept apart.

    Orbit-stabilizer along a chain: pin sites one at a time; |Aut| is the
    product of the orbit sizes of each new site under the automorphisms that
    fix the earlier ones. Each orbit contains the orbit under the lattice
    symmetries that fix the earlier sites; every other candidate is tested
    with one isomorphism search. The chain stops when colour refinement
    separates every site, so that only the identity fixes the pinned sites.
    """
    n = perms.shape[1]
    neighbours = [[] for _ in range(n)]
    graph = nx.Graph()
    graph.add_nodes_from(range(n))
    for f, edges in enumerate(family_edges):
        for e in edges:
            u, v = tuple(e)
            neighbours[u].append((f, v))
            neighbours[v].append((f, u))
            graph.add_edge(u, v, kind=f)

    def refine(pins):
        colour = [0] * n
        for c, site in enumerate(pins, 1):
            colour[site] = c
        classes = len(set(colour))
        while True:
            signature = [(colour[i], tuple(sorted((f, colour[j]) for f, j in neighbours[i]))) for i in range(n)]
            names = {sig: c for c, sig in enumerate(sorted(set(signature)))}
            colour = [names[sig] for sig in signature]
            if len(names) == classes:
                return colour
            classes = len(names)

    def exists(pins, images):
        g1, g2 = graph.copy(), graph.copy()
        nx.set_node_attributes(g1, 0, 'colour')
        nx.set_node_attributes(g2, 0, 'colour')
        for c, (a, b) in enumerate(zip(pins, images), 1):
            g1.nodes[a]['colour'] = g2.nodes[b]['colour'] = c
        if len(family_edges) == 1:
            return nx.vf2pp_is_isomorphic(g1, g2, node_label='colour')
        return nx.algorithms.isomorphism.GraphMatcher(
            g1, g2, node_match=lambda a, b: a['colour'] == b['colour'],
            edge_match=lambda a, b: a['kind'] == b['kind']).is_isomorphic()

    order, pins, fixing = 1, [], perms
    while True:
        colour = refine(pins)
        cells = {}
        for site, c in enumerate(colour):
            cells.setdefault(c, []).append(site)
        big = [cell for cell in cells.values() if len(cell) > 1]
        if not big:
            return order
        site = min(big, key=len)[0]
        orbit = set(fixing[:, site].tolist())
        untested = set(cells[colour[site]]) - orbit
        while untested:
            other = min(untested)
            found = set(fixing[:, other].tolist())
            untested -= found
            if exists(pins + [site], pins + [other]):
                orbit |= found
        order *= len(orbit)
        pins.append(site)
        fixing = fixing[fixing[:, site] == site]
