"""Geometry of a certified cluster: period vectors, site positions, momenta."""
from fractions import Fraction
from itertools import product

import numpy as np

from .certify import hnf

__all__ = ['period_vectors', 'site_positions', 'momenta']


def _reduced_basis(h, gram):
    """LLL-reduce the integer rows h in the metric ``gram`` (a short, nearly orthogonal basis)."""
    b = [np.array(row, dtype=np.int64) for row in h]
    dot = lambda u, v: float(u @ gram @ v)  # noqa: E731
    k = 1
    while k < len(b):
        for j in range(k-1, -1, -1):
            star = _gram_schmidt(b, dot)
            mu = dot(b[k], star[j]) / dot(star[j], star[j])
            b[k] = b[k] - round(mu)*b[j]
        star = _gram_schmidt(b, dot)
        mu = dot(b[k], star[k-1]) / dot(star[k-1], star[k-1])
        if dot(star[k], star[k]) >= (0.99 - mu**2) * dot(star[k-1], star[k-1]):
            k += 1
        else:
            b[k-1], b[k] = b[k], b[k-1]
            k = max(k-1, 1)
    return b


def _gram_schmidt(b, dot):
    star = []
    for v in b:
        w = v.astype(float)
        for s in star:
            w = w - dot(v, s)/dot(s, s)*s
        star.append(w)
    return star


def period_vectors(periods, lattice):
    """A short basis of the translations ``periods`` (integer rows, primitive basis).

    Returns ``(integer, cartesian)``: rows in the primitive basis, and the
    same rows in Cartesian coordinates. They span the same translations as
    ``periods`` (for example ``result['periods']``) but are usually easier to
    read than the HNF.
    """
    h = hnf(periods)
    vec = np.array(lattice.vectors)
    rows = _reduced_basis(h, vec @ vec.T)
    rows = [r if (r @ vec)[np.flatnonzero(np.abs(r @ vec) > 1e-9)[0]] > 0 else -r for r in rows]
    assert hnf(rows) == h
    rows = np.array(sorted(rows, key=lambda r: (round(float(r @ vec @ vec.T @ r), 9), tuple(-r))))
    return rows, rows @ vec


def site_positions(result, lattice):
    """Cartesian position of every site, ``{label: array}``, one representative per site."""
    vec, pos = np.array(lattice.vectors), np.array(lattice.positions)
    return {site: np.array(c) @ vec + pos[a] for site, (c, a) in result['sites'].items()}


def momenta(periods, lattice=None):
    """The momenta allowed by periodic boundary conditions with these periods.

    These are the k with exp(i k.T) = 1 for every period T. There is one per
    unit cell of the cluster. Returns fractional
    coordinates in the reciprocal primitive basis, each in [0, 1)^dim, and,
    if ``lattice`` is given, also Cartesian k-vectors.
    """
    h = [[Fraction(x) for x in row] for row in hnf(periods)]
    dim = len(h)
    # y . c is an integer for every period c iff y = H^-1 m with m integer; m only matters mod H
    diag = [int(h[j][j]) for j in range(dim)]
    points = set()
    for m in product(*(range(d) for d in diag)):
        y = _solve_upper(h, m)
        points.add(tuple(x - (x.numerator // x.denominator) for x in y))
    fractional = sorted(points)
    if lattice is None:
        return fractional
    reciprocal = 2*np.pi*np.linalg.inv(np.array(lattice.vectors)).T
    return fractional, np.array([[float(x) for x in y] for y in fractional]) @ reciprocal


def _solve_upper(h, m):
    """Exact y with H y = m (H upper triangular)."""
    dim = len(h)
    y = [Fraction(0)]*dim
    for i in range(dim-1, -1, -1):
        y[i] = (Fraction(m[i]) - sum(h[i][j]*y[j] for j in range(i+1, dim))) / h[i][i]
    return y
