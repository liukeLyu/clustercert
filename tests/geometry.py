"""Clusters built from Cartesian geometry alone, independently of clustercert.

Sites are lattice points folded into the supercell; bonds join sites at the
nearest-neighbour distance under any periodic image.
"""
from itertools import combinations, product

import numpy as np


def geometric_cluster(vectors, positions, periods, bond_length, tol=1e-8):
    vec, pos = np.array(vectors, float), np.array(positions, float)
    cell = np.array(periods, float) @ vec            # supercell vectors, rows
    inv = np.linalg.inv(cell)
    reach = int(np.abs(np.array(periods)).sum()) + 1
    sites = {}
    for c in product(range(-reach, reach+1), repeat=len(vec)):
        for p in pos:
            frac = (np.array(c) @ vec + p) @ inv
            frac = frac - np.floor(frac + tol)
            sites.setdefault(tuple(np.round(frac, 6) % 1), frac)
    sites = list(sites.values())
    images = np.array([np.array(t) @ cell for t in product((-2, -1, 0, 1, 2), repeat=len(vec))])
    r = np.array(sites) @ cell
    edges = []
    for i, j in combinations(range(len(r)), 2):
        hits = np.sum(np.abs(np.linalg.norm(r[j] - r[i] + images, axis=1) - bond_length) < 1e-6)
        if hits > 1:
            raise ValueError('cell too small: two images of the same bond')
        if hits:
            edges.append((i, j))
    return len(sites), edges
