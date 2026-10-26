"""Symmetry sectors, checked by brute force on the full (fixed-S^z) Hilbert space of small clusters."""
from math import comb
from pathlib import Path

import numpy as np
import pytest

import clustercert as cc
from clustercert.symmetry import _origin, _permute_states, space_group

EXAMPLES = Path(__file__).resolve().parents[1] / 'examples' / 'clusters'

# (lattice, periods): one small cluster per built-in lattice, plus tilted (non-diagonal HNF) supercells
CLUSTERS = [
    ('square', [[4, 0], [0, 4]]), ('square', [[3, 1], [-1, 3]]), ('triangular', [[2, 2], [-2, 4]]),
    ('honeycomb', [[2, 0], [0, 2]]), ('honeycomb', [[1, 1], [0, 3]]), ('kagome', [[2, 0], [0, 2]]),
    ('kagome', [[1, 1], [0, 3]]), ('checkerboard', [[2, 2], [0, 4]]),
    ('cubic', [[2, 0, 2], [0, 2, 2], [0, 0, 4]]), ('bcc', [[2, 0, 2], [0, 2, 2], [0, 0, 4]]),
    ('fcc', [[1, 1, 1], [0, 4, 0], [0, 0, 4]]), ('diamond', [[2, 0, 0], [0, 2, 0], [0, 0, 2]]),
    ('diamond', [[1, 1, 1], [0, 2, 0], [0, 0, 2]]), ('pyrochlore', [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]),
]


def certified(name, periods, labels=None):
    """``name``: a built-in lattice name, or a Lattice."""
    lattice = cc.LATTICES[name] if isinstance(name, str) else name
    edges = list(cc.cluster(lattice, periods).edges)
    if labels:          # the user's own labels, unrelated to (cell, sublattice)
        nodes = sorted({v for e in edges for v in e})
        rename = {v: f's{i}' for i, v in enumerate(np.random.default_rng(1).permutation(len(nodes)))}
        rename = {v: rename[i] for i, v in zip(rename, nodes)}
        edges = [(rename[u], rename[v]) for u, v in edges]
    result = cc.certify(edges, lattice, periods=periods)
    assert result['certified']
    return result, edges


def heisenberg(sym, edges, states, rng):
    """XXZ model with random couplings on the bonds, symmetric because all bonds get the same couplings."""
    jxy, jz = rng.normal(size=2)
    where = {s: i for i, s in enumerate(sym.sites)}
    h = np.zeros((len(states), len(states)))
    for u, v in edges:
        i, j = where[u], where[v]
        bi, bj = (states >> i) & 1, (states >> j) & 1
        h[np.arange(len(states)), np.arange(len(states))] += jz * (bi - 0.5) * (bj - 0.5)
        flip = bi != bj
        target = np.searchsorted(states, states[flip] ^ ((1 << i) | (1 << j)))
        h[target, np.flatnonzero(flip)] += jxy / 2
    return h


def sector_space(p):
    vals, vecs = np.linalg.eigh(p)
    assert np.allclose(vals, np.round(vals), atol=1e-8)
    return vecs[:, vals > 0.5]


def brute_force(name, periods, n_up, labels=False):
    result, edges = certified(name, periods, labels)
    sym = cc.symmetry(result, None if isinstance(name, str) else name)
    n = len(sym.sites)
    sizes = sym.block_sizes(n_up)
    states = np.array([x for x in range(2 ** n) if n_up is None or bin(x).count('1') == n_up])
    h = heisenberg(sym, edges, states, np.random.default_rng(0))
    for g in range(len(sym.permutations)):          # H commutes with every symmetry
        rows = np.searchsorted(states, _permute_states(states, sym.permutations[g]))
        assert np.allclose(h[np.ix_(rows, rows)], h)
    total = np.zeros((len(states), len(states)), complex)
    spectrum = []
    for s in sym.sectors:
        for k in s.star:
            _, p = sym.projector(s.label, n_up, momentum=k)
            assert np.allclose(p @ p, p) and np.allclose(p, p.conj().T)
            assert np.isclose(np.trace(p).real, sizes[s.label] * s.dim)
            total += p
        q = sector_space(sym.projector(s.label, n_up)[1])
        assert q.shape[1] == sizes[s.label] * s.dim
        block = np.linalg.eigvalsh(q.conj().T @ h @ q)
        spectrum += list(block) * len(s.star)
    assert np.allclose(total, np.eye(len(states)))     # projectors that sum to 1 are mutually orthogonal
    assert np.allclose(np.sort(spectrum), np.linalg.eigvalsh(h))
    return sym, sizes


def small_n_up(n, limit=250):
    return max(m for m in range(n // 2 + 1) if comb(n, m) <= limit)


@pytest.mark.parametrize('name, periods', CLUSTERS)
def test_brute_force_sectors(name, periods):
    n = cc.LATTICES[name].n_basis * abs(round(np.linalg.det(periods)))
    brute_force(name, periods, small_n_up(n), labels=True)


@pytest.mark.parametrize('name, periods', [('square', [[3, 1], [-1, 3]]), ('kagome', [[1, 1], [0, 3]])])
def test_brute_force_whole_hilbert_space(name, periods):
    sym, sizes = brute_force(name, periods, None)
    assert sum(sizes[s.label] * s.dim * len(s.star) for s in sym.sectors) == 2 ** len(sym.sites)


@pytest.mark.slow
@pytest.mark.parametrize('name, periods', CLUSTERS)
def test_brute_force_sectors_larger_block(name, periods):
    n = cc.LATTICES[name].n_basis * abs(round(np.linalg.det(periods)))
    brute_force(name, periods, small_n_up(n, 2000))


@pytest.mark.parametrize('name, periods', CLUSTERS)
def test_group_sanity(name, periods):
    sym = cc.symmetry(certified(name, periods)[0], automorphisms=False)
    n_cells = len(sym.translation_cells)
    n_point = len(sym.point_operations)
    assert len(sym.permutations) == n_cells * n_point
    assert len({p.tobytes() for p in sym.permutations}) == len(sym.permutations)
    assert sum(len(star) for star in sym.stars) == n_cells == len(sym.momenta)
    for star in sym.stars:
        secs = [s for s in sym.sectors if s.momentum == star[0]]
        order = len(secs[0].elements) // n_cells           # little co-group
        assert len(star) * order == n_point
        assert sum(s.dim ** 2 for s in secs) == order
        for a in secs:              # characters orthonormal over the little group
            for b in secs:
                overlap = np.vdot(b.characters, a.characters) / len(a.elements)
                assert np.isclose(overlap, a is b)
    for n_up in (None, 1, 2):
        total = comb(len(sym.sites), n_up) if n_up is not None else 2 ** len(sym.sites)
        sizes = sym.block_sizes(n_up)
        assert sum(sizes[s.label] * s.dim * len(s.star) for s in sym.sectors) == total
    gens = sym.generators()
    assert 1 <= len(gens) <= 6


@pytest.mark.parametrize('name, order, symmorphic', [
    ('square', 8, True), ('triangular', 12, True), ('honeycomb', 12, True), ('kagome', 12, True),
    ('checkerboard', 8, True), ('cubic', 48, True), ('bcc', 48, True), ('fcc', 48, True),
    ('diamond', 48, False), ('pyrochlore', 48, False),
])
def test_space_group_of_lattice(name, order, symmorphic):
    ops = space_group(cc.LATTICES[name])
    assert len(ops) == len({op.linear.tobytes() for op in ops}) == order    # primitive cells: one t per W
    ws = [op.linear for op in ops]
    ts = [op.translation for op in ops]
    assert (_origin(ws, ts, range(len(ops)), cc.LATTICES[name].dim) is not None) == symmorphic


@pytest.mark.parametrize('periods, kept, stars, point_group', [
    ([[2, 1], [-2, 3]], 4, [1, 4, 2, 1], 'C2v'),
    ([[5, 0], [1, 2]], 2, [1, 1, 2, 2, 2, 2], 'C2'),
    ([[2, 2], [-2, 4]], 12, [1, 3, 6, 2], 'C6v'),
])
def test_kagome_clusters(periods, kept, stars, point_group):
    sym = cc.symmetry(certified('kagome', periods)[0])
    assert len(sym.point_operations) == kept and sym.point_group == point_group
    assert sorted(len(s) for s in sym.stars) == sorted(stars)
    cogroups = {s.momentum_name: s.little_cogroup for s in sym.sectors}
    if point_group == 'C6v':
        assert cogroups == {'Γ': 'C6v', 'M': 'C2v', '(1/6, 1/3)': 'Cs', 'K': 'C3v'}
        assert [s.irrep for s in sym.sectors if s.momentum_name == 'Γ'] == ['A1', 'A2', 'B1', 'B2', 'E1', 'E2']
        assert sym.n_automorphisms == 144
    if point_group == 'C2v':
        assert cogroups['Γ'] == cogroups['M'] == 'C2v'


def test_one_magnon_labels():
    """One flipped spin on kagome: the three sublattices at Γ carry A1 + E2 of C6v (C2 fixes every sublattice)."""
    sym = cc.symmetry(certified('kagome', [[2, 2], [-2, 4]])[0])
    sizes = sym.block_sizes(n_up=1)
    assert {k: v for k, v in sizes.items() if k.startswith('Γ') and v} == {'Γ.A1': 1, 'Γ.E2': 1}
    # pyrochlore: the 4 corners of a tetrahedron carry A1g + T2g of Oh (C4 permutes them cyclically: 1 - 1 = 0)
    pyro = cc.symmetry(certified('pyrochlore', [[-1, 1, 1], [1, -1, 1], [1, 1, -1]])[0])
    assert {k: v for k, v in pyro.block_sizes(n_up=1).items() if k.startswith('Γ') and v} == {'Γ.A1g': 1, 'Γ.T2g': 1}
    # diamond: the two sublattices give A1g + A2u (odd under every operation that swaps them, such as inversion)
    diamond = cc.symmetry(certified('diamond', [[2, 0, 0], [0, 2, 0], [0, 0, 2]])[0])
    assert {k: v for k, v in diamond.block_sizes(n_up=1).items() if k.startswith('Γ') and v} == {
        'Γ.A1g': 1, 'Γ.A2u': 1}
    square = cc.symmetry(certified('square', [[4, 0], [0, 4]])[0])
    assert {k for k, v in square.block_sizes(n_up=1).items() if v} == {
        'Γ.A1', 'M.A1', 'X.A1', '(0, 1/4).A\'', '(1/4, 1/4).A\'', '(1/4, 1/2).A\''}


def test_cubic_names():
    sym = cc.symmetry(certified('cubic', [[2, 0, 2], [0, 2, 2], [0, 0, 4]])[0], automorphisms=False)
    gamma = {s.irrep for s in sym.sectors if s.momentum_name == 'Γ'}
    assert sym.point_group == 'Oh' and gamma == {f'{x}{p}' for x in ('A1', 'A2', 'E', 'T1', 'T2') for p in 'gu'}
    assert {s.momentum_name for s in sym.sectors} >= {'Γ', 'X', 'M', 'R'} - {'M'}


def test_diamond_x_point_sticks_together():
    """Diamond is non-symmorphic: at X every irrep is 2-dimensional (bands stick together in pairs)."""
    sym = cc.symmetry(certified('diamond', [[2, 0, 0], [0, 2, 0], [0, 0, 2]])[0], automorphisms=False)
    x = [s for s in sym.sectors if s.momentum_name == 'X']
    assert x and x[0].little_cogroup == 'D4h' and {s.dim for s in x} == {2} and len(x) == 4
    assert {s.momentum_name for s in sym.sectors} == {'Γ', 'X', 'L'}
    gamma = {s.irrep for s in sym.sectors if s.momentum_name == 'Γ'}
    assert sym.point_group == 'Oh' and 'A1g' in gamma and len(gamma) == 10


def test_spin_flip_parity():
    result, edges = certified('kagome', [[2, 0], [0, 2]])
    sym = cc.symmetry(result)
    split = sym.block_sizes(n_up=6, spin_flip=True)
    plain = sym.block_sizes(n_up=6)
    assert all(split[k + '+'] + split[k + '-'] == v for k, v in plain.items())
    states = np.array([x for x in range(2 ** 12) if bin(x).count('1') == 6])
    flip = np.searchsorted(states, states ^ (2 ** 12 - 1))
    for s in sym.sectors:
        p = sym.projector(s.label, 6)[1]
        even = np.trace(p + p[:, flip]).real / 2 / s.dim
        assert np.isclose(even, split[s.label + '+'])


def test_further_neighbour_families_restrict_the_group():
    j1 = cc.read_edges(EXAMPLES / 'kagome_12.txt')
    j2 = cc.read_edges(EXAMPLES / 'kagome_12_j2.txt')
    both = cc.symmetry(cc.certify({1: j1, 2: j2}, 'kagome'))
    assert len(both.point_operations) == 12 and both.n_automorphisms == 48


def test_summary_and_cli(capsys):
    from clustercert.cli import main
    sym = cc.symmetry(certified('kagome', [[2, 2], [-2, 4]])[0])
    text = sym.summary(n_up=1)
    assert 'C6v' in text and 'K' in text and 'block sizes' in text
    assert main([str(EXAMPLES / 'kagome_36.txt'), 'kagome', '--symmetry', '--n-up', '1']) == 0
    assert 'C6v              A1 A2 B1 B2 E1(2) E2(2)' in capsys.readouterr().out


def test_user_defined_lattice():
    """A lattice that is not built in (Lieb: a corner and two edge-centre sites per square) works the same way."""
    lieb = cc.lattice_from_geometry('Lieb', vectors=[[1, 0], [0, 1]], positions=[[0, 0], [0.5, 0], [0, 0.5]],
                                    bond_length=0.5)
    assert len(space_group(lieb)) == 8
    edges = list(cc.cluster(lieb, [[2, 0], [0, 2]]).edges)
    sym = cc.symmetry(cc.certify(edges, lieb), lieb)
    assert sym.point_group == 'C4v' and len(sym.permutations) == 32
    sizes = sym.block_sizes(n_up=1)
    # corner: A1; the two edge centres are swapped by C4 and the diagonal mirrors, kept by the axial ones: A1 + B1
    assert {k: v for k, v in sizes.items() if k.startswith('Γ') and v} == {'Γ.A1': 2, 'Γ.B1': 1}


def test_user_defined_lattices_brute_force():
    lieb = cc.lattice_from_geometry('Lieb', vectors=[[1, 0], [0, 1]], positions=[[0, 0], [0.5, 0], [0, 0.5]],
                                    bond_length=0.5)
    brute_force(lieb, [[2, 0], [0, 2]], 4)
    # a square lattice described with a doubled (non-primitive) cell: the extra translation appears as a point
    # operation with W = 1, names are numbered, and the sectors are still right
    doubled = cc.lattice_from_geometry('square, 2-site cell', vectors=[[2, 0], [0, 1]], positions=[[0, 0], [1, 0]],
                                       bond_length=1)
    sym, _ = brute_force(doubled, [[2, 0], [0, 3]], 3)
    assert sym.point_group.startswith('order') and len(sym.point_operations) == 8


def test_terminology_worked_example():
    """The one-magnon table of docs/TERMINOLOGY.md, section 6."""
    edges = cc.read_edges(EXAMPLES / 'kagome_36.txt')
    sym = cc.symmetry(cc.certify(edges, 'kagome'))
    where = {s: i for i, s in enumerate(sym.sites)}
    h = np.zeros((36, 36))           # Heisenberg model, one up spin, energy relative to the polarized state
    for u, v in edges:
        i, j = where[u], where[v]
        h[i, j] = h[j, i] = 0.5
        h[i, i] -= 0.5
        h[j, j] -= 0.5
    table = {}
    for s in sym.sectors:
        q = sector_space(sym.projector(s.label, 1)[1])
        if q.shape[1]:
            table[s.label] = sorted(set(np.round(np.linalg.eigvalsh(q.conj().T @ h @ q), 6)))
    assert table == {'Γ.A1': [0], 'Γ.E2': [-3], 'M.A1': [-2], 'M.B1': [-3], 'M.B2': [-1],
                     "(1/6, 1/3).A'": [-3, -0.5], '(1/6, 1/3).A"': [-2.5], 'K.A1': [-3], 'K.E': [-1.5]}
    doc = (Path(__file__).resolve().parents[1] / 'docs' / 'TERMINOLOGY.md').read_text()
    printed = sym.summary(n_up=1).splitlines()
    for block in doc.split('```text\n')[1:]:          # every output shown on the page is printed by summary()
        assert all(line in printed for line in block.split('```')[0].splitlines()), block


def test_more_symmetric_supercell_is_pointed_out():
    edges = cc.read_edges(EXAMPLES / 'pyrochlore_32.txt')
    sym = cc.symmetry(cc.certify(edges, 'pyrochlore', find_all=True))
    assert sym.better_supercell == (((2, 0, 0), (0, 2, 0), (0, 0, 2)), 48)
    cubic = cc.symmetry(cc.certify(edges, 'pyrochlore', periods=sym.better_supercell[0]))
    assert cubic.point_group == 'Oh' and cubic.better_supercell is None
    assert len(cubic.permutations) == cubic.n_automorphisms == 384
