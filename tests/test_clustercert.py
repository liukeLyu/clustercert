import random
from pathlib import Path

import networkx as nx
import numpy as np
import pytest

import clustercert as cc
from clustercert.cli import main, read_edges

from .geometry import geometric_cluster

CLUSTERS = Path(__file__).resolve().parents[1] / 'examples' / 'clusters'
BOND = {'square': 1, 'triangular': 1, 'honeycomb': 1/np.sqrt(3), 'kagome': 0.5, 'cubic': 1,
        'bcc': np.sqrt(3)/2, 'fcc': np.sqrt(2)/2, 'diamond': np.sqrt(3)/4, 'pyrochlore': np.sqrt(2)/4}


def check_witness(edges, result, lattice):
    model = cc.cluster(lattice, result['periods'])
    site = result['sites']
    assert sorted(site.values()) == sorted(model.nodes)
    assert {frozenset((site[a], site[b])) for a, b in edges} == {frozenset(e) for e in model.edges}


def relabel(graph, seed):
    nodes = list(graph)
    random.Random(seed).shuffle(nodes)
    label = {v: i for i, v in enumerate(nodes)}
    return [(label[a], label[b]) for a, b in graph.edges]


@pytest.mark.parametrize('dim, counts', [
    (2, [1, 3, 4, 7, 6, 12, 8, 15]),          # sigma(k), OEIS A000203
    (3, [1, 7, 13, 35, 31, 91, 57, 155]),     # OEIS A001001
])
def test_supercell_counts(dim, counts):
    for index, expected in enumerate(counts, 1):
        hnfs = list(cc.supercell_hnfs(index, dim))
        assert len(hnfs) == len(set(hnfs)) == expected
        assert all(cc.hnf(h) == h for h in hnfs)


def test_hnf_is_basis_independent():
    rng = np.random.default_rng(0)
    for _ in range(200):
        h = cc.hnf(rng.integers(-2, 3, (3, 3)) + 6*np.eye(3, dtype=int))
        u = np.eye(3, dtype=int)
        for _ in range(5):        # random unimodular row operations
            i, j = rng.choice(3, 2, replace=False)
            u[i] += rng.integers(-2, 3) * u[j]
        assert cc.hnf(u @ np.array(h)) == h
    with pytest.raises(ValueError):
        cc.hnf([[1, 2], [2, 4]])


@pytest.mark.parametrize('name', sorted(BOND))
def test_lattice_bonds_match_geometry(name):
    m = cc.LATTICES[name]
    assert cc.lattice_from_geometry(name, m.vectors, m.positions, BOND[name]).bonds == m.bonds


def test_line_graphs():
    rng = random.Random(1)
    for base, target in ((cc.HONEYCOMB, cc.KAGOME), (cc.DIAMOND, cc.PYROCHLORE)):
        line = cc.line_graph_lattice(base)
        for h in rng.sample(list(cc.supercell_hnfs(6, base.dim)), 5):
            a, b = cc.cluster(line, h), cc.cluster(target, h)
            assert (a is None) == (b is None)
            assert a is None or nx.is_isomorphic(a, b)
            if a is not None:
                assert nx.is_isomorphic(a, nx.line_graph(cc.cluster(base, h)))


def test_too_small_cells_are_none():
    assert cc.cluster(cc.PYROCHLORE, np.eye(3, dtype=int)) is None
    assert cc.cluster(cc.SQUARE, 2*np.eye(2, dtype=int)) is None
    assert cc.cluster(cc.SQUARE, 3*np.eye(2, dtype=int)).number_of_edges() == 18


@pytest.mark.parametrize('name, periods', [
    ('square', [[4, 0], [0, 4]]), ('square', [[3, 1], [-1, 3]]), ('triangular', [[4, 0], [1, 3]]),
    ('honeycomb', [[3, 0], [0, 3]]), ('kagome', [[2, 0], [0, 2]]), ('kagome', [[2, 2], [-2, 4]]),
    ('cubic', [[3, 0, 0], [0, 3, 0], [1, 1, 3]]), ('bcc', [[3, 0, 0], [0, 3, 0], [0, 0, 3]]),
    ('fcc', [[3, 0, 0], [0, 3, 0], [1, 0, 3]]), ('diamond', [[2, 0, 0], [0, 2, 0], [0, 0, 2]]),
    ('pyrochlore', [[-1, 1, 1], [1, -1, 1], [1, 1, -1]]), ('pyrochlore', [[2, 0, 0], [0, 1, 1], [0, 1, -2]]),
])
def test_geometric_clusters_certified(name, periods):
    m = cc.LATTICES[name]
    n, edges = geometric_cluster(m.vectors, m.positions, periods, BOND[name])
    result = cc.certify(edges, m, find_all=True)
    assert result['certified'] and n == result['n_sites']
    assert cc.hnf(periods) in result['matching']
    check_witness(edges, result, m)


@pytest.mark.parametrize('file, name', [
    ('pyrochlore_16.txt', 'pyrochlore'), ('pyrochlore_24.txt', 'pyrochlore'), ('pyrochlore_32.txt', 'pyrochlore'),
    ('kagome_12.txt', 'kagome'), ('kagome_27.txt', 'kagome'), ('kagome_36.txt', 'kagome'),
    ('honeycomb_24_nauru.txt', 'honeycomb'), ('triangular_3x6.txt', 'triangular'),
])
def test_example_clusters_certified(file, name):
    edges = read_edges(CLUSTERS / file)
    result = cc.certify(edges, cc.LATTICES[name])
    assert result['certified']
    check_witness(edges, result, cc.LATTICES[name])


def test_relabelled_cluster_recovers_its_periods():
    edges = relabel(cc.cluster(cc.KAGOME, [[2, 1], [0, 4]]), 3)
    result = cc.certify(edges, cc.KAGOME, find_all=True)
    assert result['certified'] and ((2, 1), (0, 4)) in result['matching']


def test_failures_have_reasons():
    slab = cc.certify(read_edges(CLUSTERS / 'pyrochlore_28_slab.txt'), cc.PYROCHLORE)
    assert not slab['certified'] and 'open boundary' in slab['reason']
    assert 'multiple of 3' in cc.certify(cc.cluster(cc.SQUARE, [[4, 0], [0, 4]]).edges, cc.KAGOME)['reason']
    two = nx.disjoint_union(cc.cluster(cc.SQUARE, [[3, 0], [0, 3]]), cc.cluster(cc.SQUARE, [[3, 0], [0, 3]]))
    assert 'disconnected' in cc.certify(two, cc.SQUARE)['reason']


def test_one_misplaced_wrap_bond_rejected():
    # 4x4 square torus with one horizontal wrap bond moved: degrees are still all 4
    g = nx.convert_node_labels_to_integers(cc.cluster(cc.SQUARE, [[4, 0], [0, 4]]), ordering='sorted')
    (a, b), (c, d) = (0, 12), (1, 13)       # two parallel wrap bonds along y
    assert g.has_edge(a, b) and g.has_edge(c, d) and not g.has_edge(a, d) and not g.has_edge(c, b)
    g.remove_edges_from([(a, b), (c, d)])
    g.add_edges_from([(a, d), (c, b)])
    result = cc.certify(g, cc.SQUARE)
    assert not result['certified'] and result['tested'] == 31


def test_locally_pyrochlore_graph_rejected():
    # Line graph of a random 4-regular bipartite graph: 6-regular, two tetrahedra per
    # site, but (for this seed) not a pyrochlore cluster.
    rng = random.Random(0)
    while True:
        stubs = [b for b in range(8) for _ in range(4)]
        rng.shuffle(stubs)
        pairs = [(a, 8+stubs[4*a+k]) for a in range(8) for k in range(4)]
        if len(set(pairs)) == 32 and nx.is_connected(nx.Graph(pairs)):
            break
    result = cc.certify(nx.line_graph(nx.Graph(pairs)), cc.PYROCHLORE)
    assert not result['certified'] and result['tested'] == 155


def test_identify():
    assert cc.identify(read_edges(CLUSTERS / 'kagome_12.txt')) == ['kagome']
    assert 'square' in cc.identify(cc.cluster(cc.SQUARE, [[4, 0], [0, 4]]).edges)


def test_period_vectors_and_momenta():
    edges = read_edges(CLUSTERS / 'pyrochlore_16.txt')
    result = cc.certify(edges, cc.PYROCHLORE, find_all=True)
    cubic = ((1, 1, 1), (0, 2, 0), (0, 0, 2))
    assert len(result['matching']) == 4 and cubic in result['matching']
    rows, cart = cc.period_vectors(cubic, cc.PYROCHLORE)
    assert cc.hnf(rows) == cubic
    assert sorted(np.round(np.sort(np.abs(v)), 6).tolist() for v in cart) == [[0, 0, 1]]*3   # cubic cell
    k, kcart = cc.momenta(cubic, cc.PYROCHLORE)
    assert len(k) == 4 and np.allclose(np.sort(np.linalg.norm(kcart, axis=1)), [0, 2*np.pi, 2*np.pi, 2*np.pi])
    square = cc.certify(cc.cluster(cc.SQUARE, [[4, 0], [0, 4]]).edges, cc.SQUARE)
    _, ks = cc.momenta(square['periods'], cc.SQUARE)
    assert len(ks) == 16 and np.allclose((ks / (np.pi/2)) % 1, 0)
    positions = cc.site_positions(square, cc.SQUARE)
    assert len(positions) == 16


def test_cli(tmp_path, capsys):
    out = tmp_path / 'r.json'
    assert main([str(CLUSTERS / 'kagome_12.txt'), 'kagome', '--json', str(out)]) == 0
    assert 'periodic kagome cluster' in capsys.readouterr().out and out.exists()
    assert main([str(CLUSTERS / 'pyrochlore_28_slab.txt')]) == 1
    assert 'not a periodic cluster of any built-in lattice' in capsys.readouterr().out


# ---- further-neighbour bonds ----

@pytest.mark.parametrize('name, coordination', [
    ('square', [4, 4, 4]), ('triangular', [6, 6, 6]), ('honeycomb', [3, 6, 3]), ('kagome', [4, 4, 6]),
    ('cubic', [6, 12, 8]), ('bcc', [8, 6, 12]), ('fcc', [12, 6, 24]), ('pyrochlore', [6, 12, 12]),
])
def test_neighbour_shell_coordination(name, coordination):
    m = cc.LATTICES[name]
    assert [sorted(set(cc.neighbour_shell(m, k).degrees())) for k in (1, 2, 3)] == [[z] for z in coordination]
    assert cc.neighbour_shell(m, 1).bonds == m.bonds


@pytest.mark.parametrize('name, periods', [
    ('kagome', [[2, 0], [0, 2]]), ('square', [[4, 1], [-1, 4]]), ('triangular', [[4, 0], [0, 4]]),
    ('honeycomb', [[3, 0], [0, 3]]), ('pyrochlore', [[2, 0, 0], [0, 2, 0], [0, 0, 2]]),
])
def test_geometric_j1_j2_clusters_certified(name, periods):
    m = cc.LATTICES[name]
    d1, d2 = cc.neighbour_distances(m, 2)
    _, j1 = geometric_cluster(m.vectors, m.positions, periods, d1)
    _, j2 = geometric_cluster(m.vectors, m.positions, periods, d2)
    result = cc.certify({1: j1, 2: j2}, m)
    assert result['certified'] and cc.hnf(periods) in result['matching']
    shell = cc.neighbour_shell(m, 2)
    assert {frozenset(e) for e in cc.bonds(result, shell)} == {frozenset(e) for e in j2}


def test_j2_singles_out_a_cell():
    j1 = read_edges(CLUSTERS / 'kagome_12.txt')
    j2 = read_edges(CLUSTERS / 'kagome_12_j2.txt')
    result = cc.certify({1: j1, 2: j2}, cc.KAGOME)
    assert result['certified'] and result['matching'] == [((2, 0), (0, 2))]
    assert sum(c['status'] == 'matches' for c in result['supercells']) == 1 and len(result['supercells']) == 4
    with pytest.raises(ValueError, match=r'periods \[\[\[2, 0\], \[0, 2\]\]\]'):
        cc.bonds(cc.certify(j1, cc.KAGOME, find_all=True), cc.neighbour_shell(cc.KAGOME, 2))
    explicit = cc.certify({'J1': j1, 'J2': j2}, {'J1': cc.KAGOME, 'J2': cc.neighbour_shell(cc.KAGOME, 2)})
    assert explicit['certified']


def test_wrong_j2_bond_rejected():
    j1 = read_edges(CLUSTERS / 'kagome_12.txt')
    j2 = read_edges(CLUSTERS / 'kagome_12_j2.txt')
    taken = {frozenset(e) for e in j1 + j2}
    rng = random.Random(0)
    while True:     # swap the ends of two J2 bonds: every degree is unchanged
        (a, b), (c, d) = rng.sample(j2, 2)
        if len({a, b, c, d}) == 4 and not {frozenset((a, d)), frozenset((c, b))} & taken:
            break
    wrong = [e for e in j2 if e not in ((a, b), (c, d))] + [(a, d), (c, b)]
    result = cc.certify({1: j1, 2: wrong}, cc.KAGOME)
    assert not result['certified'] and 'bonded differently' in [c for c in result['supercells']
                                                                 if c['periods'] == ((2, 0), (0, 2))][0]['status']


def test_family_failures_have_reasons():
    p16 = read_edges(CLUSTERS / 'pyrochlore_16.txt')
    with pytest.raises(ValueError, match='double up'):
        cc.bonds(cc.certify(p16, cc.PYROCHLORE, find_all=True), cc.neighbour_shell(cc.PYROCHLORE, 2))
    # 3x3 square torus: third neighbours (two steps along a row) fold onto nearest neighbours
    j1 = list(cc.cluster(cc.SQUARE, [[3, 0], [0, 3]]).edges)
    j2 = list(cc.cluster(cc.neighbour_shell(cc.SQUARE, 2), [[3, 0], [0, 3]]).edges)
    assert 'listed in both' in cc.certify({1: j1, 2: j2 + j1[:1]}, cc.SQUARE)['reason']
    result = cc.certify({1: j1, 3: j2}, cc.SQUARE)          # a 4-regular family 3 that passes the screen
    assert not result['certified']
    assert result['supercells'][0]['status'] == 'family 3 bonds land on family 1 bonds in this supercell'
    assert 'family 3: degrees [0]' in cc.certify({1: j1, 3: []}, cc.SQUARE)['reason']


def test_periods_option():
    j1 = read_edges(CLUSTERS / 'kagome_12.txt')
    assert cc.certify(j1, cc.KAGOME, periods=[[2, 0], [0, 2]])['certified']
    assert cc.certify(j1, cc.KAGOME, periods=[[0, 2], [2, 2]])['certified']        # same cell, other basis
    assert not cc.certify(read_edges(CLUSTERS / 'kagome_27.txt'), cc.KAGOME, periods=[[9, 0], [0, 1]])['certified']
    with pytest.raises(ValueError):
        cc.certify(j1, cc.KAGOME, periods=[[3, 0], [0, 3]])


def test_cli_shells(capsys):
    assert main([str(CLUSTERS / 'kagome_12.txt'), 'kagome', '--shell', '2', str(CLUSTERS / 'kagome_12_j2.txt')]) == 0
    assert 'supercell [[2, 0], [0, 2]]: matches' in capsys.readouterr().out


# ---- bond families up to J3 ----

def test_bond_families_split_mixed_shells():
    kagome = cc.bond_families(cc.KAGOME, 3)
    assert list(kagome) == [1, 2, '3a', '3b']
    assert kagome['3a'].degrees() == [4, 4, 4] and kagome['3b'].degrees() == [2, 2, 2]      # line, hexagon
    pyrochlore = cc.bond_families(cc.PYROCHLORE, 3)
    assert list(pyrochlore) == [1, 2, '3a', '3b']
    assert list(cc.bond_families(cc.SQUARE, 3)) == [1, 2, 3]
    for lattice in (cc.KAGOME, cc.PYROCHLORE):      # the split parts make up the whole shell
        families = cc.bond_families(lattice, 3)
        assert sorted(families['3a'].bonds + families['3b'].bonds) == sorted(cc.neighbour_shell(lattice, 3).bonds)


def test_bfg_bonds_up_to_j3():
    families = cc.bond_families(cc.KAGOME, 3)
    j1 = read_edges(CLUSTERS / 'kagome_36.txt')
    r = cc.certify(j1, 'kagome')
    bonds = {1: j1, 2: cc.bonds(r, families[2]), '3b': cc.bonds(r, families['3b'])}
    assert cc.certify(bonds, 'kagome')['certified']
    with pytest.raises(ValueError, match="'3a', '3b'"):
        cc.certify({1: j1, 3: bonds['3b']}, 'kagome')
    # the 12-site cluster cannot host the across-hexagon J3 in any of its supercells
    r12 = cc.certify(read_edges(CLUSTERS / 'kagome_12.txt'), 'kagome', find_all=True)
    assert all(cc.cluster(families['3b'], h) is None for h in r12['matching'])


def test_bond_families_match_geometry():
    # J3a and J3b of a geometric 6x6 kagome cluster, told apart by the number of J1 steps
    lattice = cc.KAGOME
    d = cc.neighbour_distances(lattice, 3)[-1]
    _, j1 = geometric_cluster(lattice.vectors, lattice.positions, [[6, 0], [0, 6]], 0.5)
    _, j3 = geometric_cluster(lattice.vectors, lattice.positions, [[6, 0], [0, 6]], d)
    g = nx.Graph(j1)
    j3a = [e for e in j3 if nx.shortest_path_length(g, *e) == 2]
    j3b = [e for e in j3 if nx.shortest_path_length(g, *e) == 3]
    assert len(j3a) == 2 * len(j3b)
    assert cc.certify({1: j1, '3a': j3a, '3b': j3b}, 'kagome')['certified']
    assert not cc.certify({1: j1, '3a': j3b, '3b': j3a}, 'kagome')['certified']      # swapped types are caught


def test_j1_relabelling_is_not_enough_for_j2():
    # 4x4 square torus: its J1 graph (the 4D hypercube) has more symmetries than the square lattice,
    # so a valid J1-J2 cluster can be numbered such that the J1-only site map gets J2 wrong.
    h = ((4, 0), (0, 4))
    g1, g2 = cc.cluster(cc.SQUARE, h), cc.cluster(cc.neighbour_shell(cc.SQUARE, 2), h)
    assert nx.is_isomorphic(g1, nx.hypercube_graph(4))
    true_j2 = {frozenset(e) for e in g2.edges}
    maps = list(nx.vf2pp_all_isomorphisms(g1, g1))
    keep = [m for m in maps if {frozenset((m[u], m[v])) for u, v in g2.edges} == true_j2]
    assert (len(maps), len(keep)) == (384, 128)
    label = {v: i for i, v in enumerate(sorted(g1))}
    bad = next(m for m in maps if m not in keep)
    j1 = [(label[u], label[v]) for u, v in g1.edges]
    j2 = [(label[bad[u]], label[bad[v]]) for u, v in g2.edges]
    site = cc.certify(j1, 'square')['sites']
    assert {frozenset((site[i], site[j])) for i, j in j2} != true_j2       # J1-only map fails for J2
    both = cc.certify({1: j1, 2: j2}, 'square')
    site = both['sites']
    assert both['certified'] and {frozenset((site[i], site[j])) for i, j in j2} == true_j2
