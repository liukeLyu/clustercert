import numpy as np
import pytest

import clustercert as cc

matplotlib = pytest.importorskip('matplotlib')
matplotlib.use('Agg')
from matplotlib import pyplot as plt  # noqa: E402

from .test_clustercert import CLUSTERS  # noqa: E402


def drawn_bonds(ax):
    """Line segments that are bonds (the supercell outline is drawn in blue)."""
    return [line for line in ax.lines if line.get_color() != '#2a78d6']


def test_plot_2d_draws_every_bond_once_with_user_labels():
    edges = cc.read_edges(CLUSTERS / 'kagome_12.txt')
    ax = cc.plot(cc.certify(edges, 'kagome'))
    assert len(drawn_bonds(ax)) == len(edges)
    assert {int(t.get_text()) for t in ax.texts} == set(range(12))
    assert sum(line.get_color() == '#2a78d6' for line in ax.lines) == 4      # parallelogram
    plt.close('all')


def test_plot_further_families_and_3d():
    j1, j2 = cc.read_edges(CLUSTERS / 'kagome_12.txt'), cc.read_edges(CLUSTERS / 'kagome_12_j2.txt')
    ax = cc.plot(cc.certify({1: j1, 2: j2}, 'kagome'), further=cc.neighbour_shell(cc.KAGOME, 2))
    assert len(drawn_bonds(ax)) == len(j1) + len(j2)
    p16 = cc.read_edges(CLUSTERS / 'pyrochlore_16.txt')
    ax = cc.plot(cc.certify(p16, 'pyrochlore', find_all=True))
    assert ax.name == '3d' and len(drawn_bonds(ax)) == 48
    assert sum(line.get_color() == '#2a78d6' for line in ax.lines) == 12     # cube
    plt.close('all')


def test_plot_rejects_what_it_cannot_draw():
    with pytest.raises(ValueError, match='certified'):
        cc.plot(cc.certify(cc.read_edges(CLUSTERS / 'kagome_12_impostor.txt'), 'kagome'))
    r16 = cc.certify(cc.read_edges(CLUSTERS / 'pyrochlore_16.txt'), 'pyrochlore')
    with pytest.raises(ValueError, match='double up'):
        cc.plot(r16, further=cc.neighbour_shell(cc.PYROCHLORE, 2))
    ladder = cc.Lattice('ladder', 1, 2, ((0, 1, (0,)), (0, 0, (1,)), (1, 1, (1,))))
    with pytest.raises(ValueError, match='geometry'):
        cc.plot(cc.certify(cc.cluster(ladder, [[4]]).edges, ladder), ladder)


def test_details_draws_unit_cells_matching_the_site_map():
    edges = cc.read_edges(CLUSTERS / 'kagome_12.txt')
    result = cc.certify(edges, 'kagome')
    ax = cc.plot(result, details=True)
    texts = [t.get_text() for t in ax.texts]
    assert sorted(t for t in texts if t.startswith('c=')) == ['c=(0, 0)', 'c=(0, 1)', 'c=(1, 0)', 'c=(1, 1)']
    assert len(drawn_bonds(ax)) == len(edges)       # the unit-cell tiles are blue, so not counted
    with pytest.raises(ValueError, match='2D'):
        cc.plot(cc.certify(cc.read_edges(CLUSTERS / 'pyrochlore_16.txt'), 'pyrochlore'), details=True)
    plt.close('all')


def test_plot_lattice_numbers_every_bond():
    ax = cc.plot_lattice(cc.KAGOME)
    labels = ax.get_legend_handles_labels()[1] + [t.get_text() for t in ax.get_legend().get_texts()]
    assert sum(text.startswith('bond ') for text in labels) >= len(cc.KAGOME.bonds)
    assert {t.get_text() for t in ax.texts} >= {str(k) for k in range(1, 7)}
    with pytest.raises(ValueError):
        cc.plot_lattice(cc.PYROCHLORE)
    plt.close('all')


def test_line_graph_sites_lie_in_unit_cell_zero():
    for lattice in (cc.CHECKERBOARD, cc.line_graph_lattice(cc.HONEYCOMB)):
        frac = np.linalg.solve(np.array(lattice.vectors).T, np.array(lattice.positions).T).T
        assert np.all(frac > -1e-9) and np.all(frac < 1 - 1e-9)


def test_brillouin_zones():
    from clustercert.plot import _brillouin_zone
    for lattice, n_corners, n_edges in ((cc.KAGOME, 6, 6), (cc.SQUARE, 4, 4), (cc.CUBIC, 8, 12),
                                        (cc.FCC, 24, 36), (cc.BCC, 14, 24)):
        b = 2 * np.pi * np.linalg.inv(np.array(lattice.vectors, float)).T
        corners, edges = _brillouin_zone(b)
        assert (len(corners), len(edges)) == (n_corners, n_edges), lattice.name


def test_plot_symmetry_bars_and_mosaic():
    sym = cc.symmetry(cc.certify(cc.read_edges(CLUSTERS / 'kagome_12.txt'), 'kagome'))
    zone, bars = sym.plot(6)
    assert len(bars.patches) == len(sym.sectors) == 10
    assert sorted(int(t.get_text().replace(',', '')) for t in bars.texts if t.get_text()[0].isdigit()) == \
        sorted(sym.block_sizes(6).values())
    assert {t.get_text() for t in zone.texts} == {'Γ', 'M'}
    _, tiles = cc.plot_symmetry(sym, 6, spin_flip=True, mosaic=True)
    assert len(tiles.patches) == 2 * len(sym.sectors)
    area = sum(p.get_width() + 0.008 for p in tiles.patches if p.get_y() < 0.01)        # bottom row spans [0, 1]
    assert np.isclose(area, 1)
    zone, _ = cc.symmetry(cc.certify(cc.read_edges(CLUSTERS / 'pyrochlore_16.txt'), 'pyrochlore')).plot()
    assert zone.name == '3d' and len(zone.lines) == 36         # truncated octahedron
    plt.close('all')


def test_plot_symmetry_elements_kept_and_lost():
    def solid_and_faint(ax):
        mirrors = [line for line in ax.lines if line.get_color() in ('#2a78d6', '#b4b2ad')]
        return sum(line.get_color() == '#2a78d6' for line in mirrors), sum(line.get_color() == '#b4b2ad'
                                                                           for line in mirrors)
    sym = cc.symmetry(cc.certify(cc.read_edges(CLUSTERS / 'kagome_12.txt'), 'kagome'))
    assert solid_and_faint(cc.plot_symmetry_elements(sym)) == (6, 0)
    tilted = [[3, 1], [-1, 3]]                                      # keeps C4, loses every mirror
    sym = cc.symmetry(cc.certify(list(cc.cluster(cc.SQUARE, tilted).edges), 'square', periods=tilted))
    ax = cc.plot_symmetry_elements(sym)
    assert sym.point_group == 'C4' and solid_and_faint(ax) == (0, 4)
    with pytest.raises(ValueError, match='2D'):
        cc.plot_symmetry_elements(cc.symmetry(cc.certify(cc.read_edges(CLUSTERS / 'pyrochlore_16.txt'),
                                                         'pyrochlore')))
    plt.close('all')
