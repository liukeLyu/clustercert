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
