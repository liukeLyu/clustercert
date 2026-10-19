"""Draw a certified cluster with your own site labels (needs matplotlib)."""
from itertools import product

import networkx as nx
import numpy as np

from .certify import _reduce, cluster, hnf
from .lattice import LATTICES, Lattice
from .supercell import period_vectors

__all__ = ['plot', 'plot_lattice']

INK, GREY, ACCENT = '#0b0b0b', '#52514e', '#2a78d6'
FAMILY_COLOURS = ['#eb6834', '#1baf7a', '#4a3aa7', '#e87ba4']     # further bond families, in fixed order
SUBLATTICE_COLOURS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
TILE = '#eef4fc'


def _outside(frac):
    """How far fractional supercell coordinates lie outside [0, 1]^d."""
    return float(np.sum(np.maximum(0, -frac) + np.maximum(0, frac - 1)))


def _tint(colour, amount=0.75):
    """``colour`` mixed with white, for a light fill behind dark text."""
    rgb = np.array([int(colour[i:i + 2], 16) for i in (1, 3, 5)]) / 255
    return tuple(rgb + (1 - rgb) * amount)


def plot(result, lattice=None, *, further=(), details=False, ax=None, labels=None, title=None):
    """Draw the cluster certified by ``result``, labelled with your own site labels.

    Every site is drawn once, folded into the supercell (blue outline). Every
    bond is drawn once, from a site to the true position of its partner;
    when the bond wraps around the boundary, that end is a hollow marker
    carrying the partner's label. ``further``: lattices for further bond
    families (for example ``neighbour_shell(KAGOME, 2)``), drawn dashed.
    2D clusters are drawn flat, 3D clusters in a 3D view. Returns the axes.

    ``result`` can also be a graph made by :func:`cluster`, to see a folding
    before certifying anything; its sites are then unlabelled by default.

    ``details=True`` shows the coordinates behind the certificate. Every
    site sits at its true position in its unit cell c, exactly as in
    ``result['sites']``: the supercell is drawn as its unit cells, each
    labelled with c, sites are coloured by sublattice, and arrows show the
    primitive vectors.
    """
    import matplotlib.pyplot as plt

    if isinstance(result, nx.Graph):       # a folded lattice from cluster(): it is its own site map
        if 'periods' not in result.graph:
            raise ValueError('plot takes a certify() result or a graph made by cluster()')
        result = dict(certified=True, lattice=result.graph['lattice'], periods=result.graph['periods'],
                      sites={node: node for node in result})
        labels = False if labels is None else labels
    labels = True if labels is None else labels
    if not result.get('certified'):
        raise ValueError(f'only a certified cluster can be drawn: {result.get("reason")}')
    lattice = lattice or LATTICES.get(result['lattice'])
    if lattice is None or lattice.vectors is None:
        raise ValueError('drawing needs a lattice with geometry (vectors and positions)')
    if lattice.dim not in (2, 3):
        raise ValueError('only 2D and 3D clusters can be drawn')
    if details and lattice.dim != 2:
        raise ValueError('details=True is for 2D clusters; a 3D view with every unit cell labelled is unreadable')
    further = [further] if isinstance(further, Lattice) else list(further)
    h = hnf(result['periods'])
    for m in further:
        if cluster(m, h) is None:
            raise ValueError(f'{m.name} bonds close on themselves or double up in this supercell')

    vec, pos = np.array(lattice.vectors, float), np.array(lattice.positions, float)
    _, cell = period_vectors(h, lattice)
    inverse = np.linalg.inv(cell)
    label = {node: site for site, node in result['sites'].items()}
    box = [h[j][j] for j in range(lattice.dim)]       # unit cells c of the site map: prod range(box_j)
    where = {}
    for c, a in label:
        r = np.array(c) @ vec + pos[a]
        where[(c, a)] = r if details else r - np.floor(r @ inverse + 1e-9) @ cell

    def outside(end, cell_of_end):
        """How far a bond's far end lies outside the drawn region (supercell, or the unit cells)."""
        if details:
            return float(sum(max(0, -x) + max(0, x - d + 1) for x, d in zip(cell_of_end, box)))
        return _outside(end @ inverse)

    three_d = lattice.dim == 3
    if ax is None:
        fig = plt.figure(figsize=(6, 5.5) if three_d else (5.5, 5))
        ax = fig.add_subplot(projection='3d') if three_d else fig.add_subplot()
    if three_d:
        ax.computed_zorder = False      # draw in our order: bonds, then sites, then labels on top
        ax.view_init(elev=22, azim=-58)

    def parallelotope(origin, edges, **style):
        corners = [origin + np.array(k) @ edges for k in product((0, 1), repeat=lattice.dim)]
        for i, j in product(range(len(corners)), repeat=2):    # edges join corners that differ in one index
            if i < j and bin(i ^ j).count('1') == 1:
                ax.plot(*zip(corners[i], corners[j]), zorder=0, **style)
        return corners

    if details:
        corners = []
        for c in product(*(range(d) for d in box)):
            corners += parallelotope(np.array(c) @ vec, vec, color=ACCENT, lw=0.8)
            if not three_d:
                ax.fill(*np.array([np.array(c) @ vec + np.array(k) @ vec
                                   for k in ((0, 0), (1, 0), (1, 1), (0, 1))]).T, color=TILE, zorder=-1)
            centre = np.array(c) @ vec + vec.sum(0) / 2
            ax.text(*centre, f'c={c}', color=ACCENT, fontsize=7, ha='center', va='center', zorder=0)
    else:
        corners = parallelotope(0 * cell[0], cell, color=ACCENT, lw=1.6)

    images = {}
    families = [(lattice, dict(color=GREY, lw=1.6), 'nearest-neighbour bonds')]
    families += [(m, dict(color=FAMILY_COLOURS[i % len(FAMILY_COLOURS)], lw=1.3, ls='--'), m.name)
                 for i, m in enumerate(further)]
    for m, style, name in families:
        first = True
        for (c, a), start in where.items():
            for a0, b, t in m.bonds:
                if a0 != a:
                    continue
                bond = np.array(t) @ vec + pos[b] - pos[a]
                node, partner = (c, a), (_reduce([x + y for x, y in zip(c, t)], h), b)
                # draw a wrapped bond from whichever end puts its image closest to the supercell
                ends = [(start, start + bond, partner, tuple(np.add(c, t))),
                        (where[partner], where[partner] - bond, node, tuple(np.subtract(partner[0], t)))]
                origin, end, partner, _ = min(ends, key=lambda e: outside(e[1], e[3]))
                ax.plot(*zip(origin, end), zorder=1, label=name if first and further else None, **style)
                first = False
                if np.abs(end - where[partner]).max() > 1e-9:
                    images[tuple(np.round(end, 6))] = (end, partner)

    size = (150 if three_d else 220) if labels else 40
    shade = dict(depthshade=False) if three_d else {}
    for end, partner in images.values():
        ax.scatter(*end, s=size, facecolor='white', edgecolor=GREY, lw=1, linestyle='--', zorder=2, **shade)
        if labels:
            ax.text(*end, str(label[partner]), color=GREY, fontsize=7, ha='center', va='center', zorder=4)
    points = np.array(list(where.values()))
    if details:
        for a in range(lattice.n_basis):
            colour = SUBLATTICE_COLOURS[a % len(SUBLATTICE_COLOURS)]
            mine = np.array([r for (c, b), r in where.items() if b == a])
            ax.scatter(*mine.T, s=size, facecolor=_tint(colour) if labels else colour, edgecolor=colour, lw=1.6,
                       zorder=3, label=f'sublattice {a}', **shade)
    else:
        ax.scatter(*points.T, s=size, facecolor='white' if labels else INK, edgecolor=INK, lw=1.2, zorder=3,
                   **shade)
    if labels:
        for node, r in where.items():
            ax.text(*r, str(label[node]), color=INK, fontsize=8, ha='center', va='center', zorder=5)

    if details:
        # the periods as dimension lines along the supercell, just outside everything drawn
        reach = max(np.linalg.norm(v) for v in vec)
        rows = [np.array(row) @ vec for row in h]
        drawn = np.vstack([points, *[e for e, _ in images.values()]])
        tips = []
        for j, v in enumerate(rows):
            other = sum(rows) - v
            normal = np.array([-v[1], v[0]]) / np.linalg.norm(v)
            normal = -normal if normal @ other > 0 else normal          # away from the supercell
            gap = max(0.0, float(np.max(drawn @ normal)))         # how far the drawing reaches past the edge
            base = (gap + 0.35 * reach) * normal
            ax.annotate('', base + v, base, arrowprops=dict(arrowstyle='<|-|>', shrinkA=0, shrinkB=0,
                                                            color=ACCENT, lw=1.6))
            for end in (0 * v, v):         # thin ticks back to the supercell corners
                ax.plot(*zip(end + 0.1 * reach * normal, base + end + 0.1 * reach * normal), color=ACCENT,
                        lw=0.7, ls=':')
            ax.text(*(base + v / 2 + 0.28 * reach * normal), f'$T_{j + 1}$ = {tuple(h[j])}',
                    color=ACCENT, fontsize=9, ha='center', va='center', rotation=0)
            tips += [base, base + v, base + v / 2 + 0.5 * reach * normal]
        # a key for the primitive vectors, below and left of everything drawn
        drawn = np.vstack([points, *[e for e, _ in images.values()], *corners, *tips])
        start = drawn.min(0) - 0.3 * reach
        start[1 if not three_d else 2] -= 1.1 * reach
        for i, v in enumerate(vec):
            if three_d:
                ax.quiver(*start, *v, arrow_length_ratio=0.15, color=INK, lw=1.4)
            else:
                ax.annotate('', start + v, start, arrowprops=dict(arrowstyle='-|>', shrinkA=0, shrinkB=0,
                                                                  color=INK, lw=1.4))
            ax.text(*(start + 1.18 * v), f'$a_{i + 1}$', color=INK, fontsize=10, ha='center', va='center')
        ax.text(*(start - 0.12 * reach * np.eye(lattice.dim)[0]), 'primitive\nvectors', color=GREY, fontsize=7,
                ha='right', va='center')
        ax.update_datalim(np.vstack([start - 0.6 * reach * np.eye(2)[0], start + vec * 1.25, *tips]))
        ax.autoscale_view()

    n = len(where)
    default = (f'{lattice.name}, {n} sites = {n // lattice.n_basis} unit cells x {lattice.n_basis} sublattices'
               if details else f'{lattice.name}, {n} sites, periods {[list(r) for r in h]}')
    ax.set_title(title if title is not None else default, fontsize=10)
    if further or details:
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=8, frameon=False)
    if three_d:
        span = np.ptp(np.vstack([points, *corners]), axis=0)
        ax.set_box_aspect(span)
        ax.set_axis_off()
    else:
        ax.set_aspect('equal')
        ax.axis('off')
    return ax


def plot_lattice(lattice, *, ax=None, title=None):
    """Draw one unit cell of a 2D lattice: what ``lattice.describe()`` says, as a picture.

    Unit cell 0 is shaded, with arrows for the primitive vectors. Its
    sublattices are numbered and coloured (the legend gives their positions).
    Every bond ``(a, b, t)`` of the unit cell is drawn from sublattice a in
    unit cell 0 to sublattice b in unit cell t, and labelled with its
    ``(a, b, t)``. Partner sites outside unit cell 0 are hollow, marked with
    their unit cell t. Returns the axes.
    """
    import matplotlib.pyplot as plt

    if lattice.vectors is None or lattice.dim != 2:
        raise ValueError('plot_lattice draws 2D lattices with geometry (vectors and positions)')
    vec, pos = np.array(lattice.vectors, float), np.array(lattice.positions, float)
    if ax is None:
        ax = plt.figure(figsize=(6, 5)).add_subplot()
    colour = [SUBLATTICE_COLOURS[a % len(SUBLATTICE_COLOURS)] for a in range(lattice.n_basis)]
    reach = max(np.linalg.norm(v) for v in vec)

    for c in product(range(-1, 3), repeat=2):          # the surrounding lattice, faint
        for a, b, t in lattice.bonds:
            ax.plot(*zip(np.array(c) @ vec + pos[a], (np.add(c, t)) @ vec + pos[b]), color='#d9d8d3', lw=1,
                    zorder=0)
    ax.fill(*np.array([np.array(k) @ vec for k in ((0, 0), (1, 0), (1, 1), (0, 1))]).T, color=TILE, zorder=0)
    ax.plot(*np.array([np.array(k) @ vec for k in ((0, 0), (1, 0), (1, 1), (0, 1), (0, 0))]).T, color=ACCENT,
            lw=1.4, zorder=0)
    ax.text(*(vec.sum(0) * 0.62), 'unit cell\nc = (0, 0)', color=ACCENT, fontsize=8, ha='center', va='center')

    partners, handles, placed = {}, [], [np.asarray(p) for p in pos]
    for k, (a, b, t) in enumerate(lattice.bonds, 1):
        start, end = pos[a], np.array(t) @ vec + pos[b]
        ax.plot(*zip(start, end), color=GREY, lw=2.2, zorder=1)
        for f in (0.5, 0.32, 0.68, 0.2, 0.8):          # number at the midpoint, unless another number is there
            spot = start + f * (end - start)
            if all(np.linalg.norm(spot - q) > 0.12 * reach for q in placed):
                break
        placed.append(spot)
        ax.text(*spot, str(k), color='white', fontsize=7, ha='center', va='center', zorder=7,
                bbox=dict(boxstyle='circle,pad=0.25', facecolor=GREY, edgecolor='none'))
        handles.append(plt.Line2D([], [], ls='none', label=f'bond {k}: ({a}, {b}, {tuple(t)})'))
        if any(t):
            partners[(tuple(t), b)] = end
    for (t, b), r in partners.items():
        ax.scatter(*r, s=170, facecolor='white', edgecolor=colour[b], lw=1.4, linestyle='--', zorder=2)
        ax.text(*r, str(b), color=INK, fontsize=8, ha='center', va='center', zorder=3)
        ax.text(*(r + np.array([0, -0.11 * reach])), f'c={t}', color=GREY, fontsize=6.5, ha='center', va='top',
                zorder=3)
    for a in range(lattice.n_basis):
        ax.scatter(*pos[a], s=220, facecolor=_tint(colour[a]), edgecolor=colour[a], lw=1.8, zorder=3,
                   label=f'sublattice {a} at ({", ".join(f"{x:.3g}" for x in pos[a])})')
        ax.text(*pos[a], str(a), color=INK, fontsize=9, ha='center', va='center', zorder=4)
    for i, v in enumerate(vec):
        ax.annotate('', v, (0, 0), arrowprops=dict(arrowstyle='-|>', shrinkA=0, shrinkB=0, color=INK, lw=1.6),
                    zorder=2)
        ax.text(*(1.1 * v), f'$a_{i + 1}$', color=INK, fontsize=11, ha='center', va='center')

    ax.set_title(title if title is not None else
                 f'{lattice.name}: one unit cell, {lattice.n_basis} sublattices, {len(lattice.bonds)} bonds',
                 fontsize=10)
    sites = ax.get_legend_handles_labels()[0]
    ax.legend(handles=sites + handles, loc='upper left', bbox_to_anchor=(1, 1), fontsize=8, frameon=False,
              handletextpad=0.4)
    ax.text(1.02, 0.02, 'bond (a, b, t): sublattice a in unit cell c\nto sublattice b in unit cell c + t',
            transform=ax.transAxes, fontsize=7.5, color=GREY, va='bottom')
    every = np.vstack([pos, *partners.values(), vec, [0, 0], vec.sum(0)])
    low, high = every.min(0) - 0.35 * reach, every.max(0) + 0.35 * reach
    ax.set_xlim(low[0], high[0])
    ax.set_ylim(low[1], high[1])
    ax.set_aspect('equal')
    ax.axis('off')
    return ax
