"""Draw a certified cluster with your own site labels (needs matplotlib)."""
from itertools import product

import networkx as nx
import numpy as np

from .certify import _reduce, cluster, hnf
from .lattice import LATTICES, Lattice
from .supercell import period_vectors

__all__ = ['plot', 'plot_lattice', 'plot_symmetry', 'plot_symmetry_elements']

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
            # the label grows away from the arrow, so a slanted arrow does not cut through it
            ha = 'right' if normal[0] < -0.3 else 'left' if normal[0] > 0.3 else 'center'
            va = 'top' if normal[1] < -0.3 else 'bottom' if normal[1] > 0.3 else 'center'
            ax.text(*(base + v / 2 + 0.12 * reach * normal), f'$T_{j + 1}$ = {tuple(h[j])}',
                    color=ACCENT, fontsize=9, ha=ha, va=va, rotation=0)
            tips += [base, base + v, base + v / 2 + 0.5 * reach * normal]
        # a key for the primitive vectors, in its own small axes below the legend
        key = ax.inset_axes([1.04, 0.32, 0.3, 0.35])
        for i, v in enumerate(vec):
            key.annotate('', v, (0, 0), arrowprops=dict(arrowstyle='-|>', shrinkA=0, shrinkB=0, color=INK,
                                                        lw=1.4))
            key.text(*(1.18 * v), f'$a_{i + 1}$', color=INK, fontsize=10, ha='center', va='center')
        key.text(0.5 * reach, -0.3 * reach, 'primitive vectors', color=GREY, fontsize=7, ha='center',
                 va='top')
        ends = np.vstack([np.zeros(2), 1.3 * vec, [0, -0.5 * reach], [reach, -0.5 * reach]])
        key.set_xlim(ends[:, 0].min() - 0.1 * reach, ends[:, 0].max() + 0.1 * reach)
        key.set_ylim(ends[:, 1].min(), ends[:, 1].max())
        key.set_aspect('equal')
        key.axis('off')
        ax.update_datalim(np.vstack(tips))
        ax.autoscale_view()

    n = len(where)
    default = (f'{lattice.name}, {n} sites = {n // lattice.n_basis} unit cells x {lattice.n_basis} sublattices'
               if details else f'{lattice.name}, {n} sites, periods {[list(r) for r in h]}')
    ax.set_title(title if title is not None else default, fontsize=10)
    if further or details:
        ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=8, frameon=False, labelspacing=1.2)
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


def _brillouin_zone(b):
    """Corners of the first Brillouin zone of the reciprocal vectors ``b`` (rows), and its edges as index pairs.

    The zone is the set of k at least as close to 0 as to any reciprocal
    lattice vector G: the intersection of the half-spaces G.k <= |G|^2 / 2.
    Its corners are the points where ``dim`` of these planes meet.
    """
    from itertools import combinations

    dim = len(b)
    gs = sorted((np.array(m) @ b for m in product(range(-2, 3), repeat=dim) if any(m)), key=np.linalg.norm)[:30]
    gs = np.array(gs)
    rhs = np.einsum('ij,ij->i', gs, gs) / 2
    tol = 1e-9 * rhs.max()
    corners, faces = [], []
    for planes in combinations(range(len(gs)), dim):
        m = gs[list(planes)]
        if abs(np.linalg.det(m)) < 1e-9 * np.linalg.norm(m) ** dim:
            continue
        k = np.linalg.solve(m, rhs[list(planes)])
        if np.all(gs @ k <= rhs + tol) and not any(np.allclose(k, c) for c in corners):
            corners.append(k)
            faces.append(set(np.flatnonzero(np.abs(gs @ k - rhs) < tol)))
    corners = np.array(corners)
    edges = [(i, j) for i, j in combinations(range(len(corners)), 2) if len(faces[i] & faces[j]) >= dim - 1]
    return corners, edges


def plot_symmetry(sym, n_up=None, *, spin_flip=False, mosaic=False, axes=None, title=None):
    """Draw the symmetry analysis made by :func:`symmetry`: the momenta and the size of every sector.

    Left: the first Brillouin zone with the allowed momenta, coloured by star
    (momenta of one star have the same spectrum), the named ones labelled.
    Right: one bar per sector, ``momentum.irrep``, as long as its block (see
    :meth:`Symmetry.block_sizes`, with the same ``n_up`` and ``spin_flip``),
    in the colour of its star. The number on the right is the degeneracy of
    every level in that block, ``len(star) * dim``. With ``mosaic=True`` the
    right panel is the whole Hilbert space as one rectangle instead: one
    column per star, as wide as its share of the states, and one tile per
    sector, as tall as its share of the column. 2D lattices are drawn flat,
    3D ones in a 3D view. Returns the two axes.
    """
    import matplotlib.pyplot as plt

    if sym.vectors is None:
        raise ValueError('this Symmetry has no lattice vectors; make it again with clustercert.symmetry')
    b = 2 * np.pi * np.linalg.inv(np.array(sym.vectors, float)).T       # reciprocal primitive vectors, rows
    dim = len(b)
    sizes = sym.block_sizes(n_up, spin_flip)
    rows = []                       # (label, irrep, star index, block size, degeneracy)
    for i, star in enumerate(sym.stars):
        for s in (x for x in sym.sectors if x.momentum == star[0]):
            for sign in (('+', '-') if spin_flip else ('',)):
                rows.append((s.label + sign, s.irrep + sign, i, sizes[s.label + sign], len(star) * s.dim))
    many = len(sym.stars) > len(SUBLATTICE_COLOURS)         # then one colour, and the labels carry the star
    colour = [ACCENT if many else SUBLATTICE_COLOURS[i] for i in range(len(sym.stars))]

    if axes is None:
        height = 4.6 if mosaic else max(4.2, 0.26 * len(rows) + 1.4)
        fig = plt.figure(figsize=(11 if mosaic else 10, height))
        zone_ax = fig.add_subplot(1, 2, 1, projection='3d' if dim == 3 else None)
        info_ax = fig.add_subplot(1, 2, 2)
        if mosaic:
            zone_ax.set_position([0.0, 0.08, 0.3, 0.8])
            info_ax.set_position([0.33, 0.04, 0.66, 0.84])
    else:
        zone_ax, info_ax = axes

    # the zone, and each momentum at its copy closest to Γ
    corners, edges = _brillouin_zone(b)
    if dim == 2:
        order = np.argsort(np.arctan2(corners[:, 1], corners[:, 0]))
        outline = corners[np.append(order, order[0])]
        zone_ax.fill(*outline.T, color=TILE, zorder=0)
        zone_ax.plot(*outline.T, color=GREY, lw=1.2, zorder=1)
    else:
        for i, j in edges:
            zone_ax.plot(*zip(corners[i], corners[j]), color=GREY, lw=1)
    shifts = [np.array(m) @ b for m in product(range(-2, 3), repeat=dim)]
    reach = np.abs(corners).max()
    for i, star in enumerate(sym.stars):
        points = np.array([min((np.array([float(v) for v in y]) @ b + g for g in shifts), key=np.linalg.norm)
                           for y in star])
        style = dict(depthshade=False) if dim == 3 else {}
        zone_ax.scatter(*points.T, s=75, color=colour[i], edgecolor='white', lw=1.5, zorder=3, **style)
        name = next(s.momentum_name for s in sym.sectors if s.momentum == star[0])
        if name.startswith('('):
            continue
        p = points[0]
        away = p / np.linalg.norm(p) if np.linalg.norm(p) > 1e-9 else -np.ones(dim) / np.sqrt(dim)
        zone_ax.text(*(p + 0.14 * reach * away), name, fontsize=11, color=INK, ha='center', va='center')
    zone_ax.set_title(f'{len(sym.momenta)} momenta in {len(sym.stars)} stars\n'
                      'one colour per star: same spectrum', fontsize=9, color=GREY)
    if dim == 3:
        zone_ax.set_box_aspect(np.ptp(corners, axis=0))
        zone_ax.set_axis_off()
    else:
        zone_ax.set_aspect('equal')
        zone_ax.axis('off')

    total = sum(size * degeneracy for *_, size, degeneracy in rows)
    largest = max(r[3] for r in rows) or 1
    space = 'all states' if n_up is None else f'n_up = {n_up}'
    heading = (f'{total:,} states ({space}) in {len(rows)} blocks, the largest {largest:,}\n'
               f'{len(sym.permutations)} symmetries, point group {sym.point_group}')
    if mosaic:
        heading += '; area of a tile = its states, block size × degeneracy'
        _mosaic(info_ax, sym, rows, colour, total)
    else:
        _bars(info_ax, rows, colour, largest)
    info_ax.set_title(title if title is not None else heading, fontsize=9, color=GREY, pad=26 if mosaic else 6)
    return zone_ax, info_ax


def _bars(ax, rows, colour, largest):
    """One bar per sector, with its block size, and the degeneracy of its levels on the right."""
    y = np.arange(len(rows))[::-1]
    for yi, (_, _, i, size, degeneracy) in zip(y, rows):
        ax.barh(yi, size, height=0.72, color=colour[i])
        ax.text(size + 0.015 * largest, yi, f'{size:,}', va='center', fontsize=8, color=INK)
        ax.text(1.45 * largest, yi, f'×{degeneracy}', va='center', ha='right', fontsize=8, color=GREY)
    ax.text(1.45 * largest, y[0] + 0.8, 'degeneracy', va='bottom', ha='right', fontsize=7.5, color=GREY)
    ax.set_yticks(y, [r[0] for r in rows], fontsize=8.5)
    ax.set_ylim(-0.7, len(rows) - 0.1)
    ax.set_xlim(0, 1.45 * largest)
    ax.set_xticks([])
    for side in ax.spines.values():
        side.set_visible(False)
    ax.tick_params(axis='y', length=0)
    ax.set_xlabel('block size: states to diagonalize together', fontsize=8.5, color=GREY)


def _mosaic(ax, sym, rows, colour, total):
    """The Hilbert space as a unit square: a column per star, a tile per sector, areas = numbers of states."""
    from matplotlib.patches import Rectangle

    x, gap = 0.0, 0.004
    for i, star in enumerate(sym.stars):
        mine = [r for r in rows if r[2] == i]
        states = [size * degeneracy for *_, size, degeneracy in mine]
        width = sum(states) / total
        y = 0.0
        for (_, irrep, _, size, degeneracy), n in zip(mine, states):
            height = n / sum(states)
            ax.add_patch(Rectangle((x + gap, y + gap), width - 2 * gap, height - 2 * gap,
                                   facecolor=_tint(colour[i], 0.55), edgecolor=colour[i], lw=1))
            text = f'{irrep}\n{size:,} ×{degeneracy}' if height > 0.09 else f'{irrep}  {size:,} ×{degeneracy}'
            if height > 0.035 and width > 0.05:
                ax.text(x + width / 2, y + height / 2, text, ha='center', va='center', fontsize=7.5, color=INK)
            y += height
        name = mine[0][0].rsplit('.', 1)[0]
        if width > 0.03:
            ax.text(x + width / 2, 1.02, f'{name}\n{len(star)} ' + ('momentum' if len(star) == 1 else 'momenta'),
                    ha='center', va='bottom', fontsize=8.5, color=INK)
        x += width
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.13)
    ax.axis('off')


def plot_symmetry_elements(sym, lattice=None, *, ax=None, title=None):
    """Draw which symmetry elements of a 2D lattice the cluster keeps.

    The lattice is drawn around its most symmetric point. Mirror lines and
    the rotation axis that the cluster keeps are solid; the lattice's other
    ones, lost by the folding (or by further bond families), are faint and
    dashed. ``lattice``: needed only for a lattice that is not built in.
    Returns the axes.
    """
    import matplotlib.pyplot as plt

    from .symmetry import _cartesian, _classify, _origin, space_group

    lattice = lattice if lattice is not None else LATTICES.get(sym.lattice)
    if lattice is None or lattice.dim != 2 or lattice.vectors is None:
        raise ValueError('plot_symmetry_elements draws 2D lattices with geometry; pass the lattice if it is '
                         'not built in')
    full = space_group(lattice)
    centre = _origin([op.linear for op in full], [op.translation for op in full], range(len(full)), 2)
    if centre is None:
        raise ValueError(f'{lattice.name} has glide lines: no point is fixed by every operation')
    a, pos = np.array(lattice.vectors, float), np.array(lattice.positions, float)
    centre = np.array([float(v) for v in centre]) @ a
    kept = {np.array(w).tobytes() for w in sym.point_linear}
    if ax is None:
        ax = plt.figure(figsize=(5.4, 5.4)).add_subplot()

    reach = 2.6 * max(np.linalg.norm(v) for v in a)
    span = int(np.ceil(reach / min(np.linalg.svd(a, compute_uv=False)))) + 2
    sites = []
    for c in product(range(-span, span + 1), repeat=2):
        for i, j, t in lattice.bonds:
            u, v = np.array(c) @ a + pos[i], np.add(c, t) @ a + pos[j]
            if max(np.linalg.norm(u - centre), np.linalg.norm(v - centre)) < reach:
                ax.plot(*zip(u, v), color='#e4e3df', lw=1.4, zorder=0)
        sites += [r for r in np.array(c) @ a + pos if np.linalg.norm(r - centre) < reach]
    ax.scatter(*np.array(sites).T, s=14, color='#c9c8c4', zorder=1)

    rotation = {True: 1, False: 1}
    mirrors = {True: 0, False: 0}
    for w in {op.linear.tobytes(): op.linear for op in full}.values():
        kind, order, axis = _classify(_cartesian(w, a))
        on = w.astype(int).tobytes() in kept
        if kind == 'σ':
            mirrors[on] += 1
            line = np.array([-axis[1], axis[0]])
            ax.plot(*zip(centre - 0.92 * reach * line, centre + 0.92 * reach * line), lw=1.8 if on else 1,
                    color=ACCENT if on else '#b4b2ad', ls='-' if on else (0, (4, 3)), zorder=3 if on else 2)
        elif kind == 'C':
            rotation[on] = max(rotation[on], order)
    marker = {2: 'D', 3: '^', 4: 's', 6: 'h'}
    if rotation[False] > rotation[True]:
        ax.scatter(*centre, marker=marker[rotation[False]], s=420, facecolor='white', edgecolor='#b4b2ad', lw=1.5,
                   zorder=4)
    if rotation[True] > 1:
        ax.scatter(*centre, marker=marker[rotation[True]], s=200, color=INK, zorder=5)
    kept_text = (f'C{rotation[True]} axis' if rotation[True] > 1 else 'no rotation') + \
        f', {mirrors[True]} mirror' + ('' if mirrors[True] == 1 else 's')
    ax.set_title(title if title is not None else
                 f'the cluster keeps {sym.point_group} of {sym.lattice_point_group}: {kept_text}\n'
                 'solid: kept, faint: lost', fontsize=9, color=GREY)
    ax.set_aspect('equal')
    ax.axis('off')
    return ax
