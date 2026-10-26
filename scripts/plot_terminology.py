"""Figures for docs/TERMINOLOGY.md, all computed from the 36-site kagome cluster in examples/clusters.

Run from the repository root: python scripts/plot_terminology.py -> docs/terminology_*.png (a few seconds).
"""
from itertools import product
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import clustercert as cc  # noqa: E402
from clustercert.symmetry import _cartesian, _classify  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
INK, GREY, LIGHT = '#0b0b0b', '#52514e', '#c9c8c4'
STAR_COLOURS = ['#2a78d6', '#eb6834', '#1baf7a', '#4a3aa7']      # fixed order: one per star
DPI = 110

sym = cc.symmetry(cc.certify(cc.read_edges(ROOT / 'examples' / 'clusters' / 'kagome_36.txt'), 'kagome'))
a = np.array(cc.KAGOME.vectors)
b = 2 * np.pi * np.linalg.inv(a).T          # reciprocal primitive vectors, rows


def shortest(y):
    """Cartesian k + G closest to Γ (the copy inside the first Brillouin zone)."""
    k = np.array([float(v) for v in y]) @ b
    return min((k + np.array(m) @ b for m in product(range(-2, 3), repeat=2)), key=np.linalg.norm)


def brillouin_zone(ax):
    corner = np.linalg.norm(shortest((1 / 3, 2 / 3)))
    angles = np.pi / 6 + np.arange(7) * np.pi / 3 + np.pi / 6
    hexagon = corner * np.array([np.cos(angles), np.sin(angles)]).T
    ax.fill(*hexagon.T, color='#f3f6fb', zorder=0)
    ax.plot(*hexagon.T, color=GREY, lw=1.2, zorder=1)
    ax.set_aspect('equal')
    ax.axis('off')
    return corner


def save(fig, name):
    fig.savefig(ROOT / 'docs' / name, dpi=DPI, bbox_inches='tight', facecolor='white')
    print(f'wrote docs/{name}')


# 1. The Brillouin zone and the 12 allowed momenta, coloured by star
fig, ax = plt.subplots(figsize=(4.6, 4.2))
corner = brillouin_zone(ax)
for colour, star in zip(STAR_COLOURS, sym.stars):
    name = next(s.momentum_name for s in sym.sectors if s.momentum == star[0])
    points = np.array([shortest(y) for y in star])
    ax.scatter(*points.T, s=70, color=colour, edgecolor='white', lw=1.5, zorder=3,
               label=f'{name}: star of {len(star)}')
for label, y in (('M', (0, 1 / 2)), ('K', (1 / 3, 2 / 3))):     # labels just outside the zone
    k = shortest(y)
    ax.text(*(k * (1 + 0.45 / np.linalg.norm(k))), label, fontsize=13, color=INK, ha='center', va='center')
ax.text(0.3, -0.4, 'Γ', fontsize=13, color=INK, ha='center', va='center')
ax.legend(loc='upper left', bbox_to_anchor=(0.98, 1), frameon=False, fontsize=9)
ax.set_title('Brillouin zone of kagome, allowed momenta of the 36-site cluster', fontsize=9, color=GREY)
save(fig, 'terminology_bz.png')

# 2. The point group C6v of the cluster, drawn on the lattice around a hexagon centre
fig, ax = plt.subplots(figsize=(4.4, 4.2))
pos = np.array(cc.KAGOME.positions)
centre = np.array([0.75, np.sqrt(3) / 4])           # a hexagon centre: fixed by all 12 point operations
sites = [np.array(c) @ a + p for c in product(range(-3, 4), repeat=2) for p in pos]
sites = [r for r in sites if np.linalg.norm(r - centre) < 1.6]
for i, r in enumerate(sites):
    for s in sites[i + 1:]:
        if abs(np.linalg.norm(r - s) - 0.5) < 1e-6:
            ax.plot(*zip(r, s), color=LIGHT, lw=1.6, zorder=1)
ax.scatter(*np.array(sites).T, s=28, color=GREY, zorder=2)
reach = 1.45
for w in {op.linear.tobytes(): op.linear for op in cc.space_group(cc.KAGOME)}.values():
    kind, _, normal = _classify(_cartesian(w, a))
    if kind != 'σ':
        continue
    line = np.array([-normal[1], normal[0]])
    corners = [r - centre for r in sites if np.linalg.norm(r - centre) < 0.6]     # the hexagon around the centre
    through_site = any(abs(line[0] * d[1] - line[1] * d[0]) < 1e-6 for d in corners)
    style = dict(color=STAR_COLOURS[0], ls='-') if through_site else dict(color=STAR_COLOURS[1], ls='--')
    ax.plot(*zip(centre - reach * line, centre + reach * line), lw=1.4, zorder=3, **style)
ax.plot([], [], color=STAR_COLOURS[0], label='σv: mirror through two corners of the hexagon (x3)')
ax.plot([], [], color=STAR_COLOURS[1], ls='--', label='σd: mirror through two bond centres (x3)')
ax.scatter(*centre, marker='h', s=220, color=INK, zorder=4, label='C6 axis: rotations by 60°, 120°, 180°')
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 0.02), frameon=False, fontsize=8.5)
ax.set_aspect('equal')
ax.axis('off')
ax.set_title('C6v: E, 2 C6, 2 C3, C2, 3 σv, 3 σd', fontsize=10, color=GREY)
save(fig, 'terminology_c6v.png')

# 3. A star and a little group: the 6 images of a generic k, and the operations that keep M
fig, axes = plt.subplots(1, 2, figsize=(8, 3.9))
ax = axes[0]
corner = brillouin_zone(ax)
generic = next(s for s in sym.stars if len(s) == 6)
points = np.array([shortest(y) for y in generic])
for p in points:
    ax.annotate('', p, (0, 0), arrowprops=dict(arrowstyle='->', color=LIGHT, lw=1))
ax.scatter(*points.T, s=60, color=STAR_COLOURS[2], edgecolor='white', lw=1.5, zorder=3)
ax.scatter(*points[0], s=150, facecolor='none', edgecolor=INK, lw=1.2, zorder=4)
ax.text(*(points[0] + [0.35, 0.1]), 'k', fontsize=12)
ax.scatter(0, 0, s=20, color=INK, zorder=3)
ax.set_title('the star of k = (1/6, 1/3): 6 momenta\nthe 12 operations move k onto these 6',
             fontsize=9, color=GREY)
ax = axes[1]
corner = brillouin_zone(ax)
m_star = next(s for s in sym.stars if len(s) == 3)
points = np.array([shortest(y) for y in m_star])
ax.scatter(*points.T, s=60, color=STAR_COLOURS[1], edgecolor='white', lw=1.5, zorder=3)
m = points[0]
along = m / np.linalg.norm(m)
across = np.array([-along[1], along[0]])
ax.plot(*zip(-1.25 * corner * along, 1.25 * corner * along), color=STAR_COLOURS[0], lw=1.3)
ax.plot(*zip(-1.1 * corner * across, 1.1 * corner * across), color=STAR_COLOURS[0], lw=1.3, ls='--')
ax.scatter(-m[0], -m[1], s=60, facecolor='white', edgecolor=STAR_COLOURS[1], lw=1.5, zorder=3)
ax.annotate('', -m * 0.93, m * 0.93, arrowprops=dict(arrowstyle='<->', color=GREY, lw=0.8,
                                                     connectionstyle='arc3,rad=0.35'))
ax.text(*(m + 0.5 * along + 0.45 * across), 'M', fontsize=12, ha='center', va='center')
ax.text(*(-m - 0.5 * along - 1.1 * across), '−M = M − G', fontsize=8.5, ha='center', va='center', color=GREY)
ax.scatter(0, 0, s=20, color=INK, zorder=3)
ax.set_title('the little co-group of M is C2v: E, the mirror along ΓM (solid),\n'
             'the mirror across it and C2 (both send M to −M, the same momentum)', fontsize=9, color=GREY)
save(fig, 'terminology_star.png')
