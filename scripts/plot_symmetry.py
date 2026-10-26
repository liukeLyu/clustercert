"""Figures for the README symmetry pictures, all drawn by clustercert itself.

Run from the repository root: python scripts/plot_symmetry.py -> docs/symmetry_*.png (a few seconds).
"""
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

import clustercert as cc  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CLUSTERS = ROOT / 'examples' / 'clusters'


def save(name):
    plt.savefig(ROOT / 'docs' / name, dpi=130, bbox_inches='tight', facecolor='white')
    plt.close('all')
    print(f'wrote docs/{name}')


cc.symmetry(cc.certify(cc.read_edges(CLUSTERS / 'kagome_12.txt'), 'kagome')).plot(n_up=6)
save('symmetry_quickstart.png')
cc.symmetry(cc.certify(cc.read_edges(CLUSTERS / 'kagome_36.txt'), 'kagome')).plot(n_up=18, mosaic=True)
save('symmetry_mosaic.png')
tilted = [[3, 1], [-1, 3]]
cc.plot_symmetry_elements(cc.symmetry(cc.certify(list(cc.cluster(cc.SQUARE, tilted).edges), 'square',
                                                 periods=tilted)))
save('symmetry_elements.png')
