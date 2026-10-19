"""Figure for the README: one kagome unit cell, and the 12-site cluster drawn in its unit cells.

Both panels are made by the library itself (cc.plot_lattice and cc.plot(..., details=True)).
Run from the repository root: python scripts/plot_definitions.py -> docs/definitions.png (a second or two).
"""
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

import clustercert as cc  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
edges = cc.read_edges(ROOT / 'examples' / 'clusters' / 'kagome_12.txt')
fig, (left, right) = plt.subplots(1, 2, figsize=(15, 5.2), gridspec_kw=dict(wspace=0.55))
cc.plot_lattice(cc.KAGOME, ax=left, title='(a) the lattice: one kagome unit cell')
cc.plot(cc.certify(edges, 'kagome', periods=[[2, 0], [0, 2]]), details=True, ax=right,
        title='(b) periods (2, 0), (0, 2): a supercell of 4 unit cells, 12 sites')
fig.savefig(ROOT / 'docs' / 'definitions.png', dpi=140, bbox_inches='tight', facecolor='white')
print('wrote docs/definitions.png')
