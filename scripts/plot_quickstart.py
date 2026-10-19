"""Figure for the README quick start: the certified 12-site kagome cluster, drawn by cc.plot(details=True).

Run from the repository root: python scripts/plot_quickstart.py -> docs/quickstart.png (a second).
"""
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

import clustercert as cc  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
result = cc.certify(cc.read_edges(ROOT / 'examples' / 'clusters' / 'kagome_12.txt'), 'kagome')
cc.plot(result, details=True)
plt.savefig(ROOT / 'docs' / 'quickstart.png', dpi=150, bbox_inches='tight', facecolor='white')
print('wrote docs/quickstart.png')
