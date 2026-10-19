"""The tour notebook runs end to end."""
from pathlib import Path

import pytest

nbformat = pytest.importorskip('nbformat')
nbclient = pytest.importorskip('nbclient')

EXAMPLES = Path(__file__).resolve().parents[1] / 'examples'


def test_tour_runs(monkeypatch, tmp_path):
    monkeypatch.setenv('MPLBACKEND', 'Agg')
    monkeypatch.setenv('MPLCONFIGDIR', str(tmp_path))
    nb = nbformat.read(EXAMPLES / 'tour.ipynb', as_version=4)
    nbclient.NotebookClient(nb, timeout=120, resources={'metadata': {'path': str(EXAMPLES)}}).execute()
