"""Run the Python snippets in README.md so the documented usage stays correct."""
import re
from pathlib import Path

import matplotlib

matplotlib.use('Agg')     # the quick start draws a figure

ROOT = Path(__file__).resolve().parents[1]
BLOCKS = re.findall(r"```python\n(.*?)```", (ROOT / 'README.md').read_text(), flags=re.S)
# which of the equally valid isomorphisms networkx returns depends on its version, so site-map lines are not compared
SITE_LINE = re.compile(r"^ +\d+ -> \(")


def _stable(lines):
    return [line for line in lines if not SITE_LINE.match(line)]


def test_readme_snippets_run(monkeypatch):
    monkeypatch.chdir(ROOT)
    namespace: dict = {}
    for block in BLOCKS:
        exec(compile(block, 'README.md', 'exec'), namespace)
    assert namespace['result']['certified'] and len(namespace['r16']['matching']) == 4
    assert not namespace['cc'].certify(namespace['impostor'], 'kagome')['certified']
    assert sorted(map(tuple, abs(namespace['cartesian']).round(9))) == [(0, 0, 1), (0, 1, 0), (1, 0, 0)]
    assert namespace['lieb'].degrees() == [4, 2, 2]


def test_readme_shown_output_is_real(monkeypatch, capsys):
    """Each python block followed by a text block must print (at least) the lines shown, '...' aside."""
    import subprocess
    monkeypatch.chdir(ROOT)
    text = (ROOT / 'README.md').read_text()
    namespace: dict = {}
    for code, shown in re.findall(r"```python\n(.*?)```\n\n```text\n(.*?)```", text, flags=re.S):
        exec(compile(code, 'README.md', 'exec'), namespace)
        printed = capsys.readouterr().out.splitlines()
        assert all(line in printed for line in _stable(shown.splitlines()) if line.strip() != '...'), shown
    for command, shown in re.findall(r"```text\n\$ (clustercert .*?)\n(.*?)```", text, flags=re.S):
        out = subprocess.run([str(ROOT / '.venv/bin/clustercert') if (ROOT / '.venv').exists() else 'clustercert',
                              *command.split()[1:]], capture_output=True, text=True).stdout
        assert _stable(out.splitlines()) == _stable(shown.splitlines()), command
