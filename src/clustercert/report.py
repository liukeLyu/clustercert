"""A readable certificate for the result of :func:`certify`."""
from .lattice import LATTICES
from .supercell import momenta, period_vectors

__all__ = ['summary']


def _fmt(row):
    return '(' + ', '.join(f'{x:.4g}' if abs(x) > 1e-12 else '0' for x in row) + ')'


def summary(result, lattice=None, max_sites=12):
    """The result of :func:`certify` as text: the verdict, and for a yes the certificate.

    The certificate is the periods and the site map (your label -> unit cell,
    sublattice), of which the first ``max_sites`` entries are shown.
    ``lattice`` is needed only for a lattice that is not built in.
    """
    lattice = lattice or LATTICES.get(result['lattice'])
    lines = [('YES: ' if result['certified'] else 'NO: ') + result['reason']]
    supercells = result.get('supercells', [])
    for cell in supercells if len(supercells) > 1 or not result['certified'] else []:
        lines.append(f'  supercell {[list(r) for r in cell["periods"]]}: {cell["status"]}')
    if not result['certified']:
        return '\n'.join(lines)
    if lattice is not None and lattice.vectors is not None:
        rows, cart = period_vectors(result['periods'], lattice)
        lines.append('  periods (in primitive vectors -> Cartesian):')
        lines += [f'    {tuple(int(x) for x in r)} -> {_fmt(v)}' for r, v in zip(rows, cart)]
    else:
        lines.append(f'  periods (in primitive vectors): {[list(r) for r in result["periods"]]}')
    if len(result['matching']) > 1:
        lines.append(f'  {len(result["matching"])} supercells give this same graph; the periods shown are one of them')
    try:
        sites = sorted(result['sites'].items())
    except TypeError:       # labels of mixed types
        sites = sorted(result['sites'].items(), key=lambda item: str(item[0]))
    if max_sites <= 0:
        lines.append(f"  site map: {len(sites)} sites, in result['sites']")
    else:
        lines.append('  site map, your label -> (unit cell, sublattice):')
        lines += [f'    {s} -> {c}, {a}' for s, (c, a) in sites[:max_sites]]
        if len(sites) > max_sites:
            lines.append(f'    ... ({len(sites)} sites in all)')
    k = momenta(result['periods'])
    points = ', '.join('(' + ', '.join(map(str, y)) + ')' for y in k)
    lines.append(f'  {len(k)} allowed momenta (in reciprocal primitive vectors): {points}')
    return '\n'.join(lines)
