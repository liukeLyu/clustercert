"""Decide whether a finite graph is a periodic cluster of a given lattice.

A periodic cluster is the infinite lattice L folded by a full-rank lattice
Lambda of periods (translations): L/Lambda. For N sites and n_basis sites per unit cell,
Lambda has index k = N/n_basis, and there are finitely many such Lambda. Each
is listed exactly once by its row Hermite normal form (HNF). For each one,
L/Lambda is built from integer rules and compared with the input graph:
exact invariants first (edge count, degrees, distance shells, closed walks), then an
isomorphism search (networkx VF2++). A found isomorphism is re-checked edge
by edge, so a positive answer never relies on networkx.
"""
import warnings
from itertools import product

import networkx as nx
import numpy as np

__all__ = ['hnf', 'supercell_hnfs', 'cluster', 'certify', 'bonds', 'identify']


def hnf(periods):
    """Row Hermite normal form of a nonsingular integer matrix.

    Upper triangular with positive diagonal d_j and every entry above a pivot
    in range(d_j). Two matrices span the same period lattice in Z^dim if and only
    if they have the same HNF.
    """
    h = [[int(x) for x in row] for row in periods]
    dim = len(h)
    if any(len(row) != dim for row in h):
        raise ValueError('periods must be a square matrix')
    for j in range(dim):
        while True:     # Euclid on column j among rows j..dim-1
            rows = [i for i in range(j, dim) if h[i][j]]
            if not rows:
                raise ValueError('periods are linearly dependent')
            p = min(rows, key=lambda i: abs(h[i][j]))
            h[j], h[p] = h[p], h[j]
            done = True
            for i in range(j+1, dim):
                q = h[i][j] // h[j][j]
                h[i] = [x-q*y for x, y in zip(h[i], h[j])]
                done &= h[i][j] == 0
            if done:
                break
        if h[j][j] < 0:
            h[j] = [-x for x in h[j]]
        for i in range(j):
            q = h[i][j] // h[j][j]
            h[i] = [x-q*y for x, y in zip(h[i], h[j])]
    return tuple(map(tuple, h))


def _factorizations(index, dim):
    if dim == 1:
        yield (index,)
        return
    for d in range(1, index+1):
        if index % d == 0:
            for rest in _factorizations(index//d, dim-1):
                yield (d,)+rest


def supercell_hnfs(index, dim):
    """Every period lattice in Z^dim of the given index (supercell of ``index`` unit cells), once each, as its HNF.

    Supercells with the most balanced diagonal come first, so the first match is usually the natural one.
    """
    above = [(i, j) for i in range(dim) for j in range(i+1, dim)]
    for diag in sorted(_factorizations(index, dim), key=lambda d: (max(d), d)):     # most balanced cells first
        for values in product(*(range(diag[j]) for _, j in above)):
            h = [[0]*dim for _ in range(dim)]
            for k in range(dim):
                h[k][k] = diag[k]
            for (i, j), v in zip(above, values):
                h[i][j] = v
            yield tuple(map(tuple, h))


def _reduce(c, h):
    """Representative of unit cell c modulo the rows of h, in the box prod range(h_jj)."""
    c = list(c)
    for j in range(len(h)):
        q = c[j] // h[j][j]
        c = [x-q*y for x, y in zip(c, h[j])]
    return tuple(c)


def cluster(lattice, periods):
    """The periodic cluster L/Lambda as an nx.Graph on (unit cell, sublattice) nodes.

    ``periods``: integer rows generating Lambda in the primitive basis (any
    basis; it is brought to HNF). Returns None if the cluster is too small to
    be a simple graph, that is, if a bond closes on itself or two bonds join
    the same pair of sites.
    """
    h = hnf(periods)
    if len(h) != lattice.dim:
        raise ValueError(f'{lattice.name} needs {lattice.dim}x{lattice.dim} periods')
    cells = list(product(*(range(h[j][j]) for j in range(lattice.dim))))
    edges = set()
    for c in cells:
        for a, b, t in lattice.bonds:
            edge = frozenset(((c, a), (_reduce([x+y for x, y in zip(c, t)], h), b)))
            if len(edge) < 2 or edge in edges:
                return None
            edges.add(edge)
    g = nx.Graph(lattice=lattice.name, periods=h)      # recorded so that plot(g) can draw the supercell
    g.add_nodes_from((c, a) for c in cells for a in range(lattice.n_basis))
    g.add_edges_from(tuple(e) for e in edges)
    return g


def _invariant(g):
    """An exact isomorphism invariant that separates most periodic clusters.

    The multiset over sites of (distance-shell sizes, closed walks of length
    3..6). Colour refinement alone cannot tell regular graphs apart, and
    periodic clusters are regular.
    """
    nodes = list(g)
    a = nx.to_numpy_array(g, nodelist=nodes, dtype=np.int64)
    walks, power = [], a
    for _ in range(5):
        power = power @ a
        walks.append(np.diag(power))
    shells = [np.bincount(list(nx.single_source_shortest_path_length(g, v).values())) for v in nodes]
    with warnings.catch_warnings():     # networkx notes its hash values changed in v3.5; we only compare within a run
        warnings.simplefilter('ignore', UserWarning)
        wl = nx.weisfeiler_lehman_graph_hash(g, iterations=3)
    return wl, sorted((tuple(shell), *(int(w[i]) for w in walks[1:])) for i, shell in enumerate(shells))


def _compactness(periods, lattice):
    """Sort key: longest, then total, length of a reduced period basis (0 without geometry)."""
    if lattice.vectors is None:
        return (0, 0)
    from .supercell import period_vectors
    lengths = sorted(round(float((v @ v)**0.5), 9) for v in period_vectors(periods, lattice)[1])
    return (lengths[-1], sum(lengths))


def _as_graph(edges, n_sites):
    if isinstance(edges, nx.Graph):
        g = nx.Graph(edges)
    else:
        g = nx.Graph()
        for edge in edges:
            u, v = edge
            if u == v:
                raise ValueError(f'self-bond {edge} in the edge list')
            g.add_edge(u, v)
    if n_sites is not None:
        g.add_nodes_from(range(n_sites))
    return g


def _screen(g, lattice, connected=True):
    """Necessary conditions with a readable reason, or None if all pass."""
    n = g.number_of_nodes()
    if n == 0:
        return 'empty graph'
    if n % lattice.n_basis:
        return f'N={n} is not a multiple of {lattice.n_basis} (sites per {lattice.name} unit cell)'
    k = n // lattice.n_basis
    want = sorted(lattice.degrees()*k)
    have = sorted(d for _, d in g.degree)
    if have != want:
        bulk = sorted(set(lattice.degrees()))
        odd = sorted({d for d in have} - set(bulk))
        hint = ' (sites with too few bonds: open boundary or a missing wrap-around bond?)' \
            if odd and min(odd) < min(bulk) else ''
        return f'degrees {sorted(set(have))}, but every {lattice.name} site has {bulk}{hint}'
    if connected and not nx.is_connected(g):
        return f'graph has {nx.number_connected_components(g)} disconnected pieces'
    return None


def certify(edges, lattice, *, n_sites=None, find_all=False, periods=None, progress=None):
    """Is the graph a periodic cluster of ``lattice``, i.e. the lattice folded by some translations?

    ``lattice``: a :class:`Lattice`, or the name of a built-in one ('kagome').
    ``edges``: pairs of site labels (any hashable), or an nx.Graph. For
    several kinds of bond, a dict ``{family: edges}``; see
    :func:`_certify_families`.
    ``n_sites``: also include sites 0..n_sites-1 if some carry no bond.
    ``find_all``: keep searching after the first match, to list every
    matching supercell; the witness kept is then the most compact supercell.
    ``periods``: test only this supercell (integer rows in the primitive basis).
    ``progress(done, total, matches)``: optional callback.

    Returns a dict with ``certified`` (bool) and ``reason``; when certified,
    also ``periods`` (HNF rows of Lambda, primitive basis) and ``sites``
    (site label -> (unit cell, sublattice)) such that the edge set is exactly the image
    of the lattice's bonds. ``matching`` lists the HNFs found.
    """
    lattice = _lookup(lattice)
    if isinstance(edges, dict):
        return _certify_families(edges, lattice, n_sites, periods, progress)
    g = _as_graph(edges, n_sites)
    n = g.number_of_nodes()
    result = dict(lattice=lattice.name, n_sites=n, n_edges=g.number_of_edges(), certified=False,
                  index=None, tested=0, too_small=0, matching=[])
    reason = _screen(g, lattice)
    if reason:
        result['reason'] = reason
        return result
    index = n // lattice.n_basis
    result['index'] = index
    target_invariant = _invariant(g)
    target_degrees = sorted(d for _, d in g.degree)
    target_edges = {frozenset(e) for e in g.edges}
    candidates = list(supercell_hnfs(index, lattice.dim)) if periods is None else [hnf(periods)]
    if len(candidates[0]) != lattice.dim or np.prod([r[i] for i, r in enumerate(candidates[0])]) != index:
        raise ValueError(f'periods must span a supercell of {index} unit cells of {lattice.name}')
    best = None
    for count, h in enumerate(candidates, 1):
        q = cluster(lattice, h)
        result['tested'] = count
        if q is None:
            result['too_small'] += 1
        elif (q.number_of_edges() == g.number_of_edges()
              and sorted(d for _, d in q.degree) == target_degrees
              and _invariant(q) == target_invariant):
            mapping = nx.vf2pp_isomorphism(q, g)
            if mapping is not None:
                if (len(set(mapping.values())) != n
                        or {frozenset((mapping[u], mapping[v])) for u, v in q.edges} != target_edges):
                    raise RuntimeError('isomorphism failed the explicit edge check')
                result['matching'].append(h)
                key = _compactness(h, lattice)
                if not result['certified'] or key < best:
                    best = key
                    result.update(certified=True, periods=h,
                                  sites={site: node for node, site in mapping.items()})
        if progress:
            progress(count, len(candidates), len(result['matching']))
        if result['certified'] and not find_all:
            break
    if result['certified']:
        result['reason'] = f'periodic {lattice.name} cluster, periods {[list(r) for r in result["periods"]]}'
    else:
        result['reason'] = (f'not a periodic {lattice.name} cluster with periods {[list(r) for r in candidates[0]]}'
                            if periods is not None else
                            f'not a periodic {lattice.name} cluster: none of the {len(candidates)} '
                            f'possible supercells of {index} unit cells matches')
    return result


def _lookup(lattice):
    """A Lattice, or a dict of them, given by object or by built-in name ('kagome', 'pyrochlore', ...)."""
    if isinstance(lattice, dict):
        return {k: _lookup(v) for k, v in lattice.items()}
    if isinstance(lattice, str):
        from .lattice import LATTICES
        if lattice.lower() not in LATTICES:
            raise ValueError(f'unknown lattice {lattice!r}; built in: {", ".join(LATTICES)}')
        return LATTICES[lattice.lower()]
    return lattice


def _family_lattices(edges, lattice):
    if isinstance(lattice, dict):
        if set(lattice) != set(edges):
            raise ValueError(f'bond families {sorted(edges)} and lattices {sorted(lattice)} differ')
        return {k: lattice[k] for k in edges}
    from .lattice import bond_families
    shell = {k: int(str(k).lstrip('J').rstrip('abcdefgh') or 0) for k in edges}
    if not all(s >= 1 for s in shell.values()):
        raise ValueError("with a single lattice, bond families are shells 1, 2, 3, ... or split ones like '3a'")
    known = bond_families(lattice, max(shell.values()))
    out = {}
    for k in edges:
        key = str(k).lstrip('J')
        key = int(key) if key.isdigit() else key
        if key not in known:
            split = [f for f in known if str(f).startswith(str(key))]
            raise ValueError(f'shell {key} of {lattice.name} holds inequivalent bonds; give them separately as '
                             f'{split}' if split else f'no bond family {k!r} for {lattice.name}')
        out[k] = known[key]
    return out


def _certify_families(edges, lattice, n_sites, periods, progress):
    """Certify several kinds of bond at once, for example ``{1: j1_edges, 2: j2_edges}``.

    ``lattice`` is either one lattice, whose k-th neighbour shell is then the
    lattice of family k, or a dict ``{family: Lattice}`` on the same unit cell.
    The first family must be connected (normally the nearest neighbours). The
    supercells that reproduce it are found first; each is then tested with every
    family, and the result lists what goes wrong in each supercell
    (``result['supercells']``). Certified means that one supercell reproduces every
    family at once, with ``sites`` mapping all of them exactly.
    """
    lattices = _family_lattices(edges, lattice)
    keys = list(edges)
    first = keys[0]
    graphs = {k: _as_graph(edges[k], n_sites) for k in keys}
    nodes = set().union(*graphs.values())
    for g in graphs.values():
        g.add_nodes_from(nodes)
    n = len(nodes)
    name = lattices[first].name
    result = dict(lattice=name, families={k: m.name for k, m in lattices.items()}, n_sites=n,
                  n_edges={k: g.number_of_edges() for k, g in graphs.items()}, certified=False,
                  index=None, matching=[], supercells=[])
    seen = {}
    for k in keys:
        for e in map(frozenset, graphs[k].edges):
            if e in seen:
                result['reason'] = f'bond {tuple(e)} is listed in both family {seen[e]} and family {k}'
                return result
            seen[e] = k
    for k in keys:
        reason = _screen(graphs[k], lattices[k], connected=k == first)
        if reason:
            result['reason'] = f'family {k}: {reason}'
            return result
    base = certify(graphs[first], lattices[first], find_all=True, periods=periods, progress=progress)
    result.update(index=base['index'], tested=base['tested'], too_small=base['too_small'])
    if not base['certified']:
        result['reason'] = f'family {first}: {base["reason"]}'
        return result
    target = nx.Graph()
    target.add_nodes_from(nodes)
    for k in keys:
        target.add_edges_from(graphs[k].edges, kind=k)
    # networkx has no "edge colour" argument: each bond carries its family as the attribute 'kind', and
    # GraphMatcher's edge_match decides which bonds may be paired, here only bonds of the same family
    same_kind = lambda a, b: a['kind'] == b['kind']  # noqa: E731
    best = None
    for h in base['matching']:
        folded = {k: cluster(lattices[k], h) for k in keys}
        status = next((f'family {k} bonds close on themselves or double up in this supercell'
                       for k in keys if folded[k] is None), None)
        if status is None:
            owner = {}
            for k in keys:
                for e in map(frozenset, folded[k].edges):
                    if e in owner:
                        status = f'family {k} bonds land on family {owner[e]} bonds in this supercell'
                        break
                    owner[e] = k
                if status:
                    break
        if status is None:
            model = nx.Graph()
            for k in keys:
                model.add_nodes_from(folded[k])
                model.add_edges_from(folded[k].edges, kind=k)
            matcher = nx.algorithms.isomorphism.GraphMatcher(model, target, edge_match=same_kind)
            if matcher.is_isomorphic():
                mapping = matcher.mapping
                for k in keys:
                    image = {frozenset((mapping[u], mapping[v])) for u, v in folded[k].edges}
                    if image != {frozenset(e) for e in graphs[k].edges}:
                        raise RuntimeError('isomorphism failed the explicit edge check')
                status = 'matches'
                result['matching'].append(h)
                key = _compactness(h, lattices[first])
                if best is None or key < best:
                    best = key
                    result.update(certified=True, periods=h, sites={site: node for node, site in mapping.items()})
            else:
                status = f'family {first} matches, but the other families are bonded differently'
        result['supercells'].append(dict(periods=h, status=status))
    families = ', '.join(map(str, keys))
    if result['certified']:
        result['reason'] = (f'periodic {name} cluster with bond families {families}, '
                            f'periods {[list(r) for r in result["periods"]]}')
    else:
        others = ', '.join(map(str, keys[1:]))
        result['reason'] = (f'family {first} forms a periodic {name} cluster ({len(base["matching"])} matching '
                            f'supercells), but in no supercell do the bonds of family {others} match as listed')
    return result


def bonds(result, lattice):
    """The bonds of ``lattice`` in a certified cluster, in the cluster's own site labels.

    Use it to generate further-neighbour bonds that are consistent with a
    certified nearest-neighbour cluster, for example
    ``bonds(result, neighbour_shell(KAGOME, 2))``. Raises ValueError if these
    bonds close on themselves or double up in the cluster's supercell.
    """
    folded = cluster(lattice, result['periods'])
    if folded is None:
        usable = [[list(r) for r in h] for h in result.get('matching', []) if cluster(lattice, h) is not None]
        hint = (f'; the same graph also comes from supercells where they do not, with periods {usable}: '
                'pick one with certify(..., periods=...)' if usable else '')
        raise ValueError(f'{lattice.name} bonds close on themselves or double up in this supercell{hint}')
    label = {node: site for site, node in result['sites'].items()}
    pairs = [tuple(sorted((label[u], label[v]), key=str)) for u, v in folded.edges]
    return sorted(pairs, key=lambda e: tuple(map(str, e)))


def identify(edges, lattices=None, *, n_sites=None):
    """Names of the lattices (default: every built-in) the graph is a periodic cluster of."""
    from .lattice import LATTICES
    lattices = LATTICES.values() if lattices is None else lattices
    g = _as_graph(edges, n_sites)
    return [m.name for m in lattices if certify(g, m)['certified']]
