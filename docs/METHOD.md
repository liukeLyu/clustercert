# Method

Why `certify` gives the right answer. The [README](../README.md) defines the terms: unit cells c, sublattices a, sites (c, a), bonds (a, b, t), periods and supercells.

## 1. The question

The infinite lattice has a site (c, a) for every unit cell c (an integer vector) and every sublattice a, and a bond from (c, a) to (c + t, b) for every c and every listed bond (a, b, t). Folding by periods identifies any two unit cells that differ by an integer combination of the periods. The result is a finite graph, the *folded lattice*.

The question is whether your bond list, with your own site labels, is a folded lattice: is there a choice of periods, and a relabelling of your sites, under which your bonds are exactly the folded bonds? That is a graph isomorphism between your cluster and some folded lattice.

The search uses only the integers in the lattice's bond list. Geometry (primitive vectors, sublattice positions) is used only to report results.

## 2. There are only finitely many foldings

A cluster of N sites, with n sublattices, has k = N/n unit cells, so the periods must leave exactly k distinct unit cells. Many period matrices describe the same folding: adding one period to another, or flipping a sign, changes nothing. Row reduction with integers (Euclid's algorithm on each column) brings any of them to one canonical form, the Hermite normal form: an upper-triangular matrix whose diagonal entries d₁, …, d_d are positive and multiply to k, and whose entries above each diagonal entry d_j lie in 0, …, d_j − 1. In two dimensions the periods are (d₁, x) and (0, d₂), with d₁d₂ = k and 0 ≤ x < d₂.

So every folding with k unit cells appears exactly once in a finite list:

| dimension | k = 4 | k = 6 | k = 8 |
|---|---:|---:|---:|
| 2 | 7 | 12 | 15 |
| 3 | 35 | 91 | 155 |

These are σ(k) in 2D (OEIS A000203) and OEIS A001001 in 3D; the tests check both counts.

The triangular form also labels the k unit cells of the supercell: the cells c with 0 ≤ c_j < d_j. To fold any unit cell into this box, subtract ⌊c_j / d_j⌋ times period j, for j = 1, then 2, and so on. Period j is zero before position j, so each step leaves the coordinates already fixed alone. Two cells of the box never differ by a combination of periods: in the first coordinate where they differ, that combination would have to be a multiple of d_j smaller than d_j. Folding is therefore pure integer arithmetic, with no coordinates and no tolerances.

Hence your cluster is a periodic cluster exactly when it is isomorphic to one of these finitely many folded lattices. One exception is skipped: in a very small supercell a bond can close on itself, or two bonds can join the same pair of sites, and such a folded lattice can never match a bond list. `hnf` gives the canonical form of any period matrix, so two period matrices describe the same supercell exactly when their forms agree.

## 3. Deciding isomorphism, and why to trust the answer

**A yes** comes with the relabelling: a (unit cell, sublattice) for each of your sites. Your bond set is then checked to be exactly the folded bond set. Anyone can repeat that check, so a yes does not depend on the search code.

**A no** means that no candidate matched. Most candidates are ruled out quickly, by quantities that every isomorphism preserves:
- the number of bonds and the degrees;
- for each site, the number of sites at each graph distance, and the number of closed walks of lengths 3 to 6 that start there;
- a Weisfeiler–Lehman hash.

The distance and walk counts are needed because periodic clusters are regular graphs, and the Weisfeiler–Lehman hash alone cannot tell regular graphs apart. Each surviving candidate is decided by VF2++ (A. Jüttner and P. Madarasi, Discrete Appl. Math. 242, 69 (2018)), a complete isomorphism search. A no is therefore exhaustive, but unlike a yes it relies on networkx's VF2++ being correct.

Before any of this, cheap necessary conditions give readable reasons: n divides N, the degrees are those of the lattice, and the cluster is connected.

## 4. Local checks are not enough

Every pyrochlore site has six neighbours and lies in two corner-sharing tetrahedra. The reason is that pyrochlore is the line graph of the diamond lattice: one site per diamond bond, with two sites joined when their bonds meet ([Chiu et al., arXiv:2109.13940](https://arxiv.org/abs/2109.13940); [Phys. Rev. Research 4, 023063 (2022)](https://doi.org/10.1103/PhysRevResearch.4.023063)). The four bonds that meet at a diamond vertex become a tetrahedron.

The same holds for the line graph of *any* graph in which every vertex has four neighbours and there are no triangles. Its only cliques are the stars of its vertices, so its sites also have six neighbours and two tetrahedra each. It is a pyrochlore cluster only if the starting graph is a folded diamond lattice, a global property that no local check can see. The tests include a 32-site example that passes every local check and is rejected for all 155 candidates. Kagome (the line graph of honeycomb) and checkerboard (of square) work the same way; the README's 12-site kagome impostor is the line graph of the Wagner graph.

## 5. One graph, several supercells

Usually several supercells give the same graph. Rotations and reflections of the lattice account for most of them, but geometrically different supercells can also match:
- 16-site pyrochlore: the cubic cell, and three supercells whose periods have lengths 1, 1 and √(3/2);
- 12-site kagome: the 2×2 rhombus, and supercells whose periods have lengths √3 and 2.

A nearest-neighbour Hamiltonian is fixed by the graph, so its spectrum is the same for all of them. The translation group, and with it the momentum labels, belongs to the supercell you choose. `find_all=True` lists every matching supercell and keeps the one whose periods are shortest after LLL reduction ([Lenstra, Lenstra and Lovász, Math. Ann. 261, 515 (1982)](https://doi.org/10.1007/BF01457454)).

## 6. Several kinds of bond: J₁, J₂, J₃, …

**The question.** Colour every bond by its family: J₁, J₂, J₃, …. Your cluster is then an edge-coloured graph, and so is every folded lattice. A J₁–J₂–J₃ certificate is one relabelling of your sites that maps bonds to bonds *of the same colour*: an isomorphism of edge-coloured graphs. It must be one relabelling for all colours, because the Hamiltonian has one set of sites. Matching each family on its own would allow a different site map for each term.

**How it is found, on an example.** Take the 12-site kagome cluster with its J₁ and J₂ lists.

1. *J₁ alone.* The search of §2–3 finds the supercells that reproduce the J₁ bonds: four of them, with periods (2, 0), (0, 2) (the 2×2 rhombus), (2, 1), (0, 2), (1, 1), (0, 4) and (1, 2), (0, 4). A J₁–J₂ relabelling is in particular a J₁ relabelling, once you ignore the J₂ colour, so the answer must be one of these four. No other supercell needs to be tried.
2. *Fold J₂ in each.* In three of them, two J₂ bonds land on the same pair of sites (see below), so they cannot host J₂. The rhombus remains.
3. *Match both colours at once.* In the rhombus there are 48 relabellings of your sites that map J₁ onto J₁. Any two of them differ by a symmetry of the folded J₁ graph, and here those are the 4 translations times the 12 rotations and reflections of a hexagon. A J₁–J₂ certificate is one of the 48 that also maps your J₂ list onto the folded J₂ bonds. With the correct list, all 48 do. With a list in which two J₂ bonds have their ends swapped, none does, although every site still has four J₂ neighbours. The verdict is then "family 1 matches, but the other families are bonded differently".

   The search does not list the 48 and then test them. It builds one relabelling site by site (networkx VF2, comparing edge colours). Each time it assigns one of your sites, it checks the bonds to the sites already assigned: each bond must be present on both sides or on neither, with the same colour. When a check fails, it undoes the last assignment and tries the next one. A relabelling that is right for J₁ but wrong for J₂ is therefore dropped at its first wrong J₂ bond. The search stops at the first relabelling that passes every check, and finding none means every possibility was ruled out. (The count of 48 comes from listing them separately, with `networkx.vf2pp_all_isomorphisms`, to illustrate the point.)

   *Colour in networkx.* networkx has no "edge colour" argument. Instead, every bond carries an attribute `kind`, its family key (`1`, `2`, `'3b'`, …). [`GraphMatcher`](https://networkx.org/documentation/stable/reference/algorithms/isomorphism.vf2.html) is passed a function `edge_match` that says when two bonds may be paired: here, when their `kind` is equal. VF2 calls this function for every bond between a newly assigned site and the sites already assigned, which is the colour check above. In the code this is `GraphMatcher(model, target, edge_match=same_kind)` in `src/clustercert/certify.py`, with `same_kind = lambda a, b: a['kind'] == b['kind']`. networkx's ready-made `categorical_edge_match('kind', None)` does the same.

Step 3 has to be a search: the relabelling found from J₁ alone may not work for J₂. On the 4×4 square torus, the J₁ graph is the four-dimensional hypercube, with 384 relabellings that preserve J₁. Only 128 of them, the symmetries of the square lattice, also preserve the diagonal J₂ bonds. A correct J₁–J₂ cluster whose sites were numbered by one of the other 256 is still a genuine square-lattice cluster, but the J₁-only site map sends its J₂ list to the wrong pairs. The colour-matching search finds a relabelling that works for both colours and certifies it. As before, a yes is re-checked bond by bond, family by family, and a no is exhaustive.

**When a supercell cannot host a family.** Folding can make two J₂ bonds join the same pair of sites, which is then coupled by 2J₂. It can also put a J₂ bond on top of a J₁ bond, giving J₁ + J₂. Either way the cluster no longer carries the lattice's model, whatever bond list you give. `certify` reports this for each supercell, separately from "the bonds are wired differently". Often that report is the useful answer, since a J₁ certificate says nothing about J₂:

| Cluster | supercells matching J₁ | … with clean J₂ | … with clean J₂ and across-hexagon J₃ |
|---|---:|---:|---:|
| kagome 12 | 4 | 1 (2×2 rhombus) | 0 |
| kagome 27 | 1 | 1 | 1 |
| kagome 36 | 1 | 1 | 1 |
| pyrochlore 16 | 4 | 0 | – |
| pyrochlore 32 | 10 | 1 (2×2×2 FCC) | – |

**Where the families come from.** `neighbour_shell(lattice, k)` bonds every pair of sites at the k-th smallest distance. It searches enough periodic images to miss none: two sites whose unit cells differ by more than R in some coordinate are at least σ·R − 2p apart, where σ is the smallest singular value of the matrix of primitive vectors and p is the largest distance of a sublattice from its cell's origin.

A distance shell can mix bond types. On kagome, two sites two bond lengths apart lie either along a straight line or across a hexagon, and the BFG model couples only the second. `bond_families` splits each shell by two numbers: how many J₁ steps separate the two sites, and along how many shortest J₁ paths. A lattice symmetry maps paths to paths of the same length, so equivalent bonds always agree on both numbers, and the split never breaks a bond type. It is not guaranteed to separate every pair of inequivalent types; if two types agree on both numbers, give one `Lattice` per type instead. Up to the third shell:

| Lattice | third shell |
|---|---|
| kagome | `3a`: 2 J₁ steps (along a line), 4 per site; `3b`: 3 J₁ steps via 2 paths (across a hexagon), 2 per site |
| pyrochlore | `3a`: 2 J₁ steps, 6 per site; `3b`: 3 J₁ steps via 2 paths (across a hexagon), 6 per site |
| square, triangular, honeycomb | one family |

The tests check the kagome split independently: the third-neighbour pairs of a 6×6 cluster built from Cartesian geometry, sorted by J₁ distance, certify as `3a` and `3b`, and swapping the two is caught.

## 7. Allowed momenta

A momentum k is allowed if e^{ik·T} = 1 for every period T. Write k in units of the reciprocal primitive vectors, k = y₁b₁ + … + y_d b_d, where bᵢ·aⱼ = 2π when i = j and 0 otherwise. A period with integer row r is T = r₁a₁ + … + r_d a_d, so k·T = 2π (r·y). The condition is that r·y is an integer for every period row r, that is, H y is an integer vector, where H is the period matrix. So y = H⁻¹m for integer vectors m. Choosing m with 0 ≤ m_j < d_j gives each allowed momentum once, k of them in all, one per unit cell. `momenta` computes them in exact rational arithmetic.

## 8. Scope

Covered: translation boundary conditions, and several kinds of bond as in §6 (with a connected J₁ family). Not covered: twisted, screw or glide boundary conditions, open or slab boundaries, and supercells so small that bonds double up. Bond orientations (chiral or Dzyaloshinskii–Moriya terms) are not checked, but distinct bond types, such as Kitaev's x, y and z bonds, can be certified as separate families.
