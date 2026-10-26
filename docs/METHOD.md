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

## 9. Symmetry sectors

`symmetry(result)` finds the symmetries of a certified cluster and the sectors they define. The terms (star, little group, irrep, Mulliken names) are explained from scratch in [`TERMINOLOGY.md`](TERMINOLOGY.md). Every step below uses only the lattice definition (primitive vectors, sublattice positions, bonds) and the certificate (periods, site map). Nothing is tabulated per lattice.

**Step 1: the space group of the infinite lattice.** Write a point in fractional coordinates f, so that its Cartesian position is x = Aᵀf, with the primitive vectors as the rows of A. A symmetry is f → Wf + t.
- W maps lattice translations to lattice translations, so it is an integer matrix. It keeps lengths, so WᵀGW = G with the Gram matrix G = AAᵀ. Column j of W is then an integer vector of the same length as a_j. There are finitely many such vectors, hence finitely many W: the *holohedry* of the lattice of translations.
- For a given W, the image of sublattice 0 must be some sublattice b, which fixes t = f_b − Wf_0 up to a lattice vector. That leaves n_basis candidates. A candidate is kept if every sublattice a lands on a sublattice σ(a), Wf_a + t = f_σ(a) + s_a with s_a integer, and every bond (a, b, τ) lands on a bond (σ(a), σ(b), Wτ + s_b − s_a).

The operation sends site (c, a) to (Wc + s_a, σ(a)), in integers. Nothing in the search assumes that t = 0, so glides and screw axes are found like any other operation. The tests check the orders: 8 for square and checkerboard, 12 for triangular, honeycomb and kagome, 48 for cubic, bcc, fcc, diamond and pyrochlore. They also check that diamond and pyrochlore are non-symmorphic: no point is fixed by every operation.

**Step 2: the symmetries of the cluster.** An operation is a symmetry of the folded lattice exactly when WΛ = Λ, that is, when W maps every period row into the period lattice. This is tested exactly, by folding W r back into the supercell (§2) and checking that it lands on cell 0. Each kept operation, followed by each translation of the supercell, becomes a permutation of your site labels through the site map. Every permutation is checked to map your bonds, family by family, onto themselves. Operations that act identically on the cluster are merged, so the group is the set of distinct permutations; its size is the number of translations (unit cells) times the number of kept point operations.

**Step 3: momenta, stars and little groups.** A translation by the cell vector c acts on a state of momentum k as e^{−ik·c}, with k·c = 2π y·c for k in reciprocal coordinates y. If g has linear part W, then g T_c g⁻¹ = T_{Wc}, so g sends momentum y to yW⁻¹. The stars are the orbits of the allowed momenta (§7) under the kept W, computed in exact rational arithmetic. The little co-group of k is the set of W with yW⁻¹ − y an integer vector; the little group G_k adds the translations. Orbit counting gives |star| × |little co-group| = |point group| for every star, which the tests check.

**Step 4: the irreps of the little group.** The sectors at k are the irreps of G_k in which every translation acts as the scalar e^{−ik·c}. They are found without assuming anything about G_k:
- Choose one element r_j of G_k for each point operation (one per coset of the translations), m of them in all, m = |little co-group|. The *induced representation* has basis e_1, …, e_m and D(g) e_i = e^{−ik·c} e_j whenever g r_i = r_j T_c.
- By Frobenius reciprocity, every irrep of G_k on which the translations act as e^{−ik·c} occurs in D exactly as often as its dimension d. So Σ d² = m, the projective analogue of the familiar rule Σ d² = |G|, and the tests check it.
- A random Hermitian matrix X, averaged into C = Σ_j D(r_j) X D(r_j)†, commutes with D. Its eigenspaces are irreducible, one per copy of each irrep, with probability one. Their characters, χ(r_j) = tr(V† D(r_j) V), identify the irreps; repeated copies are dropped. The result is accepted only if the characters are orthonormal and Σ d² = m; otherwise another X is drawn. This is the idea of Dixon's method (J. D. Dixon, *Computing irreducible representations of groups*, Math. Comp. 24, 707 (1970)), applied to a representation of dimension at most 48.
- On the rest of G_k, χ(r_j T_c) = e^{−ik·c} χ(r_j).

Nothing here depends on G_k being a direct product of translations and point operations, so non-symmorphic little groups, whose irreps are projective representations of the co-group, are handled like any other. At the X point of diamond, for example, the four irreps all come out 2-dimensional.

**Step 5: the size of each sector.** The projector onto sector (k, ρ) is P = (d/|G_k|) Σ_{g∈G_k} χ_ρ(g)* U(g). Its trace is the dimension of the sector, and divided by d it is the number of times ρ occurs: the block size for one row of ρ. Only tr U(g) is needed. U(g) permutes basis states, so its trace counts the basis states it fixes. A spin-½ configuration is fixed exactly when it is constant on every cycle of g, which gives 2^(number of cycles); with n_up up spins fixed, it is the coefficient of x^(n_up) in Π_cycles (1 + x^length). With the spin flip F, a configuration is fixed by g F when it alternates along every cycle, so tr U(g)F = 2^(number of cycles) if every cycle is even and 0 otherwise. The sizes are summed in floating point and rounded, which is exact up to about 10¹³ states per sector; the tests check that the sizes, with their degeneracies, add up to the dimension of the space exactly.

**Validation by brute force.** For one small cluster of every built-in lattice, including tilted supercells and the non-symmorphic diamond and pyrochlore, the tests build every projector, for every momentum of every star, on the full Hilbert space or a fixed-S^z part of it. They check that each one is Hermitian and idempotent with the predicted trace, and that together they sum to the identity; projectors that sum to the identity are mutually orthogonal. An XXZ Hamiltonian with random couplings on the bonds commutes with every permutation, and its eigenvalues, computed block by block in the sectors and repeated over each star, reproduce the eigenvalues of the full Hamiltonian. This is the symmetry-adapted basis used in exact diagonalization ([Wietek and Läuchli, arXiv:1804.05028](https://arxiv.org/abs/1804.05028); [Phys. Rev. E 98, 033309 (2018)](https://doi.org/10.1103/PhysRevE.98.033309)), built here as explicit matrices.

**Accidental symmetries.** The bond graph can have automorphisms that are not lattice symmetries. Their number is found by orbit and stabilizer: pin a site, count its orbit under the automorphisms that fix the sites pinned before, and repeat until colour refinement separates every site, at which point only the identity fixes the pinned sites. Each orbit contains the orbit under the lattice symmetries, so only one site per remaining lattice orbit needs an isomorphism search (VF2++).

**Names.** The names are computed from the operations, after the sectors are found, and are never used to find them.
- *Point groups* (Schoenflies) follow from the kinds of operations present: the order of the highest rotation axis, perpendicular C2 axes, horizontal and vertical mirrors, inversion, S_2n, and several C3 axes for the cubic groups.
- *Irreps* (Mulliken) follow from the characters: A or B by the main rotation (or S4 for D2d), E and T by dimension, 1 or 2 by the C2 axis perpendicular to the main axis, or the vertical mirror, closest to the x axis, g or u by inversion, ′ or ″ by the horizontal mirror. The characters must be those of a genuine representation of the co-group. At Γ they are. At other k, the rotations are taken about the most symmetric point fixed by every operation of the little co-group (searched on a grid of fractional coordinates with denominator 24); if there is no such point, or the irrep is complex, it is numbered instead.
- *Momenta* are named by the type of Bravais lattice, read off the order of the holohedry and the number of shortest translations, and by the distance of k from Γ and the order of its stabilizer in the holohedry. This does not depend on the choice of primitive vectors or on the orientation.

**Scope.** Spatial symmetries of the lattice, and the spin flip, are included; time reversal, spin rotation and spin–orbit coupling (double groups) are not. A model whose couplings depend on the sublattice, or whose bonds carry an orientation, can have less symmetry than its lattice. The unit cell must be primitive: an operation that does not map the chosen Bravais lattice onto itself is not found. With a non-primitive cell, several translations t belong to the same W; the analysis still works, but point groups and irreps are then numbered rather than named.
