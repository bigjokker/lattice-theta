# Six balanced minima force a symmetry: a finite classification

2026-10-23, Claude Opus 5.5. This implements priority 4 of
`RANK4-NIGHT-REVIEW-20261022.md`, the finite marked-cell classifier specified
in section 3 of `RANK4-SIX-MINIMA-COMMON-LATTICE.md`.

**Status: accepted by Astra's final audit as a computer-assisted proof, 2026-10-23.**
See [RANK4-THEOREM-S-FINAL-AUDIT-ASTRA.md](RANK4-THEOREM-S-FINAL-AUDIT-ASTRA.md)
for the checked scope, proof chain, source corrections and reproducible artifacts.
The audit checked all 502,740 positive and negative dispositions, recomputed
the parity-class minima on every positive definite face, verified every
uniform symmetry certificate, and matched both saved fan files to the
published cones. The broader Theorem R/R-prime remains a separate candidate.
Everything is exact rational arithmetic, except the
floating-point Delaunay triangulation, which only proposes candidates that are
then verified exactly. Every claimed step and every file is listed so it can
be checked.

## 1. Results

**Theorem S (computer-assisted).** Let H be any positive definite real quadratic form
on Z^4, and τ = I + N an involution of F2^4 with N^2 = 0 and rank N = 2
(210 such τ). Suppose that for each of the six pairs {k, τk} of classes moved
by τ, the classes k and τk of Z^4/2Z^4 have the same minimal H-norm and the
same number of minimal vectors. Then some T in GL4(Z) satisfies T^t H T = H
and T ≡ τ (mod 2).

**Corollary 1 (six balanced minima imply a symmetry zero).** Let L be even,
positive definite, of rank four and odd determinant, and let P be an
isotropic plane of L/2L. If all six sub-coset minimal shells of P are
balanced, then some h in Aut(L) is a symmetry witness (in the sense of
`RANK4-MINUS-TYPE.md`) for every ordered pair (a, c) of distinct nonzero
classes of P. So the three signed theta series of P vanish identically.

**Corollary 2 (the minus-type conjecture).** Every even, orthogonally
indecomposable, positive definite rank-four lattice with det = 5 mod 8 has
property P: no **nontrivial even characteristic** (a, b), meaning a and b
both nonzero mod 2 and a.b even, has an identically vanishing theta series.
More generally, for an even rank-four lattice of odd determinant, every
vanishing theta series of a nontrivial even characteristic is a symmetry
zero (level 1). In indecomposable lattices these are exactly the zeros of
Theorem 5: totally singular planes of plus type with an h with h^2 = -I.
Odd characteristics, which vanish automatically, are outside these
statements.

### Proof outline

1. **Plane invariance** (`RANK4-STAGE2-STRATEGY.md` section 2). If
   Theta[a, b] vanishes, so do all three shifts of P = span(a, c).
2. **Six series** (exact; `RANK4-SIX-MINIMA-COMMON-LATTICE.md` section 1).
   Each shift splits by norm mod 8 into two halves. Each half equals
   theta(v + 2K) - theta(w + 2K) at every coefficient, so each of the six
   halves has a balanced minimal shell. In the language of H = (1/2)B^t G B
   on K, with the marking τ of section 2 of that note, this is exactly the
   hypothesis of Theorem S.
3. **Theorem S** gives T in Aut(H) with T ≡ τ (mod 2).
4. **Extension to L** (section 2 of the common-lattice note). T ≡ τ fixes W =
   2L/2K pointwise mod 2K, so T(2L) = 2L. Then h = B T B^-1 is in Aut(L).
   - τ induces the identity on V/W = P = K/2L, so h fixes a, c and a + c
     mod 2L.
   - τ swaps the two cosets of every half, which carry opposite signs. So
     h reverses the character on every shift: ε(h) = 1.
   - Hence h is a symmetry witness.
5. **Theorem 5** (`RANK4-MINUS-TYPE.md`, indecomposable, odd det): a symmetry
   witness exists only when there is h' with h'^2 = -I, and that forces plus
   type. So an indecomposable minus-type lattice has no zero. ∎

Steps 1, 2, 4 and 5 are existing, reviewed project results (step 4 is
Astra's), so the new ingredient is Theorem S.

## 2. How Theorem S was computed

`rank4_secondary_fan.py` and `rank4_six_cells.py`; outputs in
`results/rank4-secondary-fan/` and `results/rank4-six-cells/`.

**(a) The secondary fan (Voronoi's L-type domains) of quaternary forms.**
1. For a generic rational form, the Delaunay simplices at 0 are proposed by a
   floating triangulation of a box of lattice points.
2. Each is verified **exactly**: no lattice point lies inside or on its
   circumsphere apart from its vertices, by an exact enumeration around the
   centre.
3. The star is checked complete: total |det| = 120. Each simplex is counted at
   its 5 vertices, and covol(Z^4) = 1.
4. For each facet, the affine dependence of the six points (simplex plus the
   neighbour's opposite vertex) gives the regulator inequality r . q >= 0.
   The irredundant inequalities and the extreme rays come from cdd with exact
   rationals.
5. Voronoi's algorithm: start at A4* (Gram 5I - J), cross every facet that
   meets the positive definite cone, and identify cones up to GL4(Z) by an
   exact isometry test of their central forms (sum of primitive rays).

Result: **3 orbits of maximal cones, each simplicial with 10 facets and 10
rays.** This matches the classical count of three primitive L-types in
dimension four. Central forms:

    [[4,-1,-1,-1],[-1,4,-1,-1],[-1,-1,4,-1],[-1,-1,-1,4]]   (principal, A4*)
    [[5,-2,-2,-2],[-2,6,-1,-1],[-2,-1,6,-1],[-2,-1,-1,6]]
    [[6,-2,-3,-3],[-2,6,-1,-1],[-3,-1,6,1],[-3,-1,1,6]]

**(b) Why faces suffice.** The closures of the fan's cones cover the positive
definite cone, and every cone is a GL4(Z)-image of one of the three
representatives. The cones are simplicial, so every subset of rays is a face.
On the relative interior of a face the Delaunay subdivision is constant, and
so is the set of minimal vectors of every class of Z^4/2Z^4. The reason: v is
minimal in v + 2Z^4 iff v/2 lies in the Voronoi cell of 0, iff 0 and v are
vertices of a common Delaunay cell, because Vor(0) ∩ Vor(v) is symmetric
about v/2. This constancy was also checked at two different interior points
of every face: 0 exceptions in 2,394 positive definite faces. Theorem S's
hypothesis and conclusion are GL4(Z)-equivariant, with τ conjugated mod 2.
So it suffices to treat each face of the three representatives with **all
210** τ.

**(c) The cells.** For a face S (a ray subset whose sum is positive definite)
and a τ:
1. Require equal multiplicities |M_k| = |M_τk| for the six pairs.
2. The six equalities H(v_k) = H(v_τk) are linear in the ray weights.
3. The cell is {weights > 0 on S, the six equalities}. Its closure's rays come
   from cdd exactly, and the cell meets relint(S) iff each ray of S has
   positive weight on some cell ray. Otherwise the cell lies in a smaller
   face and is treated there.

**(d) Certification.** At a generic point of each cell, the complete
automorphism group is computed and searched for a T with T ≡ τ (mod 2) that
preserves **every generator** of the cell, hence every form in it.

| | Count |
|---|---:|
| faces (ray subsets) | 3,069 |
| positive definite faces | 2,394 |
| face/marking cells with six balanced minima | 5,256 |
| ... certified | **5,256** |
| ... uncertified | **0** |

Cells by dimension: 1: 84; 2: 1,185; 3: 2,118; 4: 1,422; 5: 408; 6: 39.
The witnesses T split as T^2 = I in 4,093 cells, T^2 = -I in 302, and other
in 861. Run time 9.3 s on 12 workers.

## 3. Validation

| Check | Result |
|---|---|
| Independent recheck of all 5,256 cells (`results/rank4-six-cells/validate.py recheck`) | at every generic point the six pairs are balanced, by an independent complete box (`verify_lattice_census.box_vectors`); every T is integral, unimodular, ≡ τ mod 2, and preserves every saved generator |
| Negative control: generic forms of the open cones, all markings tried (no balance) | 0 of 120 certified |
| Negative control: balanced cells with a different marking | 0 of 300 certified |
| Multiplicity filter | with only the minimum equalities imposed, no extra cells appear: equal minima never come with unequal multiplicities on a face |
| Transvection markings (N of rank 1, 105 markings), the same pipeline | 10,446 cells, all certified. This suggests a broader rigidity: parity-class minima with multiplicities may determine a quaternary form up to the level-2 congruence subgroup. Astra's note mentions minimum-data rigidity for unmarked lattices; this marked version may be known or provable directly. |
| Recomputing the fan from another generic start | again 3 orbits from the start [[7,2,-3,1],[2,9,1,-4],[-3,1,8,2],[1,-4,2,11]]; each matches exactly one saved orbit by an exact isometry of central forms (`validate-fan2.json`) |
| Consistency with all earlier data | every six-of-six example (5,868 + 1,559 census planes, 63,132 box planes, 1,895 MILP points) has a symmetry; no indecomposable minus-type six-of-six exists through det 1997 or in the box; the det-9373 family and all near-misses are 5 of 6 |

## 4. What needs independent review

1. **The fan computation.** The three orbits match the classical count. The
   cones being simplicial is to be checked against the literature (Delone and
   Štogrin; Vallentin's thesis, chapter on L-types of dimension 4). A
   cross-check of the 52 Delaunay types, or the 51 iso-edge types, from the
   face enumeration would be stronger.
2. **The constancy argument in 2(b)**, minimal class vectors versus Delaunay
   cells: standard, but it should be stated with a citation.
3. **Coverage.** Every face of every fan cone is GL4(Z)-equivalent to a face
   of one of the three representatives, and the treatment of cells lying in
   smaller faces.
4. **The software trust base:** `verify_rank4_feasibility.box_isometries`
   (complete automorphism groups), `rank4_six_minima.class_minimum` and
   `rank4_pairing_search._enumerate` (complete class minima), and cdd's exact
   ray enumeration.
5. **Steps 4 and 5 of the proof outline:** that h meets the definition of a
   symmetry witness used in Theorem 5, and that Theorem 5's hypotheses
   (even, indecomposable, odd det) are exactly the ones needed.

If this holds, the remaining census and discovery work for minus type is
moot: no example exists at any determinant. The computation also gives a
uniform certificate format: for any even odd-det lattice with a zero, the
symmetry witness can be produced directly.

## Files

| File | Content |
|---|---|
| `rank4_secondary_fan.py` | exact Delaunay stars, regulators, secondary cones, Voronoi's algorithm |
| `results/rank4-secondary-fan/fan.json` | the three cone orbits: facets, rays, central forms |
| `rank4_six_cells.py` | the classification |
| `results/rank4-six-cells/cells.json` | all 5,256 cells: face, τ, dimension, generic point, generators, T |
| `results/rank4-six-cells/validate.py`, `validate-*.json` | the checks of section 3 |
| `rank4_six_minima.py`, `test_rank4_six_minima.py` | common lattice, six-minimum records, certificate checker |
