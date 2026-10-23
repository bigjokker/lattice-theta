# Theorem S: the minus-type parity conjecture holds at every determinant

Release 1.4.0, 2026-10-23. This is a computer-assisted proof, accepted in a
final audit by a Codex session (Astra), with an independent exact replay
by a second AI system (Grok Bot) and a literature and proof audit by a third
(Grok). Claude Opus 5.5 sessions designed and ran the computation. As
elsewhere in this repository, these are AI-assisted reviews with independent
computational replication, not peer review, and no priority claim is made.

The audited files are preserved byte for byte in
[`reports/rank4-theorem-s/package/`](../reports/rank4-theorem-s/package/);
see its `README.md` and `MANIFEST.json`.

## 1. Statements

Notation is as in the [follow-up note](../papers/rank-four-follow-up-results.md),
section 6:
- L is even, positive definite, of rank four, with odd determinant;
- V = L/2L, with Q(x) = N(x)/2 mod 2 and B(x, y) = x^t G y mod 2;
- a characteristic (a, b) is **nontrivial even** if a, b are nonzero mod 2
  and a.b is even;
- c = G^-1 b mod 2.

**Theorem S.** Let H be a positive definite real quadratic form on Z^4, and
τ = I + N an involution of F2^4 with N^2 = 0 and rank N = 2. Suppose that,
for each of the six pairs {k, τk} of classes of Z^4/2Z^4 moved by τ, the two
classes have equal minimal H-norm and equal numbers of minimal vectors. Then
some T in GL4(Z) satisfies T^t H T = H and T ≡ τ (mod 2).

**Corollary S1 (the minus-type parity conjecture).** An even, orthogonally
indecomposable, positive definite rank-four lattice with determinant ≡ 5
(mod 8) has no identically vanishing theta series of a nontrivial even
characteristic, at any determinant.

**Corollary S2 (odd determinant).**
- In an even positive definite rank-four lattice of odd determinant, every
  identically vanishing theta series of a nontrivial even characteristic is
  a symmetry zero (level 1).
- With Theorem E (section 6), in an indecomposable such lattice every zero
  comes from an automorphism h with h^2 = -I and lies on the totally singular
  plane (h + I)L/2L. So such lattices have zeros only in plus type.

The scope is exactly what is stated:
- odd characteristics are excluded;
- decomposable minus-type lattices can have zeros (the det-21 control);
- the results say nothing about odd lattices or even determinants, where
  genuine level-4 zeros exist (determinants 280, 1729 and 4480);
- they do not classify isospectral quaternary pairs.

## 2. Proof

1. **Plane invariance** ([stage-two reductions](rank-four-stage-two-reductions.md);
   Theorem F). For an even rank-four Gram G of odd determinant d we have
   d ≡ 1 (mod 4). With A = adj(G) and δ = (3d + 1)/2, the integral
   symplectic matrix with blocks (2I, 3G; A, δI) and the translation z ↦ z + 1
   permute the ordered bases of the plane P = span(a, c). So the vanishing of
   one characteristic of P implies the vanishing of all three shifts of P.
   The automorphy multiplier depends on the characteristic. It is nonzero
   because, at a fixed positive definite period, the theta function in its
   elliptic variable is not identically zero, and the transformation changes
   that variable invertibly. (The case a = c cannot vanish: the sign is
   constant on each norm shell.)
2. **Six series in a common lattice.** Let K = {y in L : y mod 2L in P}
   (index 4, containing 2L), with basis columns B and H = B^t G B / 2
   (integral, det H = det G).
   - For a shift (a, c) and y = a + 2x, N(y) ≡ N(a) + 4B(a, x) (mod 8), so
     the shift splits into two halves by B(a, x) mod 2.
   - Each half equals θ(v + 2K) - θ(w + 2K) at **every** coefficient, for two
     classes v, w of K/2K.
   - The twelve classes involved form six such pairs, swapped by a linear
     involution τ of K/2K with τ - I of rank 2 and square zero.
   - A zero makes all six differences vanish. In particular each pair has
     equal minimal norm and multiplicity: Theorem S's hypothesis for (H, τ).
3. **Theorem S** gives T in Aut(H) with T ≡ τ (mod 2).
4. **Extension.** T fixes W = 2L/2K pointwise mod 2K, so T(2L) = 2L, and
   h = B T B^-1 is in Aut(L).
   - h fixes P mod 2L.
   - For distinct nonzero a, c in P, B((ha - a)/2, c) ≡ 1, so h reverses the
     character on every shift.
   - Hence h is a symmetry witness (ε(h) = 1), for every ordered basis of P.
5. **Theorem E** ([minus type](rank-four-minus-type.md), Theorem 5). In an
   even, indecomposable rank-four lattice of odd determinant, a symmetry
   witness forces an h' with h'^2 = -I, whose fixed plane (h' + I)L/2L is
   totally singular. That forces plus type (det ≡ 1 mod 8). So minus type has
   none. ∎

## 3. Proof of Theorem S (computation)

**Reduction to finitely many cells.**
1. **Cones.** Voronoi's secondary fan of positive definite quaternary forms
   has three GL4(Z)-orbits of maximal cones (L-types), all simplicial with 10
   rays.
   - The computed cones are matched, by explicit unimodular maps checked on
     all ten rays, to Voronoi's Δ, Δ′, Δ″. The rays are as printed in
     Vallentin, Weißbach and Zimmermann, *The chromatic number of
     4-dimensional lattices*, arXiv:2407.03513, section 2.1 (Indag. Math. 36
     (2025), p. 992).
   - Coverage of the positive definite cone by their faces and GL4(Z)-images
     is Voronoi's theorem; see the modern statement, Theorem 2.3 of
     arXiv:math/0601084, and Delone 1929, p. 161: "il n'existe que 3 domaines
     de Dirichlet primitifs dans l'espace à 4 dimensions".
2. **Faces.** The rays are positive semidefinite, so a face has positive
   definite relative interior iff its ray sum is positive definite. There are
   2,394 such faces, counted on each of the three representatives with
   repetitions.
3. **Minimal sets.** v is minimal in v + 2Z^4 iff v/2 lies in the Voronoi cell
   of 0 (Voronoi, Crelle 134 (1908), §48; Conway and Sloane, *Low-dimensional
   lattices VI*, Theorem 2), iff 0 and v share a Delaunay cell. So the
   complete minimal sets, with multiplicities, are constant on each face's
   relative interior.
4. **Markings.** All 210 markings τ are treated (35 choices of
   W = ker N = im N, times 6 isomorphisms V/W → W). The statement is
   GL4(Z)-equivariant, with τ conjugated.

**Disposition of every case.** For a face with ray sum H(μ) = Σ μ_j R_j
(μ_j > 0), the six minimum equalities are linear: E μ = 0.

| Disposition (2,394 faces × 210 markings = 502,740) | Count | Certificate |
|---|---:|---|
| multiplicity mismatch | 477,492 | the two complete minimizing sets |
| lower face (no strictly positive solution) | 19,992 | an exact y with E^t y >= 0 and nonzero: a Stiemke certificate |
| certified | 5,256 | a rational μ > 0 with E μ = 0, an integral unimodular T ≡ τ, and a rational Z with **D = Z E**, where column j of D is T^t R_j T - R_j |
| uncertified | **0** | |

D = Z E shows that T preserves H(μ) for every real solution of E μ = 0,
irrational weights included. So the certificates cover all positive definite
real forms.

## 4. Verification in this release

The checks below use the Python standard library only, as everywhere in this
repository.

```sh
python code/rank4_theorem_s_check_published.py --output scratch/theorem-s/published.json
python code/rank4_theorem_s_check_bundle.py --workers 2 --output scratch/theorem-s/bundle-audit.json
```

**The published-cone check** verifies the explicit maps from both saved fans
onto the published Δ, Δ′, Δ″ (all ten rays each). It also checks that the
package sources and the exact fan have the hashes recorded with the audited
classification.

**The bundle check** recomputes the complete minimal sets of all 2,394 faces
with Grok Bot's exact enumerator. For all 502,740 cases it checks the stored
disposition:
- the minimizing sets;
- the Stiemke vectors;
- the positive solutions;
- T and D = Z E;
- that the certified set equals the original classification.

It does not use the classifier that found the certificates, cdd, or the
simplex solver.

Both run in CI. Locally: under a second, and about 36 s on two workers.

The package also contains:
- the generating code: `sources/rank4_six_cells.py` and
  `sources/rank4_secondary_fan.py`. These need numpy, scipy and pycddlib,
  unlike the rest of the repository. They regenerate the classification in
  about 10 s on 12 workers;
- Grok Bot's independent standard-library replay
  (`independent_replay_grokbot/replay.py`), which regenerates every
  disposition and certificate in about 7 minutes on one thread.

## 5. Provenance

| Item | Where |
|---|---|
| final audit (accepted) | [reviews/rank-four-theorem-s-final-audit.md](../reviews/rank-four-theorem-s-final-audit.md); the original in the package, `final_audit_astra/` |
| quick review, its resolution, the proof-status note | package `documents/` |
| independent replay | package `independent_replay_grokbot/` |
| literature and proof audit | package `literature_audit_grok/` |

**The fan files.** `data/fan.json` in the package is the first fan run. It
was **overwritten by mistake** and **regenerated** by the original
deterministic algorithm, so it is not the untouched original file. Its
consistency is shown in four ways:
- its central forms equal the recorded ones;
- every certified cell of `data/cells.json` lies on its stated face;
- the independent replay, made on a copy taken before the overwrite, gives
  identical results on it;
- its rays map onto the published cones.

`data/fan-exact.json` is the rerun with exact bounds and per-facet coverage
records.

**Errata recorded in the audit** (none affects the conclusion):
- the order-8 argument needs Φ8 = x^4 + 1;
- the plane-invariance multiplier is argued per characteristic (as above);
- Theorem 2.3 is in arXiv:math/0601084, while the explicit cones are in
  arXiv:2407.03513 section 2.1;
- the packaged classifier's `coverage_ok` counts facet records only. The
  proof rests on the audited fan hashes and the published-cone
  identification, both checked here.

**A related computed result**, marked vonorm rigidity in dimension four (of
which Theorem S is one case), is not part of this release. It awaits an
independent replay.
