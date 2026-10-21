# Review of the stage-two reductions

2026-10-21, Claude Opus 5.5 (the session that proposed the combined
conjecture). AI-assisted local review with independent computation. This is
not peer review; independent human review is still outstanding. The reviewed
note, [rank-four-stage-two-reductions.md](../docs/rank-four-stage-two-reductions.md), is left
unchanged apart from a pointer to this review.

## 1. Verdict

| Section | Claim | Verdict |
|---|---|---|
| 1 | plus-type test range is det = 1 mod 8, 281..497; compare all zeros with the union of all h-planes | **accepted** (gap census table confirms 202..279 is covered: 269 even odd-det classes, 238 indecomposable, 0 zeros) |
| 2 | plane invariance via Gamma = [[2I, 3G], [H, delta I]] | **correct**; it also fixes an error in my proposal (ordinary S passes to G^{-1}, not to the same lattice) |
| 3 | det-21 control (decomposable, one singular vector, zero) | **correct**; it refutes the Gauss-sum shortcut as I stated it |
| 4 | mixed order-four lemma | **correct** (proof checked line by line, plus independent computation) |
| 5 | odd nonsquare det: mapping implies swapping, so no strict level 2 | **correct** |
| 6 | minus type: O1, O2 never isometric, so every zero is level 4 | **correct**; it uses Theorem 5 of the [minus-type report](../docs/rank-four-minus-type.md), checked in section 7 below |

No errors found. The reductions are sound as written. Small presentation
points are in section 3.

## 2. What was checked, by hand

**Section 2.** A^tC = 2H and B^tD = 3 delta G are symmetric, and
A^tD - C^tB = (2 delta - 3D) I = I, so Gamma is symmetric. AOmega + B = (2tau+3)G
and COmega + D = (D tau + delta) I, so Omega = tau G maps to a scalar multiple of G.
D = 1 mod 4 makes delta even. In the standard characteristic action
(m' -> Dm' - Cm'' + diag(CD^t), m'' -> -Bm' + Am'' + diag(AB^t)) the new pair
is (delta a - Hb, -3Ga + 2b) = (c, Ga) mod 2. The diagonal terms diag(delta H)
and diag(6G) are even. The Moebius map has determinant 1, so the theta
function vanishes identically before the transformation iff it does after.

**Section 3.** Recomputed: B(a,c) = 0; Q(a), Q(c), Q(a+c) = 0, 1, 1; each A2
and [[2,1],[1,4]] factor carries an odd characteristic. My proposed mod-8
split is a valid identity: with y = a + 2x, y^2 = a^2 + 4(a,x) mod 8 because
L is even. But it uses only the plane, so it cannot separate this control
from an indecomposable lattice. Accepted.

**Section 4.** Checked in order:
- The eigenspace dimensions are (2, 2): the +-i eigenvalues pair up, and
  g^2 != -I.
- Gram(N) = cI: e.f = ge.gf = -f.e. And c is even because M is even.
- 2M lies in N perp P.
- The glue maps injectively into each discriminant group (N and P are
  primitive) and is an anti-isometry of the discriminant bilinear forms
  (because M is integral).
- 2s + t - 2k = 2.
- k = 2: (s,t) = (1,4) needs three generators in a binary discriminant
  group; (s,t) = (2,2) pairs a zero form against a nondegenerate one. Both
  cases are excluded, and s >= 3 is impossible.
- k = 1 forces (1,2). The glue's N component (e+f)/2 is the only nonzero
  2-torsion class fixed by g.
- R_+- preserve M.
- The difference formulas (alpha+beta)(e-f) and (beta-alpha)(e+f) hold, and
  the two coefficients have opposite parity.
The lemma uses only M even and v_2(det M) = 2.

**Section 5.**
- With v_2(det M) = 2, M has at most two even index-two overlattices
  (cyclic case: one subgroup of order 2; Klein case: three isotropic
  elements would force the bilinear form to vanish).
- The orbit v1 -> v1+w -> v1+w+u -> v1+u is right.
- By CRT there is j = 1 mod 4 with g^j of 2-power order, so g^j v1 = v2.
- 2-power orders in GL_4(Z) are at most 8.
- A Z[i]-structure gives det M = (det A)^2 for an integral Hermitian A (the
  realification of a Hermitian matrix has determinant |det A|^2). This is
  impossible for det M = 4d with d not a square.
- A swap fixes w = v1 + v2, so it lies in Aut(L) and fixes a mod 2L, which
  is the level-1 definition of the [level-four report](../docs/rank-four-level-four-zero.md) section 1.

**Section 6.**
- Q(a) = 1 makes every vector of Oi minus M have half-integral norm, so M
  is the integral-norm part of each Oi and every isometry O1 -> O2 lies in
  Aut(M).
- The dual Klein configuration (S^* inside L^*, O1^*, O2^* inside M^*) has
  sign vector a and coset vector b. So the dual characteristic is (b, a) and
  its c' is Ga mod 2.
- H = dG^{-1} is even, because a 3x3 even symmetric matrix is an odd-size
  alternating matrix mod 2. Also det H = d^3 = 5 mod 8, and Q_H(b) = Q_G(c).
- The dualization is needed. Plane invariance only shows that (c, a) is also
  a zero; it says nothing about the level of (a, c).

## 3. Presentation points (no change in substance)

1. Theorem 5 is the only step that uses indecomposability, and it is used
   only to exclude level 1. Theorem 5 reduces a witness to an h with
   h^2 = -I, and such an h forces det L to be a square (same Hermitian
   argument). Without indecomposability, level-1 zeros do occur at nonsquare
   determinant (the 1,440 decomposable cases in section 4). This partly
   answers the strategy note's question of where indecomposability enters:
   for levels 1-3 it enters through Theorem 5, and a level-4 argument must
   use it again. The chain therefore rests on Lemma A and Theorem 5, and
   both need human review before citation.
2. Section 6's argument also applies to plus type at nonsquare determinant
   (det = 1 mod 8, not a square), on planes that are **not totally singular**:
   there too, every zero is level 4. In plus type at nonsquare determinant,
   any zero is therefore either level 4, or level 3 or 4 on a totally
   singular plane. Section 7 of the strategy note is correct that plus type
   needs more, but this much carries over.
3. The phrase "every relevant characteristic" in section 6 should say "every
   even characteristic with a != c". Any O1 that is isometric to O2 is
   isospectral with it, so the statement is equivalent to "no zero of level
   <= 3".

## 4. Independent computation

Script `code/rank4_stage2_review_check.py` (own LLL, Fincke-Pohst and
basis-image automorphism search, importing nothing from the project; SHA-256
`0f8c8760...c68fb5ac0`). Output: `reports/rank4-stage2/review-independent-check.json`.

Universe: the 300 Gram matrices of `reports/rank4-parity/task-b-det201.json`
(SHA-256 `c462ed47...af93a904`). These are all isometry classes of even
positive-definite rank-four lattices with odd det <= 201, deduplicated by
exact isometry in that census. Decomposable classes are included; only the
Gram matrices are read. For each L, each of the 15 sign kernels M, and each
of the 7 pairs {v, v+w} (31,500 configurations), Aut(M) was recomputed from
scratch. Every recomputed |Aut(L)| spot-checked against the census matched.

| | nonsquare det (288 lattices) | square det (12 lattices) |
|---|---:|---:|
| configurations | 30,240 | 1,260 |
| some g in Aut(M) maps O1 onto O2 | 1,440 | 186 |
| ... and some g swaps them | 1,440 | 186 |
| mapping without swap (section 5 claim) | **0** | 0 |
| elements of Aut(M) with g^2 = -I | **0** | 66 |
| mixed order-four elements | 1,824 | 440 |
| four-cycle classes, involution R with Rv = gv found (section 4) | 14,592 / 14,592 | 3,520 / 3,520 |

The test is not vacuous. All 1,440 nonsquare mapping cases lie on
decomposable lattices (census flag `decomposable`). Of the 186 square cases,
180 are decomposable and 6 are the det-121 Gaussian zeros. The section 4 count covers kernels
of every det = 4d, d odd <= 201, well beyond the 32 kernels through det 64
of the original observation. The square-determinant column lies outside the
proof's hypothesis; there, too, no mapping-without-swap case occurs.

## 5. Consequence

If Lemma A, Theorem 5 and sections 4-6 also survive human review, the
minus-type conjecture is equivalent to:

> no even, indecomposable lattice L with det L = 5 mod 8 has a characteristic
> whose two coset lattices O1, O2 (Klein configuration around the sign kernel
> M) are isospectral.

Any such pair would be non-isometric, i.e. a level-4 pair of Kneser
2-neighbours. With Q(a) = 1 (allowed by plane invariance), M is the
integral-norm part of both. Equivalently, sqrt(2) O1 and sqrt(2) O2 are
odd integral lattices of determinant 16d with the same even sublattice
sqrt(2) M. That is the precise residual problem.

## 6. Plus-type test (strategy note section 1), executed

`code/rank4_plus_type_census.py` (this session) reuses `rank4_gap_census.build_chunk`
and `replay_chunk` unchanged, with a new region and output folder
`reports/rank4-plus-type/`.

Universe of the new run: every isometry class of even positive definite
rank-four lattices with det = 1 mod 8, 281 <= det <= 497 (28 determinants,
one chunk each). Complete reduced-basis generator, 24,580 candidates, dedup
by exact isometry within each determinant: 672 classes, 49 decomposable,
623 indecomposable, all 105 even characteristics of each decided. 622 have P.
1 fails: det 361 = 19^2, G = [[4,-2,0,-1],[-2,6,-1,0],[0,-1,4,2],[-1,0,2,6]],
6 zeros, all level 1. No unresolved characteristic. Build 225 CPU s
(21 s wall, 16 workers). Independent replay (`verify_rank4_parity.replay_decided`
per chunk) verified all 28 chunks (16 s wall).
Files: `even-plus-281-497/det*.json.gz`, `-summary.json`, `-replay.json`, `run.log`, `replay.log`.

Zero-set comparison (`planes` task, `reports/rank4-plus-type/planes.json`):
for each indecomposable class, the set of zero characteristics (a, c),
c = G^{-1} b mod 2, was compared with the union over **all** h in Aut(L)
with h^2 = -I of the six ordered pairs of distinct nonzero classes of
(h+I)L/2L. This is the comparison the strategy note asked for.

| Source (even lattices, odd det) | Classes | Indecomposable | With h^2 = -I | Zeros | Predicted | Unexplained zeros | Predicted but nonzero |
|---|---:|---:|---:|---:|---:|---:|---:|
| Task B, det <= 201 | 300 | 255 | 1 (det 121) | 6 | 6 | 0 | 0 |
| gap census, 202..279 (both types) | 269 | 238 | 0 | 0 | 0 | 0 | 0 |
| minus extension, det = 5 mod 8, 202..501 | 787 | 703 | 0 | 0 | 0 | 0 | 0 |
| this run, det = 1 mod 8, 281..497 | 672 | 623 | 1 (det 361) | 6 | 6 | 0 | 0 |

The gap-census and minus-extension rows overlap on minus type 202..279
(123 indecomposable classes). Without double counting, the union is
**every even positive definite rank-four class with odd det <= 501:
1,890 classes, 1,696 indecomposable.** In that universe the zero set of
every indecomposable class equals the h-plane prediction exactly. So the
combined conjecture holds through det 501. In particular no plus-type zero
outside the square-root-of-minus-I mechanism exists there. The only
indecomposable zeros are the 12 at det 121 and det 361, both square
determinants, as section 5 of the strategy note requires. Per the strategy
note, the bound is not raised further; the remaining work is the level-4
problem of section 5 above.

## 7. Check of Lemma A and Theorem 5 ([minus-type report](../docs/rank-four-minus-type.md) section 2)

The reductions of sections 5-6 use Theorem 5 to exclude level 1. This
section checks its proof line by line. It was not checked separately before.

**Lemma A.**
- N = gbar - I is nilpotent, and N^2 != 0 because gbar^2 = I + N^2.
- For x in ker N, B(x, Ny) = B(gx, gy) - B(x, y) = 0. So ker N is
  orthogonal to im N.
- In the 3 + 1 Jordan case, N^2 u lies in both, and is orthogonal to N^2 u,
  Nu and z. Nondegeneracy then gives B(N^2 u, u) = 1.
- Invariance on the pair (u, Nu) gives B(u, Nu) = B(u, Nu) + 1, a
  contradiction.
- The exhaustive check (180 elements of order 4 in the 720-element group,
  each fixing exactly one nonzero vector) is consistent: Sp_4(F2) = S_6, and
  its 180 elements of order 4 are the 90 4-cycles and the 90 products of a
  4-cycle with a disjoint transposition.

**Lemma B.** Phi_8 has degree 4, so an element of order 8 in GL_4(Z) has
characteristic polynomial x^4 + 1 and g^4 = -I.

**Theorem 5.**
- *Reduction.* The corollary to Theorem 3 reduces a witness to its 2-part g'
  of order 4 or 8. That corollary is reviewed in [rank-four-parity-review.md](rank-four-parity-review.md).
- *Case (i).* Lemma A gives a = c, and Theorem 1 then makes the
  characteristic nonzero.
- *Case (ii).* For h = g'^2 with h^2 = -I: (h+I)^2 = 2h, so mod 2 the image
  of h+I lies inside its kernel. Theorem 4 gives an image of dimension 2, so
  the fixed space ker(hbar + I) is 2-dimensional and hbar != I. So g'bar
  cannot have order <= 2. Its order is 4, since S_6 has no element of order
  8, and case (i) applies.
- *Case (iii).* If g'^2 = -I, the fixed space of g'bar is ker(g'bar + I),
  which equals the plane (h+I)L/2L by the dimension count just given. If
  instead g'^2 = s is an involution other than +-I with s = I mod 2L, then
  (x - sx)/2 and (x + sx)/2 both lie in L. So L is the orthogonal sum of the
  two eigenlattices, contradicting indecomposability.
- *Minus type.* Theorem 4 makes any h with h^2 = -I plus type.

Verdict: Lemma A, Lemma B and Theorem 5 are correct as written. The audits
in section 2 of the minus-type report support the proof (834 indecomposable
minus-type classes, every automorphism, every stabilized pair). With this
check, the chain Lemma A -> Theorem 5 -> sections 4-6 has been reviewed end
to end. It is still an AI-assisted local review, not peer review.
