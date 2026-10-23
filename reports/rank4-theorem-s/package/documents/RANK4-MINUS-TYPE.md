# Minus-type parity conjecture: stage 1 closed, stage 2 open

Item 7 of `A_what_is_next.txt`, 2026-10-21. Started by a Codex session
("Astra"), which wrote the scripts and produced the results below but froze
before writing this note; the note, the proof write-up, the replay of the
det 202-501 census, the symmetry audit of that census and an independent
check of the F2 lemma were completed by a Claude Opus 5.5 session the same
day. Python 3.11 standard library, exact arithmetic throughout. Notation as in
`RANK4-PARITY-CONJECTURE.md` (V = L/2L, Q, B, c = G^{-1} b mod 2, epsilon(g),
levels 1-4).

Conjecture (item 7). Every even, orthogonally indecomposable, positive
definite rank-four lattice with det = 5 mod 8 (minus type) has the parity
property P.

## 1. Verdict

- **Stage 1 is done (proved).** Theorem 5 below classifies symmetry zeros in
  every even, indecomposable, rank-four lattice of odd determinant: (a, b) is
  a symmetry zero iff some h in Aut(L) with h^2 = -I has a and c as distinct
  nonzero classes of its fixed plane (h+I)L/2L. Minus type has no such h, so
  **minus-type lattices have no symmetry zeros, at any determinant**. This
  closes the gap left open in section 5.4 of the parity report (witnesses of
  order 8, and order-4 witnesses acting with order 4 on L/2L).
- **Stage 2 is open.** Zeros without a symmetry witness (levels 2-4) are not
  excluded. Two natural constructions are ruled out (section 3); no proof
  route is established.
- **Census evidence extended to det 501.** All 703 indecomposable minus-type
  classes with 202 <= det <= 501 have P (section 4), independently replayed.
  With the 131 classes through det 201: all 834 indecomposable minus-type
  classes with det <= 501 have P.
- **Side result (not about minus type).** Halving the even det-4480 level-4
  example gives a certified **level-4 zero of determinant 280** (odd lattice),
  smaller than the det-1729 example the October reports cite (section 5).

## 2. Stage 1: classification of symmetry witnesses

Setting: L even, rank four, odd det, indecomposable; (a, b) an even
characteristic with a, b nonzero mod 2. Because det is odd, B is a
nondegenerate alternating form on V (B(x, x) = x.x = 0 mod 2), so every
automorphism g induces gbar in Sp(V) = Sp_4(F2), of order 720. A symmetry
witness is g in Aut(L) with ga = a mod 2L, g^{-t} b = b mod 2 (equivalently
gc = c mod 2L) and epsilon(g) = 1.

**Lemma A (order-4 symplectic maps fix a line).** Every element of order 4 in
Sp_4(F2) fixes exactly one nonzero vector.
Proof. Put N = gbar - I. As gbar has 2-power order in characteristic 2, N is
nilpotent; gbar^2 = I + N^2 != I, so N^2 != 0 and gbar has a Jordan block of
size at least 3. Size 4 gives a one-dimensional fixed space. Suppose the
blocks are 3 + 1: basis u, Nu, N^2 u, z with Nz = 0, so ker N = span(N^2 u, z),
im N = span(Nu, N^2 u). For x in ker N, B(x, gy) = B(gx, gy) = B(x, y), so
ker N is orthogonal to im N; hence N^2 u is orthogonal to N^2 u, Nu and z, and
nondegeneracy forces B(N^2 u, u) = 1. Invariance on (u, Nu):
B(u + Nu, Nu + N^2 u) = B(u, Nu) + B(u, N^2 u) + B(Nu, Nu) + B(Nu, N^2 u)
= B(u, Nu) + 1, which must equal B(u, Nu). Contradiction.
Checks: all 15^4 column choices, 720 symplectic maps, 180 of order 4, each
fixing one nonzero vector (`rank4_minus_type.py finite`); repeated
independently over all 2^16 matrices against the standard form J (same
table: 180 elements of order 4, all with one fixed nonzero vector).

**Lemma B (order 8 forces a square root of -I).** An element g of order 8 in
GL_4(Z) has characteristic polynomial x^4 + 1 (Phi_8 has degree 4 and must
divide it), so g^4 = -I and h = g^2 satisfies h^2 = -I.

**Theorem 5 (symmetry zeros for odd determinant).** Let L be even, rank four,
indecomposable, with odd determinant. An even characteristic (a, b) has a
symmetry witness iff there is h in Aut(L) with h^2 = -I such that a and c are
distinct nonzero classes of the plane (h+I)L/2L. In particular L has no
symmetry zero if it has minus type.
Proof. "If" is Theorem 4 of the parity report. Conversely let g be a witness.
By the corollary to Theorem 3 its 2-part g' = g^m (m odd) is a witness of
order 4 or 8, fixing a and c mod 2L.
(i) If g'bar has order 4, Lemma A gives a single fixed nonzero class, so
a = c, i.e. b = Ga mod 2, and Theorem 1 says Theta[a,b] != 0, contradicting
epsilon(g') = 1.
(ii) If g' has order 8, Lemma B gives h = g'^2 with h^2 = -I. If g'bar had
order <= 2, h would fix all of V, but by Theorem 4 its fixed space
(h+I)L/2L has index 4. So g'bar has order 4 and (i) applies.
(iii) If g' has order 4 and g'bar has order <= 2, then g'^2 = I mod 2L. If
g'^2 = -I we are in the stated case: a, c lie in the fixed plane of h = g'
(Theorem 4 identifies it with (h+I)L/2L) and are distinct by Theorem 1. If
not, g'^2 is an involution other than +-I congruent to I mod 2L, which splits
L orthogonally into its eigenlattices (x = (x + g'^2 x)/2 + (x - g'^2 x)/2
with both halves in L), contradicting indecomposability.
Minus type: by Theorem 4 any h with h^2 = -I forces plus type.

**Finite checks.** `rank4_minus_type.py census` re-derives, for every
automorphism of every stored class, the order, the mod-2 order and the fixed
classes. It computes every stabilized pair's sign (a primal and a dual
computation, asserted equal), and asserts the implications used above. Run on:

| Source (universe) | Classes | Indecomposable | Minus-type indec. | Symmetry zeros in minus type | All assertions |
|---|---:|---:|---:|---:|---|
| `results/rank4-parity/task-b-det201.json` (all even classes, odd det <= 201) | 300 | 255 | 131 | 0 | pass |
| `results/rank4-minus-type/minus-det202-501.json.gz` (section 4) | 787 | 703 | 703 | 0 | pass |

In the second run the 703 indecomposable groups contain 592 elements of
order 4, every one of mod-2 order 4 with one fixed nonzero class, and no
element of order 8 (as Lemma B and Theorem 4 predict for minus type).
Outputs: `odd-det201-symmetry-audit.json`, `minus-det202-501-symmetry-audit.json`,
`finite-groups.json` in `results/rank4-minus-type/`.

## 3. Stage 2: constructions tested (none can refute the conjecture)

Stage 2 asks about zeros with no witness: levels 2-4. Two constructions were
examined; neither can produce a minus-type counterexample. This is evidence
about those constructions only, not about stage 2.

**Conway-Sloane template (proved).** The common-kernel construction behind
the level-4 examples (`results/rank4-level4/task-a.json`, first pair) gives a
third lattice whose Gram, times 3, is linear in the four parameters:

    3M = [[9a, 3a, 3a, 0], [3a, a+b+d, a-b, d-b], [3a, a-b, a+b+c, b-c], [0, d-b, b-c, b+c+d]].

This identity, and the common-kernel isometry, were checked for all four
parameter unit vectors, which suffices by linearity
(`rank4_minus_type_constructions.py template`). Claim: whenever M is integral
and even, for any rational parameters (so any rational rescaling), every
entry of M is even, so 16 | det M and det M is never odd. Proof:
M12 = a is an integer, and M11 = 3a even makes a even. Integrality of
M23, M24, M34 makes b, c, d integers. Then M22 = a - 2 M23 + M24 even gives
M24 even, M33 = a - 2 M23 - M34 even gives M34 even, M44 = b + M24 - M34 even
gives b even, and M23 = (a - b)/3 is even because a - b is even and divisible
by 3.

**Quarter-turn probe (bounded, no hit).** Kernel M = diag(2A, 2A, 2B, 2C)
extended by (1,1,1,1)/2, with a fixed order-4 isometry T (T^2 != +-I). Each
coset pair (a, Ta) whose third overlattice is even was kept.
Universe: odd A, B, C in [1, 15] with 2A + B + C = 0 mod 4 and
A^2 BC = 5 mod 8: 128 parameter triples, 512 configurations with an even
third lattice, 59 distinct LLL-reduced Gram matrices (not deduplicated by
isometry). **0 indecomposable**: every hit is decomposable, so the probe says
nothing about the conjecture beyond this family
(`quarter-turn-probe15.json`).

**Mixed order-four kernel observation (finite check, unproved).** On every
kernel class with det = 4 mod 8 in `results/rank4-parity/task-c-det64.json`,
the check takes each order-4 automorphism T with T^2 != -I, and each class v
on which T acts as a 4-cycle (Tv != v, T^2 v != v). In every case some
involution of the same group also maps v to Tv: 32 kernel classes,
540 such T, 3,936 steps, no exception (`mixed-kernel-lemma-check.json`).
The universe matches the kernels of odd-det lattices (v_2(det L0) = 2).
Astra recorded no proof and no stated use. Presumably the aim was to reduce
level-2 configurations built from such T to involutive ones. Treat it as an
untested lead.

## 4. Census extension, 202 <= det <= 501

Universe: all isometry classes of even positive definite rank-four lattices
with 202 <= det <= 501 and det = 5 mod 8. The generator is the complete
reduced-basis builder of the parity report, with its determinant filter
narrowed (`rank4_minus_type.py extend`). It produced 23,612 candidates. Dedup
was by exact isometry (complete basis-image search between candidates with
equal invariants, which are det and vector counts of all norms <= 12),
giving 787 classes in 38 determinants. 84 classes are decomposable and 703
indecomposable, all of minus type.

Result: **703 of 703 indecomposable classes have P**. That is 73,815 = 703 x
105 characteristics, each proved nonzero by a complete shell coefficient
(first nonzero at norm4 <= 324), with 0 failures and 0 unresolved.
Decomposable classes were not decided. Generation and decision took 300 s.
Output: `minus-det202-501.json.gz`.

Independent replay (`rank4_minus_type.py replay-extension`, through
`verify_rank4_parity.py`): a loose-box candidate generator, box isometry
search and box coefficient enumeration, importing nothing from the builder.
Verdict `verified`: 23,612 candidates, 787 classes, 703 indecomposable,
703 hold, 0 fail, 0 unresolved, 101 s (`replay-minus-det202-501.json`,
log alongside). The first replay was launched by Astra and orphaned when that
session froze. It was stopped and rerun to completion from the same
unchanged source (SHA-256 `496c5c8c...1ec097d`).

Combined with Task B of the parity report: **every even indecomposable
minus-type rank-four lattice with det <= 501 has P (834 classes)**. Note that
item 7 advised extending the census only to test a stage-1 obstruction. This
extension is cheap and doubles as a test of Theorem 5's assertions, but a
larger table without failures adds little beyond that.

## 5. Side result: a level-4 zero of determinant 280

The even det-4480 level-4 example (`examples/rank4-level4/other-level4-certificates.json`,
row A4 CS (2,8,14,20) of `RANK4-LEVEL4-SEARCH.md` section 4) has every Gram
entry even. Dividing by 2 gives the primitive odd lattice

    G = [[3,1,1,0],[1,5,-1,2],[1,-1,4,-1],[0,2,-1,7]],  det 280,
    a = (1,0,0,1), b = (1,1,1,1),

whose characteristic theta vanishes (scaling the Gram only rescales tau).
It was certified directly rather than by appeal to scaling
(`rank4_minus_type_constructions.py det280`, `det280-rescaling.json`):
- Aut(L) of order 2 by two independent isometry searches, both elements
  epsilon 0;
- the zero proved by the modular cutoff, with the shells reproduced by an
  independent box enumeration;
- level 4, because the two coset lattices are not isometric (independent
  box search).

So the smallest known level-4 determinant is 280, not 1729. Because
Aut(L) = {I,-I} with epsilon 0 on both elements, it is also the smallest known
counterexample to the rank-four symmetry converse: the bracket drops from
24 < det <= 1156 to 24 < det <= 280. No minimality is claimed. This does not
touch the minus-type question (odd lattice, 4 | det). Dated correction notes
were added on 2026-10-21 to the level-4, Eisenstein and literature reports,
their reviews and the assessment, here and in `Github/`. The release carries
the certificate `data/rank4-level4/det280-level4.json` and the replay
`code/verify_rank4_det280.py`. A separate brute-force check agrees: all 383
signed shells through norm4 3072 cancel, over 173,956 vectors.

Addendum, 2026-10-21 (later the same day): minimality is now proved.
`RANK4-GAP-CENSUS.md` lists every isometry class with 25 <= det <= 279
(14,987 new classes in four regions, independently replayed) and finds every
zero at level 1, so 280 is the minimum determinant of a rank-four
counterexample to the symmetry converse.

## 6. Files

New scripts: `rank4_minus_type.py` (tasks finite, census, mixed-kernels,
quarter-turn, extend, replay-extension; `census` now also reads .gz) and
`rank4_minus_type_constructions.py` (template, det280). Outputs in
`results/rank4-minus-type/`. Every output file refuses to overwrite. No
earlier file was modified.

Regenerate (from the workspace root):

    python -B rank4_minus_type.py finite --output <out>.json
    python -B rank4_minus_type.py census --output <out>.json
    python -B rank4_minus_type.py census --source results/rank4-minus-type/minus-det202-501.json.gz --output <out>.json
    python -B rank4_minus_type.py mixed-kernels --output <out>.json
    python -B rank4_minus_type.py quarter-turn --output <out>.json
    python -B rank4_minus_type.py extend --dmin 202 --dmax 501 --output <out>.json.gz
    python -B rank4_minus_type.py replay-extension --source results/rank4-minus-type/minus-det202-501.json.gz --output <out>.json
    python -B rank4_minus_type_constructions.py template --output <out>.json
    python -B rank4_minus_type_constructions.py det280 --output <out>.json

## 7. What next

Stage 2 is the whole remaining problem. The route proposed in item 7 is to
show that every zero in this setting forces span(a, c) to be a totally
singular plane. It is still unproved, and stage 1 does not supply it: Theorem 5 proves it only for
symmetry zeros. The det-1156 and det-1729/280 examples show that zeros
without witnesses exist outside minus type. So a proof must use the
quadratic form on V, not automorphisms. Independent human review of
Lemma A and Theorem 5 is needed before either is cited (item 11).

Addendum, 2026-10-21 (night): `RANK4-STAGE2-STRATEGY.md` (a Codex session)
reduces stage 2 further, and `RANK4-STAGE2-REVIEW.md` reviews it and checks
Lemma A and Theorem 5 line by line (its section 7). If the reductions hold,
every minus-type zero would be level 4: two non-isometric isospectral coset
lattices. The census of even lattices of odd determinant now runs to
det 501 in both types (1,890 classes, 1,696 indecomposable). In that range
every zero of an indecomposable class lies on the plane of a square root
of -I, as Theorem 5 predicts.
