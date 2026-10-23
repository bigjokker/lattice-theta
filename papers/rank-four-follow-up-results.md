# Rank-four follow-up results: non-isometric coset lattices, the minimum counterexample determinant, and the parity property

Releases 1.1.0 to 1.4.0, 2026-10-21 to 2026-10-23. This note states the results added after
the [counterexample manuscript](rank-four-counterexample.md), with pointers to
the full reports, the review notes, the exact certificates and the
reproduction commands. The computations were specified, run and reviewed by
Claude Fable 5.1 sessions and then assessed by a separate Claude session
([reviews](../reviews/README.md)); the determinant census of section 4 was
run after the assessment and is verified by independent replay rather than
by a separate review note. The odd-determinant results of section 6 were
obtained by a Codex session and Claude Opus 5.5 sessions and are reviewed in
[rank-four-stage-two-review.md](../reviews/rank-four-stage-two-review.md).
Theorem S of section 7 (release 1.4.0) is a computer-assisted proof, audited
in [rank-four-theorem-s-final-audit.md](../reviews/rank-four-theorem-s-final-audit.md).
These are AI-assisted local reviews with
independent computational replication, not peer review. No publication,
priority or novelty claim is made. Every count below names its enumeration
universe; "class" always means an isometry class, with deduplication by a
complete basis-image isometry search.

## 1. Setting: an even zero is an isospectral pair

L is a positive definite integral lattice of rank four with Gram matrix G, and
(a, b) is an even characteristic: binary columns, both nonzero modulo two,
with a.b even. Let L0 = {x in L : x.b even}, t in L \ L0 and xi = a/2. The
sign kernel L0 has index two in L, and the two coset lattices

    L1 = L0 + Z xi,    L2 = L0 + Z (xi + t)

satisfy Theta_G[a,b] = i^(a.b) (theta_L1 - theta_L2) up to the common
theta_L0 ([foundations](../docs/foundations-and-exact-enumeration.md), section
8; [literature note](../docs/rank-four-isospectrality-literature.md), section
1). So the characteristic vanishes identically exactly when L1 and L2 are
isospectral. A symmetry witness (an automorphism of L with sign one on the
stabilizer) is an element of Aut(L0) that interchanges L1 and L2. The possible
explanations of a zero therefore form a hierarchy:

| Level | Explanation |
|---|---|
| 1 | an element of Aut(L0) interchanges L1 and L2: a symmetry witness |
| 2 | an element of Aut(L0) maps L1 onto L2 but no element interchanges them |
| 3 | an isometry L1 -> L2 exists, but none preserves L0 |
| 4 | L1 and L2 are isospectral and not isometric |

Ranks one to three are always level 1 (the symmetry converse, proved in the
[low-rank manuscript](low-rank-and-root-lattice-results.md)). The
determinant-23,716 counterexample of the manuscript is level 2. The results
below exhibit level 4 in rank four, determine the minimum determinant of a
counterexample, give a smaller level-2 example, settle two parity
conjectures negatively, and reduce the surviving minus-type parity
conjecture to level-four zeros.

## 2. Level four: a zero whose coset lattices are not isometric

**Theorem A.** Let G = [[3,1,1,0],[1,9,-2,4],[1,-2,7,-2],[0,4,-2,13]],
a = (1,0,0,1), b = (1,1,1,1). Then L is odd, primitive, positive definite of
determinant 1729 = 7 * 13 * 19, the characteristic is even, Theta_G[a,b]
vanishes identically, and the coset lattices L1, L2 are isometric to the
Conway-Sloane lattices L^-(1,7,13,19) and L^+(1,7,13,19), Schiemann's
isospectral pair, which are not isometric. Aut(L) = Aut(L0) = {I, -I}.

The vanishing is proved twice at every level: by the Conway-Sloane sign-change
bijection, checked exactly on the nine relevant cosets, and by the project's
weight-two modular cutoff (B = 17920 in norm4 units) with complete signed
shells computed by three independent enumerators. Non-isometry of L1 and L2
is a complete basis-image search replayed by integer boxes, together with an
exact isometry invariant; it is also Cervino-Hein's theorem. The construction
is the observation that the Conway-Sloane pair, re-embedded, is a 2-neighbour
pair: the index-two sublattices of all-even coordinates in the two lattices
are isometric, and transporting one lattice along that isometry puts the pair
in a Klein configuration whose third intermediate lattice is L. The
isospectral pairs and their non-isometry are Conway-Sloane's and
Cervino-Hein's; the re-embedding and the characteristic zero are what this
project adds.

The same construction gives eight level-4 configurations: one for each
Theorem-1 pair tried, with parameters (1,7,13,19), (1,7,13,43), (1,7,25,31),
(1,7,13,31), (1,7,19,25), (2,8,14,20), (7,13,19,49), of determinants 1729,
3913, 5425, 2821, 3325, 4480, 84721, and the Earnest-Nipp pair, which is the
(7,13,19,49) configuration in another basis. All eight are certified
([data/rank4-level4](../data/rank4-level4/)). Universes searched for other
levels: the 450 index-two cells of the nine known pairs (Task A) and the 169
sign kernels with det(L0) <= 24 from the census (Task D: 17,745 pairs, each
either nonzero or level 1). Finite-domain statement: no level-2, level-3 or
level-4 configuration has det(L0) <= 24. No level-3 configuration is known.

**Theorem A'' (the smallest level-four zero).** Let
G = [[3,1,1,0],[1,5,-1,2],[1,-1,4,-1],[0,2,-1,7]], a = (1,0,0,1),
b = (1,1,1,1). Then L is odd, primitive, positive definite of determinant
280 = 2^3 * 5 * 7, Aut(L) = {I, -I} with epsilon 0 on both elements,
Theta_G[a,b] vanishes identically, and the coset lattices L1, L2 are not
isometric: a level-4 zero with no symmetry witness. By Theorem D no lattice of
smaller determinant has a zero of level above one, so this is the smallest
level-four zero and the smallest counterexample to the symmetry converse of
any level.

G is half the Gram matrix of the determinant-4480 configuration above
(parameters (2,8,14,20)), all of whose entries are even; scaling a Gram matrix
changes neither the vanishing nor the isometry type of the coset lattices.
The lattice is also certified directly: the modular cutoff (B = 3072, 383
complete signed shells, 173,956 vectors, every signed coefficient zero), the
full automorphism group and the non-isometry of the coset lattices by
complete searches ([det280-level4.json](../data/rank4-level4/det280-level4.json),
replay `code/verify_rank4_det280.py`, which first checks the determinant-4480
certificate with the unchanged level-four verifier). Its coset lattices are
not integral.

Precision. "The vanishing has no geometric explanation" means no single
Euclidean isometry between the coset lattices explains it; the sign-change
bijection is a piecewise norm-preserving bijection, so a geometric description
does exist.

**Corollary A' (a determinant-3844 pair).** The only agreeing index-two cell
of the Conway-Sloane remark-(x) pair (determinant 961) gives the even lattices
S1 = [[4,1,0,0],[1,8,0,0],[0,0,8,2],[0,0,2,16]] (decomposable) and
S2 = [[4,-2,1,1],[-2,8,-3,-1],[1,-3,10,-4],[1,-1,-4,18]] of determinant
3844 = 62^2. Both have level 62 and square determinant, so their theta series
are weight-two modular forms on Gamma_0(62) with trivial character; the index
is 96, the inclusive Sturm bound 16, and the coefficients 0 through 16 agree
(two independent enumerators), while no isometry exists (two independent
searches). So S1, S2 is an isospectral non-isometric pair at every level.
Whether it appears in the literature has not been checked.

Full report: [rank-four-level-four-zero.md](../docs/rank-four-level-four-zero.md);
review: [rank-four-level-four-review.md](../reviews/rank-four-level-four-review.md).

## 3. An even symmetry counterexample: determinant 1156

**Theorem B.** Let G = [[4,-2,1,-2],[-2,8,2,-2],[1,2,6,-4],[-2,-2,-4,12]],
a = (0,1,0,1), b = (1,0,1,0). Then L is even, primitive, positive definite of
determinant 1156 = 2^2 * 17^2, Aut(L) = {I, -I}, so no symmetry witness
exists and L is indecomposable, and Theta_G[a,b] vanishes identically. It is
a level-2 zero: the sign kernel L0 is the trace form of the binary Hermitian
lattice over the Eisenstein integers with parameters (p, q, delta) =
(8, 6, 4 + 10 omega); the order-three isometry of L0 given by multiplication by
omega carries L1 onto L2, and no automorphism of L0 interchanges L1 and L2.

The vanishing is proved three ways: by the order-three isometry (a complete
all-level proof), by the modular cutoff (B = 144; 16 complete signed shells,
192 vectors, computed by three enumerators, every signed coefficient zero), and
by the manuscript's own ambient-coset certificate schema, replayed by the
unchanged `code/rank4_counterexample.py`. The determinant-23,716 example of
the manuscript is the same construction with ambient determinant 77^2; this
one has ambient determinant 17^2.

Family statement. In the Eisenstein family F (line characteristics of
index-two overlattices of trace forms of binary Hermitian lattices over
Z[omega]), the complete enumeration for 5 <= D' <= 41 (498 reduced triples,
914 members) contains exactly 12 level-2 members, all of determinant 1156 and
all isometric to G; every other member is level 1. The extension to
D' <= 154 (9659 members, 133 at level 2) confirms 1156 as the family minimum
and recovers the determinant-23,716 example. Determinant 1156 is therefore the
minimum of the Eisenstein family and the smallest known even counterexample.
It is not the global minimum: by Theorem D that is 280, attained by the odd
lattice of Theorem A''.

Full report: [rank-four-eisenstein-family.md](../docs/rank-four-eisenstein-family.md);
review: [rank-four-eisenstein-review.md](../reviews/rank-four-eisenstein-review.md);
certificate: [det1156-level2.json](../data/rank4-eisenstein/det1156-level2.json).

## 4. The minimum counterexample determinant is 280

**Theorem D.** Every positive definite integral rank-four lattice with
det <= 279 satisfies the symmetry converse: each of its identically vanishing
even characteristics has a symmetry witness. Hence 280 is the minimum
determinant of a rank-four counterexample. At determinant 280 exactly one
isometry class fails the converse, the lattice of Theorem A'', and it fails
on exactly one characteristic.

Proof. A lattice is even or odd, and an even rank-four lattice has det = 0 or
1 mod 4. Decomposable lattices satisfy the converse at every determinant
([decomposition appendix](../docs/rank-four-decomposition-and-fibres.md),
from the rank one to three converse). The determinant-24 census settles
det <= 24 (169 classes, every zero with a checked witness). For
25 <= det <= 279 the indecomposable classes are covered by seven complete
class lists, three from release 1.1.0 (odd lattices with 4 !| det through
168; even lattices of odd determinant through 201; even lattices with 4 | det
through 64; 183 zeros, all witnessed) and four new ones:

| Universe | Candidates | Classes | Indecomposable | Even zeros |
|---|---:|---:|---:|---:|
| odd lattices, 4 !\| det, 169..279 | 106,794 | 6,658 | 2,900 | 157 |
| odd lattices, 4 \| det, 25..279 | 98,629 | 7,043 | 2,577 | 4,042 |
| even lattices, odd det, 202..279 | 8,972 | 269 | 238 | 0 |
| even lattices, 4 \| det, 65..279 | 22,080 | 1,017 | 524 | 1,131 |

Each new list is generated one determinant at a time by the complete
reduced-basis builder of the parity report and deduplicated by exact
isometry. All 655,095 even characteristics of the 6,239 new indecomposable
classes are decided: 5,330 are zeros, and every zero has a symmetry witness
in the full automorphism group; the other 649,765 are nonzero by a complete
shell coefficient (the first nonzero shell has norm4 at most 238, so the
modular cutoff was never needed); nothing is unresolved. The independent
verifier `code/verify_rank4_parity.py` (a separate candidate generator, box
isometry searches, box enumeration, about 810,000 pairwise non-isometry
searches) replays all 239 chunks as `verified`. The census of determinant 280
itself (243 odd and 27 even classes, 187 zeros) finds exactly one zero without
a witness, of level 4, in a class isometric to the lattice of Theorem A''.

The universe of the claim is all positive definite integral rank-four
lattices with det <= 279, with no primitivity or indecomposability
restriction. Cost: about 50 CPU-minutes to build and 60 to replay, about ten
minutes on 16 cores.

Full report: [rank-four-minimum-determinant.md](../docs/rank-four-minimum-determinant.md);
script `code/rank4_gap_census.py`; preserved chunks, summaries, replays and
closure: `reports/rank4-gap-census/`.

## 5. The parity property in rank four

L has the parity property P when every even characteristic (a, b nonzero
modulo two, a.b even) has a theta series that is not identically zero. A
decomposable lattice never has P ([foundations](../docs/foundations-and-exact-enumeration.md),
section 7). The determinant-24 census suggested two conjectures:

- (C0) an odd (not even) indecomposable rank-four lattice whose determinant is
  not divisible by four has P;
- (C1) an even indecomposable rank-four lattice of odd determinant has P.

**Theorem C (both conjectures are false, with smallest determinants).**

(i) Universe: isometry classes of odd positive definite rank-four lattices with
25 <= det <= 168 and 4 not dividing det (46,494 reduced-basis candidates,
3,487 classes, 2,343 decomposable, 1,144 indecomposable). Exactly 13
indecomposable classes fail P, at determinants 49, 70, 91, 91, 119, 126, 130,
130, 133, 133, 154, 154, 161. The smallest is

    G = [[3,-1,-1,0],[-1,3,0,-1],[-1,0,3,1],[0,-1,1,3]],   det 49,

with automorphism group of order 16 and seven zeros, all level 1 with
involution witnesses; all 55 zeros of the 13 classes have involution
witnesses. With the census (every odd-lattice failure with det <= 24 has
4 | det), 49 is the smallest determinant of an indecomposable odd lattice
without P and with 4 not dividing det.

(ii) Universe: isometry classes of even positive definite rank-four lattices
with odd det <= 201 (10,188 candidates, 300 classes, 45 decomposable, 255
indecomposable). Exactly one indecomposable class fails P:

    G = [[4,-2,-1,0],[-2,4,0,1],[-1,0,4,-2],[0,1,-2,4]],   det 121 = 11^2,

the trace form of the Hermitian lattice [[4, 2+i],[2-i, 4]] over the Gaussian
integers. Its automorphism group has order 24; it has exactly six zeros, the
ordered pairs of distinct nonzero classes of the totally singular plane
(1+i)L/2L, each witnessed by an element g with g^2 = -I (multiplication by i).
The other 254 indecomposable classes hold P with every one of their 105
characteristics proved nonzero by a complete shell.

**Structural theorems** (proved in the report, checked line by line in the
review). For L even:

1. If b = G a modulo 2, every shell of a + 2L has constant sign, so (a, b) is
   never a zero (any rank).
2. For odd determinant and a different from G^{-1} b mod 2, neither coset
   lattice is integral, so classification results for integral lattices do not
   apply directly to the pair L1, L2.
3. For rank four, odd determinant and L indecomposable, no involution of L is
   a symmetry witness: every symmetry zero needs an automorphism whose 2-part
   has order 4 or 8.
4. For rank four and odd determinant, an automorphism g with g^2 = -I has
   (g + I)L / 2L a totally singular plane, L is of plus type (det = 1 mod 8),
   and all six ordered pairs on that plane are zeros witnessed by g. This is
   the rank-four analogue of the E8 mechanism, checked exactly on E8.

Consequently a minus-type lattice (det = 5 mod 8) has neither an involution
witness nor a square root of -I. Data: all 131 distinct indecomposable
minus-type classes with det <= 201 hold P. Whether every even indecomposable
rank-four lattice with det = 5 mod 8 has P is open; it is the first item of
[deferred research](../docs/deferred-research.md). Theorem E below excludes
the remaining witnesses (order 4 acting faithfully on L/2L, order 8), so only
zeros without any witness remain.

Degenerate case: among the 69 classes of even lattices with 4 | det <= 64,
27 are indecomposable, 7 hold P and 20 fail; their 122 zeros are all level 1
and all have involution witnesses. Of the 18 members of the nine known
rank-four isospectral pairs, 17 hold P; the exception, the first lattice of
the Conway-Sloane remark-(x) pair, is decomposable.

**Corollary C' (the ordinary theta series does not determine P in rank
four).** The two lattices of the Conway-Sloane remark-(x) pair of determinant
961 have the same theta series at every level (their doubles are even of
level 124 and square determinant; the index of Gamma_0(124) is 192, the
inclusive Sturm bound 32, and the coefficients 0 through 32 agree by two
independent enumerators) and are not isometric (two independent searches).
The first member is decomposable and has 36 vanishing even characteristics;
the second is indecomposable and has none. The statement concerns rank-four
lattices in general; it does not settle the question restricted to pairs of
indecomposable lattices.

Full report: [rank-four-parity-property.md](../docs/rank-four-parity-property.md);
review: [rank-four-parity-review.md](../reviews/rank-four-parity-review.md);
certificates: [c1-counterexamples.json](../data/rank4-parity/c1-counterexamples.json),
[pair-members-parity.json](../data/rank4-parity/pair-members-parity.json.gz);
corollary script: [code/check_small_pairs.py](../code/check_small_pairs.py).

## 6. Even lattices of odd determinant: symmetry zeros and the level reduction

Setting: L even, rank four, odd determinant d (so d = 1 mod 4), V = L/2L with
the nondegenerate alternating form B(x, y) = x.y mod 2 and the quadratic form
Q(x) = x.x/2 mod 2; L has plus type for d = 1 mod 8 and minus type for
d = 5 mod 8. For an even characteristic (a, b) put c = G^{-1} b mod 2; then
B(a, c) = 0, and Theorem 1 makes (a, b) nonzero when a = c.

**Theorem E (symmetry zeros).** Let L be even, orthogonally indecomposable,
of rank four and odd determinant. An even characteristic (a, b) has a
symmetry witness if and only if some h in Aut(L) with h^2 = -I has a and c as
distinct nonzero classes of the plane (h+I)L/2L. In particular a minus-type
lattice has no symmetry zero, at any determinant.

Proof outline. Every element of order four in Sp_4(F2) fixes exactly one
nonzero vector (Lemma A: a Jordan-block argument for the symplectic form,
also checked on all 720 elements). An element of order eight in GL_4(Z) has
characteristic polynomial x^4 + 1 (Lemma B). By Theorem 3 a witness may be
replaced by its 2-part, of order 4 or 8. If that part acts with order four on
V it fixes only one nonzero class, forcing a = c. An order-eight witness
squares to a root of -I, which acts nontrivially on V by Theorem 4. A witness
of order four acting with order at most two on V either squares to -I (the
stated case) or squares to an involution congruent to I mod 2L, which splits
L. The "if" direction is Theorem 4. Exhaustive audits of every automorphism
of every class in the even odd-determinant lists below agree.

**Theorem F (level reduction).** Let L be even, of rank four and odd
determinant d, and let (a, b) be an even characteristic with a != c, coset
lattices L1, L2 as in section 1.

(i) Whether Theta_G[a, b] vanishes depends only on the plane span(a, c): the
translation tau -> tau + 1 sends (a, c) to (a, a + c), and the integral
symplectic matrix [[2I, 3G], [H, delta I]], with H = d G^{-1} and
delta = (3d + 1)/2, maps tau G to ((2 tau + 3)/(d tau + delta)) G and sends
(a, c) to (c, a). Together they act transitively on the six ordered bases.

(ii) If d is not a square, every element of Aut(L0) mapping L1 onto L2 can be
replaced by one interchanging them. So no zero has level 2, decomposable or
not. The key step is a lemma on order-four automorphisms g with g^2 != -I of
an even lattice M with v_2(det M) = 2: whenever g moves a class of M/2M in a
four-cycle, an explicit involution of M (built from the eigenlattices of g^2)
does the same. A square root of -I cannot occur, since it would make
det L0 = 4d a square.

(iii) If L is indecomposable of minus type, L1 and L2 are never isometric. If
Q(a) = 1, L0 is the set of integral-norm vectors of both L1 and L2, so an
isometry preserves L0; (ii) then gives an interchanging element, a symmetry
witness, which Theorem E excludes. If Q(a) = 0 then Q(c) = 1, and the dual
configuration rescaled by d (third lattice d G^{-1}, even, of minus type)
falls under the first case.

Consequently, any zero of an even indecomposable minus-type lattice would
have level 4: its coset lattices L1, L2 would be isospectral, non-isometric
Kneser 2-neighbours whose third intermediate lattice is even, indecomposable
and of minus type. The same holds in plus type at nonsquare determinant on
planes that are not totally singular. The minus-type parity conjecture is
thus equivalent to the absence of such isospectral pairs. A decomposable
control (A2 (+) [[2,1],[1,4]], det 21, a zero whose plane has exactly one
singular vector) shows that indecomposability must enter any proof; in
Theorems E and F it enters only through Theorem E.

**Proposition G (even lattices of odd determinant through 501).** Universe:
every isometry class of even positive definite rank-four lattices with odd
det <= 501, deduplicated by exact isometry: 1,890 classes, 1,696 of them
indecomposable, from four complete lists:

| Universe | Candidates | Classes | Indecomposable | Indecomposable classes with zeros |
|---|---:|---:|---:|---|
| odd det <= 201 (section 5, Task B) | 10,188 | 300 | 255 | det 121 (6 zeros) |
| odd det 202..279 (section 4) | 8,972 | 269 | 238 | none |
| det = 5 mod 8, 202..501 | 23,612 | 787 | 703 | none |
| det = 1 mod 8, 281..497 | 24,580 | 672 | 623 | det 361 (6 zeros) |

The second and third lists overlap in minus type for 202..279 (123
indecomposable classes), which are counted once in the totals. For every
indecomposable class, the set of zero characteristics equals the union over
all h in Aut(L) with h^2 = -I of the six ordered pairs on (h+I)L/2L, with
nothing unresolved. The only indecomposable classes with zeros are the
Gaussian lattice of determinant 121 and

    G = [[4,-2,0,-1],[-2,6,-1,0],[0,-1,4,2],[-1,0,2,6]],   det 361 = 19^2,

each with six zeros witnessed by a square root of -I. In particular all 834
indecomposable minus-type classes with det <= 501 have P, and in this range
every zero of an even indecomposable lattice of odd determinant has a
symmetry witness. The two new lists are independently replayed by
`code/verify_rank4_parity.py`.

Two constructions cannot produce a minus-type counterexample: the
Conway-Sloane common-kernel template never yields an even lattice of odd
determinant, and a bounded quarter-turn probe found only decomposable
lattices.

Reports: [minus-type report](../docs/rank-four-minus-type.md) (Theorem E, the
census through 501 in minus type, constructions),
[stage-two reductions](../docs/rank-four-stage-two-reductions.md) (Theorem F);
review: [rank-four-stage-two-review.md](../reviews/rank-four-stage-two-review.md)
(line-by-line check of Lemma A, Theorem E and Theorem F, an independent
computation over 31,500 configurations, and the plus-type list); preserved
outputs: `reports/rank4-minus-type/`, `reports/rank4-stage2/`,
`reports/rank4-plus-type/`.

## 7. The minus-type conjecture is proved: Theorem S (release 1.4.0)

**Theorem S.** Let H be a positive definite real quadratic form on Z^4 and
tau = I + N an involution of F2^4 with N^2 = 0 and rank N = 2. If for each
of the six pairs {k, tau k} of classes of Z^4/2Z^4 moved by tau the two
classes have equal minimal norm and equal minimal multiplicity, then some T
in GL4(Z) satisfies T^t H T = H and T = tau (mod 2).

**Corollary.** An even, orthogonally indecomposable, positive definite
rank-four lattice with det = 5 mod 8 has no identically vanishing
nontrivial even characteristic, at any determinant. In an even rank-four
lattice of odd determinant every such zero is a symmetry zero, so (Theorem E)
in an indecomposable lattice it comes from a square root of -I on a totally
singular plane, which occurs only in plus type.

The deduction: by plane invariance (Theorem F's setting) a zero makes all
three shifts of its plane P vanish; in the common lattice K = {y : y mod 2L
in P} each shift splits by norm modulo 8 into two halves, each equal to
theta(v + 2K) - theta(w + 2K) at every coefficient, so the six halves have
balanced minima in exactly the sense of Theorem S for H = B^t G B / 2 and a
square-zero rank-two marking tau; the lift T fixes 2L/2K, extends to an
automorphism of L, fixes P modulo 2L and reverses the character, hence is a
symmetry witness; Theorem E excludes it in indecomposable minus type.

Theorem S is computer-assisted. Voronoi's secondary fan of quaternary forms
has three orbits of maximal cones, all simplicial; the computed cones are
matched by explicit unimodular maps to the published Delta, Delta', Delta''.
On the relative interior of each of the 2,394 positive definite faces the
complete minimal sets of the fifteen classes are constant, and the six
minimum equalities are linear. Each of the 2,394 x 210 = 502,740 face and
marking cases has an exact disposition: 477,492 multiplicity mismatches (the
two minimizing sets), 19,992 lower-face exclusions (Stiemke vectors), and
5,256 cells with a rational interior point, an integral T = tau mod 2 and
D = Z E, which proves invariance on the whole equality subspace. None is
uncertified.

The classification was designed and run by Claude Opus 5.5 sessions,
replayed independently in standard-library Python by Grok Bot, audited
against the primary literature by Grok, and accepted in a final audit by a
Codex session (Astra), whose checker recomputes every disposition without the
classifier or its linear-programming backend. These are AI-assisted reviews
with independent computational replication, not peer review.

Report: [Theorem S](../docs/rank-four-theorem-s.md); audit:
[rank-four-theorem-s-final-audit.md](../reviews/rank-four-theorem-s-final-audit.md);
package: `reports/rank4-theorem-s/package/`; checks:
`code/rank4_theorem_s_check_published.py`, `code/rank4_theorem_s_check_bundle.py`.

## 8. Reporting conventions and corrections carried into this release

Every class count names its universe and is deduplicated by exact isometry.
The reviews record three precision points that this note already applies:
the plus- and minus-type counts are counts of distinct classes (134 and 131,
not 135 and 134, which repeated three census classes); witness orders quoted
from the stored certificates describe the first stored witness, not a
partition (the D4 zeros, for instance, also have involution witnesses); and
the decision counts are 26,775 decisions in the even list (26,769 nonzero)
and 120,120 in the odd list (120,065 nonzero, 55 zeros). Release 1.1.0 stated
the bracket 24 < det <= 1156 for the minimum counterexample determinant;
release 1.1.1 lowered it to 24 < det <= 280 (Theorem A'') and release 1.2.0
closed it (Theorem D). Release 1.3.0 adds section 6 and release 1.4.0 section 7. The archival reports and reviews keep their dated
correction notes.

## 9. Attribution

The isospectral pairs are Schiemann's (1990), Conway-Sloane's (1992) and
Earnest-Nipp's (1991); their non-isometry is Cervino-Hein's theorem (2011);
the theta-series determination of ternary forms is Schiemann's theorem (1997).
The [literature note](../docs/rank-four-isospectrality-literature.md) records
which pages were read and what each source proves. This project supplies the
Klein-configuration bridge, the compatible index-two re-embedding, the
characteristic zeros, the determinant census, the parity results and the
odd-determinant results of section 6 (the symplectic plane transformation is
an instance of the general theta transformation formula, DLMF 21.5.9). The
Eisenstein-integer structure of the counterexamples was observed in this
project and has not been checked against the literature. The weight-two
coefficient criterion used in the corollaries is the Sturm bound (see section
8.2 of Best et al., *Computing classical modular forms*).

## 10. Where to find things

| Item | Location |
|---|---|
| Reports | [level four](../docs/rank-four-level-four-zero.md), [Eisenstein family](../docs/rank-four-eisenstein-family.md), [minimum determinant](../docs/rank-four-minimum-determinant.md), [parity property](../docs/rank-four-parity-property.md), [minus type](../docs/rank-four-minus-type.md), [stage-two reductions](../docs/rank-four-stage-two-reductions.md), [Theorem S](../docs/rank-four-theorem-s.md), [literature](../docs/rank-four-isospectrality-literature.md) |
| Reviews | [level four](../reviews/rank-four-level-four-review.md), [Eisenstein](../reviews/rank-four-eisenstein-review.md), [parity](../reviews/rank-four-parity-review.md), [stage two](../reviews/rank-four-stage-two-review.md), [Theorem S final audit](../reviews/rank-four-theorem-s-final-audit.md), [assessment](../reviews/project-assessment-2026-10-21.md) |
| Builders and verifiers | `code/rank4_level4_search.py`, `code/verify_rank4_level4.py`, `code/verify_rank4_det280.py`, `code/rank4_eisenstein_search.py`, `code/verify_rank4_eisenstein.py`, `code/rank4_parity_search.py`, `code/verify_rank4_parity.py`, `code/rank4_gap_census.py`, `code/rank4_minus_type.py`, `code/rank4_minus_type_constructions.py`, `code/rank4_plus_type_census.py` |
| Independent review checks | `code/rank4_level4_review_check.py`, `code/rank4_eisenstein_review_check.py`, `code/rank4_eisenstein_review_enum.py`, `code/rank4_parity_review_check.py`, `code/rank4_parity_review_enum.py`, `code/check_small_pairs.py`, `code/rank4_stage2_checks.py`, `code/rank4_stage2_review_check.py`, `code/rank4_theorem_s_check_published.py`, `code/rank4_theorem_s_check_bundle.py` |
| Certificates | `data/rank4-level4/`, `data/rank4-eisenstein/`, `data/rank4-parity/` |
| Preserved outputs and replays | `reports/rank4-level4/`, `reports/rank4-eisenstein/`, `reports/rank4-parity/`, `reports/rank4-gap-census/`, `reports/rank4-minus-type/`, `reports/rank4-stage2/`, `reports/rank4-plus-type/`, `reports/rank4-theorem-s/package/`, `reports/project-review-2026-10-21/` |
| Tests | `tests/test_rank4_level4.py`, `tests/test_rank4_det280.py`, `tests/test_rank4_eisenstein.py`, `tests/test_rank4_parity.py`, `tests/test_rank4_gap_census.py`, `tests/test_rank4_stage_two.py` |
| Commands | [reproduction guide](../docs/reproduction.md) |
