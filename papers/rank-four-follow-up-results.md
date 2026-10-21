# Rank-four follow-up results: non-isometric coset lattices, a determinant-1156 counterexample, and the parity property

Release 1.1.0, 2026-10-21. This note states the results added after the
[counterexample manuscript](rank-four-counterexample.md), with pointers to the
full reports, the review notes, the exact certificates and the reproduction
commands. The computations were specified, run and reviewed by Claude Fable
5.1 sessions and then assessed by a separate Claude session
([reviews](../reviews/README.md)); they are AI-assisted local reviews with
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
below exhibit level 4 in rank four, a smaller level-2 example, and settle two
parity conjectures negatively.

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

## 3. A smaller symmetry counterexample: determinant 1156

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
and recovers the determinant-23,716 example. Bracket for the smallest
determinant of a rank-four counterexample to the symmetry converse:

    24 < det(L) <= 1156,

the lower bound from the determinant-24 census, the upper bound from this
example. No minimality claim is made between 25 and 1155; the family minimum
is not a global minimum.

Full report: [rank-four-eisenstein-family.md](../docs/rank-four-eisenstein-family.md);
review: [rank-four-eisenstein-review.md](../reviews/rank-four-eisenstein-review.md);
certificate: [det1156-level2.json](../data/rank4-eisenstein/det1156-level2.json).

## 4. The parity property in rank four

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
[deferred research](../docs/deferred-research.md). Order-4 witnesses acting
faithfully on L/2L, order-8 witnesses, and zeros without any witness are not
excluded by the theorems.

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

## 5. Reporting conventions and corrections carried into this release

Every class count names its universe and is deduplicated by exact isometry.
The reviews record three precision points that this note already applies:
the plus- and minus-type counts are counts of distinct classes (134 and 131,
not 135 and 134, which repeated three census classes); witness orders quoted
from the stored certificates describe the first stored witness, not a
partition (the D4 zeros, for instance, also have involution witnesses); and
the decision counts are 26,775 decisions in the even list (26,769 nonzero)
and 120,120 in the odd list (120,065 nonzero, 55 zeros).

## 6. Attribution

The isospectral pairs are Schiemann's (1990), Conway-Sloane's (1992) and
Earnest-Nipp's (1991); their non-isometry is Cervino-Hein's theorem (2011);
the theta-series determination of ternary forms is Schiemann's theorem (1997).
The [literature note](../docs/rank-four-isospectrality-literature.md) records
which pages were read and what each source proves. This project supplies the
Klein-configuration bridge, the compatible index-two re-embedding, the
characteristic zeros and the parity results. The Eisenstein-integer structure
of the counterexamples was observed in this project and has not been checked
against the literature. The weight-two coefficient criterion used in the
corollaries is the Sturm bound (see section 8.2 of Best et al., *Computing
classical modular forms*).

## 7. Where to find things

| Item | Location |
|---|---|
| Reports | [level four](../docs/rank-four-level-four-zero.md), [Eisenstein family](../docs/rank-four-eisenstein-family.md), [parity property](../docs/rank-four-parity-property.md), [literature](../docs/rank-four-isospectrality-literature.md) |
| Reviews | [level four](../reviews/rank-four-level-four-review.md), [Eisenstein](../reviews/rank-four-eisenstein-review.md), [parity](../reviews/rank-four-parity-review.md), [assessment](../reviews/project-assessment-2026-10-21.md) |
| Builders and verifiers | `code/rank4_level4_search.py`, `code/verify_rank4_level4.py`, `code/rank4_eisenstein_search.py`, `code/verify_rank4_eisenstein.py`, `code/rank4_parity_search.py`, `code/verify_rank4_parity.py` |
| Independent review checks | `code/rank4_level4_review_check.py`, `code/rank4_eisenstein_review_check.py`, `code/rank4_eisenstein_review_enum.py`, `code/rank4_parity_review_check.py`, `code/rank4_parity_review_enum.py`, `code/check_small_pairs.py` |
| Certificates | `data/rank4-level4/`, `data/rank4-eisenstein/`, `data/rank4-parity/` |
| Preserved outputs and replays | `reports/rank4-level4/`, `reports/rank4-eisenstein/`, `reports/rank4-parity/`, `reports/project-review-2026-10-21/` |
| Tests | `tests/test_rank4_level4.py`, `tests/test_rank4_eisenstein.py`, `tests/test_rank4_parity.py` |
| Commands | [reproduction guide](../docs/reproduction.md) |
