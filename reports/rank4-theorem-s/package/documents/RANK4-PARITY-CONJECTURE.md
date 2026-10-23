# The rank-four parity property: (C1) refuted at determinant 121, the odd side settled at 49

Executed 2026-10-20 by a Claude Fable 5.1 session from
`Fable-rank4-parity-conjecture.md` (item 4 of `A_what_is_next.txt`), for
review by the brief-writing session. All new files are listed in section 8;
no existing file was modified. Python 3.11 standard library, exact
arithmetic throughout. Every count below names its universe and its
deduplication criterion (exact isometry: complete basis-image search), and
every negative result names what was excluded.

Notation. L is a positive definite integral rank-four lattice with Gram G.
The parity property P(L): every even characteristic (a, b binary columns,
a and b nonzero mod 2, a.b even) has Theta_L[a,b] not identically zero.
For even L write V = L/2L, Q(x) = x^t G x / 2 mod 2, B(x,y) = x^t G y mod 2,
Rad = ker(G mod 2); for odd det, c = G^{-1} b mod 2 is the primal class of b,
the sign kernel L0 is c^perp, and a.b even means B(a,c) = 0. "Level" is the
hierarchy of `LITERATURE-RANK4.md` section 1 (1 = symmetry witness, 2 = an
isometry of the sign kernel carrying one coset lattice onto the other without
interchanging them, 3 = an isometry not preserving the sign kernel, 4 =
isospectral non-isometric).

## 1. Verdict

**(C1) is false.** "L even, indecomposable, rank four, det(L) odd implies
P(L)" fails first at determinant 121. The complete list of isometry classes of
even positive definite rank-four lattices with odd determinant at most 201
(Task B: 10,188 reduced-basis candidates, 300 classes, 255 indecomposable)
contains exactly one indecomposable class without P:

    G = [[4,-2,-1,0],[-2,4,0,1],[-1,0,4,-2],[0,1,-2,4]],   det 121 = 11^2,

the trace form of the Hermitian Z[i]-lattice with Gram [[4, 2+i],[2-i, 4]]
(exact isometry checked). It has Aut(L) of order 24 and exactly six even
zeros, all of level 1: they are the six ordered pairs (a, c) of distinct
nonzero classes of the plane Pi = (1+i)L/2L, which is totally singular for Q,
and every witness is an element g of order 4 with g^2 = -I (multiplication by
i). The other 254 indecomposable classes, including all 131 of minus type
(det = 5 mod 8), hold P with every one of their 105 characteristics proved
nonzero by a complete shell coefficient; the 123 other plus-type classes all
contain 36 ordered pairs spanning totally singular planes, none a zero. The
structural content (section 5): in an even indecomposable rank-four lattice of
odd determinant no involution is a symmetry witness (Theorem 3, proved), so
every symmetry zero needs an element whose 2-part has order 4 or 8; a Z[i]
structure g (g^2 = -I) always supplies six zeros on its totally singular fixed
plane (Theorem 4, proved), which is the rank-four analogue of the E8 witness;
and b = G a mod 2 never gives a zero in any even lattice (Theorem 1, proved).
What remains open is whether symmetry zeros with witnesses of other shapes, or
identity zeros without witnesses, can occur for odd determinant; none does
through det 201, and no minus-type class has a zero.

**The odd-lattice side.** (C0) "indecomposable and 4 does not divide det
implies P" fails first at determinant 49, not 169: Task E lists every odd
rank-four lattice with 25 <= det <= 168 and 4 not dividing det (46,494
candidates, 3487 classes, 1144 indecomposable), of which 13 fail P, the
smallest being G = [[3,-1,-1,0],[-1,3,0,-1],[-1,0,3,1],[0,-1,1,3]] (det 49,
Aut of order 16, seven level-1 zeros with involution witnesses). All 55 zeros
of the 13 failing classes are level 1 with involution witnesses, which is
exactly the mechanism Theorem 3 excludes for even lattices.

## 2. The evidence table (Task A)

Control: `rank4_parity_facts_check.py` was rerun unchanged before the work and
its three sections are reproduced by `python rank4_parity_search.py A`
(`results/rank4-parity/task-a.json`, rendering `task-a.md`): 169 census
classes, 27 indecomposable, holding tally {even det=1 mod 4: 4, odd det=1: 3,
odd det=3: 5, odd det=2: 2, odd det=0: 1}, failing tally {even det=0 mod 4: 4,
odd det=0: 8}; 123 Task B Eisenstein representatives, 51 indecomposable, the 12
indecomposable odd lattices with 4 not dividing det listed by the control
(det 169 the smallest); both Conway-Sloane det-1729 lattices hold P.

**Decomposability.** Two criteria are applied to every lattice and required
to agree (they do on all 169 + 1215 + 26 lattices of the table and on every
class of Tasks B, C, E). (i) Splitting involution: L is decomposable iff some
involution T != +-I in Aut(L) has (I+T)/2 integral. Proof: if L = A (+) B then
T = (-1 on A, +1 on B) is such an involution; conversely P = (I+T)/2 is an
integral idempotent, self-adjoint for the form because T is, so L = PL (+)
(I-P)L orthogonally with both parts nonzero. The code constructs the two
summands and checks the direct sum (`decomposition_from_involution`). (ii)
Eichler components: vectors of norm at most max G_ii that are not a sum of two
nonzero orthogonal vectors, joined when their inner product is nonzero; the
components generate mutually orthogonal sublattices whose sum is L (every basis
vector is a sum of indecomposable vectors of no larger norm), so at least two
components means decomposable, and one component means indecomposable
(`eichler_components`). All 45 decomposable classes of Task B split as 2+2
(an even rank-four lattice of odd det has no rank-one even summand).

**2.1 Census (universe: the 169 isometry classes with det <= 24 of the
accepted census, exact-isometry dedup).** Even classes: 18, namely 12 with
4 | det and 6 with odd det (det 5, 9, 13, 17, 21, 21); det mod 4 of an even
class is 0 or 1 on all 18 (the classical fact, checked). The even odd-det
classes are of minus type for det 5, 13, 21, 21 and plus type for 9, 17. The
27 indecomposable rows are in `task-a.md`; P status and zero counts come from
the stored exact verdicts (every zero has a stored epsilon-one witness, every
nonzero a complete shell; the census stores all 105 even characteristics with
a, b nonzero for every class). Stored (first-found) witness orders of the even
zeros: order 8 for the nine zeros of D4 (det 4), order 2 for all other even
indecomposable census zeros; the D4 zeros also have involution witnesses
(precision note of 2026-10-21, `RANK4-PARITY-REVIEW.md` section 7). No census
class has b = Ga mod 2 at a zero.

**2.2 Eisenstein family F (universes: the 123 Task B and the 1092 Task D
representatives of `results/rank4-eisenstein/`, exact-isometry dedup there).**
All fail P (each carries its line zero). Task B: 51 indecomposable, 26 even
(all with 4 | det), 25 odd lattices (13 with 4 | det, 12 with det = 1 mod 4).
Task D: 726 indecomposable, 356 even (all 4 | det), 370 odd lattices (190
with 4 | det, 180 with det = 1 mod 4). Smallest indecomposable odd lattice
with 4 not dividing det, by level of the line zero: level 1, det 169 (D' = 13,
Aut of order 8, as in the control); level 2, det 2401 = 49^2 (D' = 49, Aut of
order 2, G = [[5,2,1,-1],[2,9,3,-1],[1,3,7,2],[-1,-1,2,11]], a = b = (1,1,1,1));
no level-2 example among the 123 Task B representatives. Both are now known
not to be smallest overall (section 6: det 49). The F2 position of the 468 + 48
even line zeros is recorded for Task D (section 5.1).

**2.3 Isospectral-pair members and level-4 lattices (universe: the 18 members
of the nine pairs of `results/rank4-level4/task-a.json` and the eight
certified level-4 lattices of `examples/rank4-level4/`).** P was decided for
all 26 lattices, 105 characteristics each, exactly.

| lattice | det | parity | type | Aut | P | zeros |
|---|---:|---|---|---:|---|---|
| CS (1,7,13,19) L^+, L^- | 1729 | even | plus | 2, 2 | holds, holds | 0 |
| CS remark (x) (11,1), G1 | 961 | odd | - | 4 | fails (decomposable) | 36 |
| CS remark (x) (11,1), G2 | 961 | odd | - | 2 | holds | 0 |
| Earnest-Nipp G1, G2 | 84721 | even | plus | 2, 2 | holds, holds | 0 |
| CS (1,7,13,43) G1, G2 | 3913 | even | plus | 2, 2 | holds, holds | 0 |
| CS (1,7,25,31) G1, G2 | 5425 | even | plus | 2, 2 | holds, holds | 0 |
| CS (1,7,13,31) G1, G2 | 2821 | odd | - | 2, 2 | holds, holds | 0 |
| CS (1,7,19,25) G1, G2 | 3325 | odd | - | 2, 2 | holds, holds | 0 |
| CS (2,8,14,20) G1, G2 | 4480 | odd | - | 2, 2 | holds, holds | 0 |
| CS (7,13,19,49) G1, G2 | 84721 | even | plus | 2, 2 | holds, holds | 0 |
| eight level-4 lattices | 1729, 2821, 3325, 3913, 4480, 5425, 84721, 84721 | odd except 4480 | - | 2 | fail | 1 each, level 4 |

Seventeen pair members hold P (the brief expected eighteen); the exception is
the first lattice of the Conway-Sloane remark-(x) pair, which is decomposable
(splitting involution found, two Eichler components) and fails by PROOFS.md
section 7. The ten even pair members are all of plus type with 36 totally
singular ordered pairs each and no zero; the largest first-nonzero index among
the 17 holding members is norm4 = 74 (Earnest-Nipp). Each level-4 lattice has
exactly one even zero, its certified level-4 characteristic, and 104 nonzeros;
the zero is re-proved here by the inclusive cutoff for det 1729 (B = 17920),
2821 (B = 28672) and 4480 (B = 6144), and for the other five (cutoffs 38400 to
125440, above the 5,000,000-vector limit of this run) it is taken from the
item-2 certificates (hash-bound; replayed in `results/rank4-level4/`).
Certificates of the 18 pair members: `examples/rank4-parity/pair-members-parity.json`.

**Corollary (added 2026-10-21 from `PROJECT-REVIEW-2026-10-21.md`): the
ordinary theta series does not determine the parity property in rank four.**
The two members G1, G2 of the Conway-Sloane remark-(x) pair of determinant
961 have the same theta series at every level: 2G1 and 2G2 are even of level
124 and square determinant 124^2, so their theta series are weight-two forms
on Gamma_0(124) with trivial character; the index is 192, the inclusive Sturm
bound 32, and the coefficients 0 through 32 agree by two independent
enumerations. They are not isometric (two independent searches). Yet G1 is
decomposable and fails P with 36 zeros, while G2 is indecomposable and holds
P. So the ordinary representation numbers cannot recover the parity property.
The statement concerns rank-four lattices in general; it does not settle the
same question restricted to pairs of indecomposable lattices. Script and
saved output: `results/project-review-20261021/check_small_pairs.py`,
`small-pairs.json` (rerun by the reviewing session, identical output).

## 3. Task B: all even classes with odd det <= 201

**Generator.** `candidate_grams(D0, even_only, det_filter)` restates the
reduced Gram box of `RANK4-FEASIBILITY.md` section 1 (text in the module
constant `COVERAGE`): a basis with Gram-Schmidt coefficients |mu_ij| <= 1/2
whose projected vectors are successively shortest satisfies d1 = A, d_(i+1) >=
(3/4) d_i, det >= (3/4)^6 A^4 and G_ii >= A; cross entries lie in integral
intervals of radius d_j/2 around sum_(l<j) mu_jl mu_kl d_l; diagonals satisfy
G_kk >= max(A, q + (3/4) d_(k-1)) and prefix * d_k^r * (3/4)^(r(r-1)/2) <= D0
for the r remaining pivots; the fourth diagonal solves prefix*(G44 - q) = det.
The even restriction enters only as "all diagonals even", which discards bases
of odd lattices and nothing else, since every basis of an even lattice has even
diagonals; the determinant filter discards bases outside the target universe.
Checks (`python rank4_parity_search.py checks`): with no restriction at D0 = 24
the generator returns exactly the frozen 1,510 candidate bases of
`rank4_feasibility.candidates(24)`; with the even/odd-det restriction at 24 it
returns 168 candidates reducing to 6 classes in exact-isometry bijection with
the census's even odd-det classes (det 5, 9, 13, 17, 21, 21, decomposable
A2 (+) A2 and A2 (+) [[2,1],[1,4]] included). The independent scanner of the
verifier (loose boxes and direct diagonal scanning, generalised from
`verify_rank4_feasibility.independent_candidates`) regenerates every list.

**Deduplication.** Candidates with equal determinant and equal vector counts
for all norms <= 12 are compared by the complete basis-image isometry search
(`rank4_feasibility.isometries`); a differing invariant is an exact proof of
non-isometry. Each class list stores every candidate with a unimodular
transport to its representative; the replay regenerates the candidates, checks
all transports and proves pairwise non-isometry of representatives of equal
determinant with the independent box search.

**Cost and counts (universe: isometry classes of even positive definite
rank-four lattices with odd det <= D0).**

| D0 | candidates | classes | decomposable | indecomposable | P holds | P fails | build time | replay time |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 49 | 636 | 22 | 6 | 16 | 16 | 0 | 9.5 s | 12 s |
| 101 | 2,576 | 79 | 14 | 65 | 65 | 0 | 27.6 s | 36 s |
| 201 | 10,188 | 300 | 45 | 255 | 254 | 1 (det 121) | 103 s | 109 s |

All three lists completed with every isometry search finished (node limit
5,000,000 never reached); 201 is the largest bound tried, not a limit of the
method. Classes per determinant at 201 range from 1 (det 5, 9, 13, 17, 25, 49,
169) to 19 (det 189). Aut(L) orders of the 255 indecomposable classes: 2 (1),
4 (54), 8 (68), 12 (46), 16 (17), 24 (39), 48 (25), 72 (2), 144 (1), 240 (2).
Types: 124 plus (det = 1 mod 8), 131 minus (det = 5 mod 8), as the classical
rule predicts on every class.

**Decisions.** 255 x 105 = 26,775 exact decisions: 26,769 nonzero by a
complete shell coefficient, 6 zeros by symmetry (det 121), none unresolved, no
cutoff needed. First-nonzero norm4 distribution: 2 (4,235), 4 (4,204), 6
(4,310), 8 (3,517), 10 (2,575), 12 (1,859), 14 (1,274), 16 (988), 18 (698), 20
(553), 22 (386), 24 (308), 26 (220), 28 (236), 30 (192), 32 (168), 34 (144),
36 (134), 38 (84), 40 (88), 42 (64), 44 (88), 46 (84), 48 (68), 50 (64), 52
(52), 54 (28), 58 (8), 60 (20), 62 (12), then 12 each at 68, 76, 84, 92, 100,
108, 116, 124, 132. The largest, 132, belongs to the minus-type class
[[2,-1,0,-1],[-1,2,0,0],[0,0,2,-1],[-1,0,-1,34]] of det 197 (Aut of order 24).
The minimal shell of a + 2L decides 24,203 of the 26,769 nonzero cases (26,775 decisions; corrected 2026-10-21)
(90.4%); 48 decisions needed the second ladder step (norm4 400), none the
third. The classes carrying an automorphism g with g^2 = -I are exactly four:
det 9 (A2 (+) A2), det 49 ([[2,1],[1,4]] (+) [[2,1],[1,4]]), det 121
([[2,1],[1,6]] (+) [[2,1],[1,6]]) and the indecomposable det-121 class; all
four fail P, as Theorem 4 predicts.

**The counterexample (certificate `examples/rank4-parity/c1-counterexamples.json`,
six documents, one per zero; replay `results/rank4-parity/replay-certificates.json`).**
Reduced primitive Gram G as in section 1 (LLL-reduced from the list
representative [[4,-2,-1,-1],[-2,4,0,1],[-1,0,4,2],[-1,1,2,4]] by
U = [[1,0,0,0],[0,1,0,0],[0,0,1,-1],[0,0,0,1]]); first zero a = (0,1,1,0),
b = (0,1,1,0), c = (1,0,0,1); Q(a) = Q(c) = Q(a+c) = 0, B(a,c) = 0, type plus
(9 nonzero singular classes). Sign kernel basis K = [[1,0,0,0],[0,2,-1,0],
[0,0,1,0],[0,0,0,1]], G0 = [[4,-4,1,0],[-4,16,-8,2],[1,-8,8,-3],[0,2,-3,4]]
(d0 = 484), complement t = e2, classes v1 = (0,1,1,0), v2 = (0,0,1,0) in L0
coordinates; coset lattices O_v1, O_v2 with 4 x Gram [[16,-6,4,0],[-6,8,0,-2],
[4,0,32,-12],[0,-2,-12,16]] and [[16,-16,2,0],[-16,64,-16,8],[2,-16,8,-6],
[0,8,-6,16]], LLL-reduced to [[8,2,-2,-2],[2,12,-2,2],[-2,-2,16,4],[-2,2,4,24]]
and [[8,2,2,2],[2,12,2,-2],[2,2,16,4],[2,-2,4,24]] (det 30976 each, isometric).
Identity proofs: (i) inclusive cutoff N0 = 22, M = 352, B = 96, complete
signed shells through norm4 96 (21 shells, 260 vectors), all zero, replayed by
the box backend; (ii) symmetry witness g = [[0,0,0,-1],[0,0,-1,0],[0,1,0,0],
[1,0,0,0]] of order 4 with g^2 = -I and epsilon = 1; (iii) the swap in
Aut(L0) (order 8; two elements map v1 to v2, both interchange them). Full
Aut(L) of order 24 (element orders 1:1, 2:13, 3:2, 4:2, 6:2, 12:4); stabilizer
of (a,b): four elements, epsilon 0 on +-I and 1 on the two elements of order 4.
Indecomposable by both criteria (six indecomposable vectors up to sign in the
norm-4 sphere, one component). Level 1. The other five certificates are the
other ordered pairs of Pi = {(0,1,1,0), (1,0,0,1), (1,1,1,1)}.

## 4. Task C: even classes with 4 | det <= 64 (data for Task D)

Universe: isometry classes of even positive definite rank-four lattices with
4 | det <= 64; 1,101 candidates, 69 classes (42 decomposable, 27
indecomposable), built in 15.8 s, replayed in 25 s. Of the 27 indecomposable
classes 7 hold P (det 28, 32, 48, 52, 52, 52, 64) and 20 fail (det 4, 12, 20,
20, 28, 32, 36, 36, 36, 44, 44, 48, 52, 52, 56, 60, 60, 60, 64, 64). Their 122
zeros are all level 1 and all have involution witnesses; the stored first-found
witness is an involution for 104 of them and has order 8 for the other 18 (the
nine zeros of D4 and of 2D4). (Corrected 2026-10-21; the original text read
"104 have involution witnesses, 18 have witnesses of order 8".) Position
relative to the radical of the
form mod 2 (every indecomposable class here has radical dimension 2 except 2D4,
whose G is 0 mod 2):

| zeros | rad dim | Q on rad | a in rad | b in im(G mod 2) | Q(a) | Q(c) | span(a,c) totally singular |
|---:|---:|---|---|---|---:|---|---|
| 38 | 2 | not zero | yes | yes | 1 | 0 or 1 | no |
| 26 + 9 | 2 | zero | yes | yes | 0 | 1 | no (includes D4) |
| 18 | 2 | zero | yes | yes | 0 | 0 | yes |
| 12 | 2 | zero | no | no | 1 | undefined | - |
| 6 + 4 | 2 | zero / not zero | yes / no | no | 0 | undefined | - |
| 9 | 4 | zero | yes | no | 0 | undefined | - (2D4) |

So in the degenerate case a lies in the radical for 112 of the 122 zeros,
and the E8-type condition fails outright (D4: Q(c) = 1) or is undefined (b
not in the image of G mod 2) for 104 of them.

## 5. Task D: structure

### 5.1 Is span(a, c) totally singular at every known even zero? (D1)

Universe: the 1,031 known even zeros of even rank-four lattices in this
workspace (census even classes 386, Eisenstein Task B 48 and Task D 468 line
zeros, the level-4 det-4480 zero, Task B 6, Task C 122), 538 of them in
indecomposable lattices. For the 493 zeros of decomposable lattices the
condition carries no information: the 36 zeros of A2 (+) A2 (det 9, plus type)
all span totally singular planes, the 36 zeros of A2 (+) [[2,1],[1,4]] (det 21,
minus type) cannot, since minus type has no totally singular plane, and the
other 421 decomposable zeros have 4 | det. Among the 538 indecomposable zeros:

- the 6 zeros of odd determinant (det 121, plus type) all span totally
  singular planes;
- the 532 zeros with 4 | det (radical of dimension 2, or 4 for det-4480 and
  2D4) split into 203 with a in Rad, b in the image, Q(a) = Q(c) = 0 (totally
  singular modulo the radical), 94 with a in Rad and Q(a) = 1 or Q(c) = 1
  (the condition fails; D4 is among them), and 235 with b not in the image of
  G mod 2, where c does not exist (203 Eisenstein line zeros, 31 Task C
  zeros and the det-4480 level-4 zero, whose Gram is twice an odd lattice).

**Answer.** In the nondegenerate case (odd det) the E8-type necessary
condition holds at every known zero, and section 5.3 proves it for the only
known mechanism. In the degenerate case it is dead: D4 itself violates it, and
the dominant feature of degenerate zeros is a in Rad (a/2 in L*), with 235
zeros having no c at all. No zero anywhere has b = Ga mod 2 (Theorem 1).

### 5.2 Plus type: totally singular pairs are almost never zeros (D2)

In a nondegenerate plus-type quadratic space of dimension 4 over F2 there are
9 nonzero singular vectors, 6 totally singular planes (the two rulings of the
hyperbolic quadric) and 36 ordered pairs (a, c) spanning them; minus type has
5 singular vectors and no totally singular plane. The data confirm 36 on every
one of the 134 distinct plus-type indecomposable classes of the table (Task B
124, which already contain the census class of det 17, and the 10 even pair
members) and 0 on the 131 distinct minus-type classes (all in Task B, which
contains the census classes of det 5, 13, 21). (The counts 135 and 134 of the
original text added the census classes a second time; corrected 2026-10-21.)
Zeros occur on exactly one plus-type class (det 121), on
one of its six planes.

**Theorem 1 (lines are never zeros).** Let L be an even lattice of any rank
and (a, b) an even characteristic with b = G a mod 2 L^* (for odd det: a = c).
Then every shell of a + 2L has constant sign and Theta_L[a,b] is not zero.
Proof. For y = a + 2x, y.y = a.a + 4 a.x + 4 x.x = a.a + 4 B(a,x) mod 8
because x.x is even; and x^t b = x^t G a mod 2 = B(a,x). So the sign
(-1)^{x.b} on the shell of norm N equals (-1)^{(N - a.a)/4}, the coefficient
is +-(number of vectors), and the first nonempty shell is nonzero.

**Theorem 2 (coset lattices of odd determinant are not integral).** For odd
det and a != c, neither coset lattice L0 + Z a/2 nor L0 + Z(a/2 + t) is
integral: integrality of L0 + Z a/2 needs a.x even for all x in L0 = c^perp,
i.e. a in (c^perp)^perp = {0, c}. So an even zero of an even lattice of odd
determinant is an isospectral pair of non-integral lattices of odd determinant
det L whose doubles are even lattices of determinant 256 det L. This records
the brief's suggestion; it blocks a direct appeal to classification results on
integral lattices and yields no contradiction by itself.

**Theorem 3 (no involution witness for odd determinant).** Let L be even,
rank four, indecomposable, with odd determinant. No involution g in Aut(L) has
epsilon(g) = 1 on an even characteristic with a, b nonzero mod 2.
Proof. Let g be an involution stabilizing (a, b) with epsilon(g) = 1, and let
L_+ = L cap V_+, L_- = L cap V_- be the eigenlattices. From g a = a mod 2L,
a_- = -(ga - a)/2 lies in L, hence a_+ = a - a_- lies in L_+ and a_- in L_-;
likewise c = c_+ + c_- with c_+ in L_+, c_- in L_- (odd det identifies the
dual condition with g c = c mod 2L). epsilon(g) = B((ga - a)/2, c) =
a_-.c_- mod 2, so a_-.c_- is odd; since a.c = a_+.c_+ + a_-.c_- is even,
a_+.c_+ is odd too. Hence both eigenlattices are nonzero, and neither has rank
one (in Z v, alpha v . gamma v = alpha gamma v.v is even because L is even).
So both have rank two. Put k with [L : L_+ (+) L_-] = 2^k; k >= 1 because L is
indecomposable, the glue group Gamma = L / (L_+ (+) L_-) is killed by 2
(2x = (x + gx) + (x - gx)), it embeds into L_+^*/L_+ and into L_-^*/L_-, and
det(L_+) det(L_-) = 4^k det L. As det L is odd, the 2-parts of det L_+ and
det L_- multiply to 4^k while each is at least 2^k (the image of Gamma), so
each equals 2^k, and k <= 2 because a rank-two discriminant group has at most
two generators. k = 1: det L_+ = 2 mod 4, impossible for an even binary form
[[2u, v],[v, 2w]] (its determinant is 0 or 3 mod 4). k = 2: the 2-part of
L_+^*/L_+ has order 4 and contains (Z/2)^2, so it is (Z/2)^2, both Smith
invariants of the Gram of L_+ are even, every inner product in L_+ is even,
and a_+.c_+ is even, a contradiction.
Corollary. epsilon(g^m) = m epsilon(g), so for odd m the 2-part of any witness
is a witness; every symmetry witness in such a lattice has 2-part of order 4
or 8 (orders of finite elements of GL_4(Z) are 1,2,3,4,5,6,8,10,12). The data
agree: the six det-121 witnesses have order 4; the 55 odd-lattice zeros of
Task E and the 104 + 18 degenerate zeros of Task C have witnesses of order 2
or 8, i.e. the hypotheses "even" and "odd det" are both used.

**Theorem 4 (Z[i]-structures give six zeros).** Let L be even, rank four,
odd det, and g in Aut(L) with g^2 = -I. Then F = (g+I)L satisfies 2L <= F <= L
with [L:F] = |det(g+I)| = 4, all inner products on F are even ((g+I)x.(g+I)y =
x.(g^t G + G g + 2G) y = 2 x.y because g^t G = G g^{-1} = -G g), so Pi = F/2L
is a totally singular plane and L has plus type. Every ordered pair (a, c) of
distinct nonzero classes of Pi is a symmetry zero with witness g.
Proof. F = {x : gx = x mod 2L} (if (g-I)x in 2L then x = -(g+I)(g-I)x/2 in
(g+I)L; conversely (g-I)(g+I) = -2I). Write a = (g+I)alpha, c = (g+I)gamma.
Then g fixes a and c mod 2L, B(a,c) = 0, and epsilon(g) = B((ga-a)/2, c) =
B(-alpha, c) = alpha.c mod 2. The functional alpha -> B(alpha, c) is well
defined on L/F' with F' = (g-I)L (B(F', F) = 0 because (g-I)x.(g+I)y = -2x.gy),
(g+I) induces an isomorphism L/F' -> Pi, the functional is nonzero (Rad = 0),
and it vanishes at gamma (gamma.(g+I)gamma = gamma.g gamma + gamma.gamma, and
gamma.g gamma = 0 since g^{-1} = -g). So its kernel is {0, gamma}, and the two
classes with B(alpha, c) = 1 are exactly those with a not in {0, c}: both give
epsilon(g) = 1.
Corollaries. (i) A minus-type even rank-four lattice of odd determinant has no
automorphism squaring to -I. (ii) In an indecomposable such lattice an element
g of order 4 whose square acts trivially on L/2L has g^2 = -I (an involution
congruent to I mod 2L splits L into its eigenlattices), so Theorem 4 covers
every order-4 witness with that property. (iii) For an even Z[i]-lattice with
Hermitian Gram [[p, delta],[conj delta, q]] the trace form has det
(pq - N(delta))^2; evenness needs p, q even, odd det needs N(delta) odd, so
det = D^2 with D = 3 mod 4: the classes with such a g through det 201 are
exactly the four listed in section 3 (D = 3, 7, 11, 11).

**The minimal-shell route and the obstruction.** The minimal shell of a + 2L
decides 90.4% of the 26,769 nonzero decisions (26,775 decisions in all) of Task B but not all (the
first nonzero index reaches norm4 = 132 at det 197), so no uniform minimal-shell
argument is available for odd determinant. What is proved: Theorems 1, 3, 4
and the corollaries. What is open: symmetry zeros whose witness has 2-part of
order 4 acting with order 4 on L/2L, or of order 8, and identity zeros without
any witness (levels 2 to 4). The smallest case where the argument does not
close is any class whose group contains such an element; none of the 254
holding classes through det 201 has a zero, so the obstruction is one of
proof, not of data.

### 5.3 Why E8 escapes, and the rank-four analogue (D3)

E8 = D8^+ carries the coordinate rotation g = diag(J, J, J, J), J = [[0,-1],
[1,0]], which preserves D8 and the glue vector class, so g in Aut(E8) with
g^2 = -I (`python rank4_parity_search.py e8`, `results/rank4-parity/e8-check.json`,
checked exactly: g^t G g = G, g^2 = -I in the integral basis of E8 produced
from the standard generators). Its fixed space F/2E8 = (1+i)E8/2E8 has
dimension 4, is totally singular, contains 35 of the 1,575 totally singular
planes, and all 210 ordered pairs of distinct nonzero classes in it are E8
zeros by the plane criterion of `E7-E8-MILESTONE.md`; g itself is the
epsilon-one witness for 120 of them (eight of the fourteen partners of each
class, because in rank eight the functional alpha -> B(alpha, c) on L/F' has a
kernel of dimension 3). The group of order 696,729,600 supplies a conjugate of
g through every plane. The rank-four analogue is exactly a Z[i]-structure
whose fixed plane is the plane itself (Theorem 4); nothing prevents it for odd
determinant in plus type, and the smallest indecomposable instance is det 121.
In minus type the mechanism is impossible (Corollary (i)); for 4 | det the
witnesses are involutions and the fixed-plane picture is replaced by a in Rad.

### 5.4 Minus type (D4)

For det = 5 mod 8 there is no totally singular plane, no automorphism squaring
to -I (Theorem 4) and no involution witness (Theorem 3). This does not prove
P: witnesses with 2-part of order 4 acting faithfully of order 4 on L/2L
(O_4^-(2) is isomorphic to S_5 and contains such elements) or of order 8 are
not excluded, and identity zeros are not excluded. The plane condition is not
known to be necessary for identity zeros. Data: all 131 distinct minus-type
indecomposable classes through det 201 hold P with every characteristic
decided by a complete shell. (The original text counted 134 by adding the
three census classes that Task B already contains; corrected 2026-10-21.)

## 6. Task E: the odd-lattice side

Known odd indecomposable failures with 4 not dividing det before this task:
det 169 at level 1 and det 2401 at level 2 (Eisenstein), det 1729 at level 4.
Cost estimate (`python rank4_parity_search.py estimate`): all-parity
candidates with 4 not dividing det at D0 = 24, 48, 72, 96 number 812, 3,447,
8,334, 15,687 (generation 0.04 to 1.3 s), so the odd-lattice listing to 168
was estimated at about 50,000 candidates and well under an hour; it was run.

Universe: isometry classes of odd (not even) positive definite rank-four
lattices with 25 <= det <= 168 and 4 not dividing det; even Grams produced by
the box are discarded before dedup (an even lattice has no odd basis). 46,494
candidates, 3,487 classes (det mod 4 = 1: 1,066; 2: 1,216; 3: 1,205), class
list 134.5 s, decisions to 453 s. 2,343 classes are decomposable (they fail P
by PROOFS.md section 7) and 1,144 indecomposable: 1,131 hold P and 13 fail,
with determinants 49, 70, 91, 91, 119, 126, 130, 130, 133, 133, 154, 154, 161
(det mod 4 = 1: 4, 2: 6, 3: 3). All 55 zeros of the 13 failing classes are
level 1 with involution witnesses; twelve of the thirteen groups have order 4
(three involutions), the det-49 group order 16. First-nonzero index at most
norm4 = 142; the minimal shell decides 107,389 of the 120,120 decisions
(1144 x 105; 120,065 nonzero and 55 zeros; corrected 2026-10-21); 32 needed
the second ladder step; no cutoff was needed.

**The true smallest odd indecomposable failure with 4 not dividing det is
det 49:**

    G = [[3,-1,-1,0],[-1,3,0,-1],[-1,0,3,1],[0,-1,1,3]],  Aut of order 16,

seven even zeros, e.g. a = b = (0,1,1,0) with involution witness
[[0,0,0,-1],[0,-1,0,0],[0,0,1,0],[-1,0,0,0]]. All 42 indecomposable odd
classes with 25 <= det < 49 and 4 not dividing det hold P. The family
G_k = [[3,-1,-1,0],[-1,3,0,-1],[-1,0,3,1],[0,-1,1,k]] (det 7(3k-2): 49, 70,
91, 133, 154 for k = 3..5, 7, 8) accounts for five of the thirteen failures.

## 7. Corrections needed elsewhere (not applied)

1. `Fable-rank4-parity-conjecture.md` section 1: "(C1) ... the even case is
   clean so far" is refuted at det 121; "the smallest odd indecomposable
   failure with 4 not dividing det lies in [25, 169]" is now exactly 49.
2. `A_what_is_next.txt` item 4: the conjecture "indecomposable with
   determinant not divisible by 4 implies the parity property" is false for
   odd lattices (det 49) and for even lattices (det 121); its remark that the
   det-23716 counterexample "rules out any involution-glue iff criterion" is
   unaffected.
3. `RANK4-EISENSTEIN-SEARCH.md` and its review: if any sentence presents the
   det-169 member as the smallest indecomposable odd lattice without P (rather
   than the smallest inside the family F), it should now cite det 49.
4. No other document found by search states (C0) or (C1).

## 8. Deliverables and regeneration

New files: `RANK4-PARITY-CONJECTURE.md` (this report), `rank4_parity_search.py`,
`verify_rank4_parity.py`, `test_rank4_parity.py`,
`examples/rank4-parity/c1-counterexamples.json` (six certificates),
`examples/rank4-parity/pair-members-parity.json`, and in `results/rank4-parity/`:
`task-a.json`, `task-a.md`, `task-a.log`, `task-b-det49.json`,
`task-b-det101.json`, `task-b-det201.json` (+ `.log`), `task-c-det64.json`
(+ `.log`), `task-d.json` (+ `.log`), `task-e-estimate.json` (+ `.log`),
`task-e-det168.json` (+ `.log`), `e8-check.json` (+ `.log`), `certify.log`,
`replay-task-a.json`, `replay-task-b-det49.json`, `replay-task-b-det101.json`,
`replay-task-b-det201.json`, `replay-task-c-det64.json`,
`replay-task-e-det168.json`, `replay-certificates.json`,
`replay-pair-members.json` (each with `.log`), `tests.txt` (11 tests, all
pass). `control-facts.txt` is untouched. Output paths must not exist before a
command runs (the builders refuse to overwrite). The decided class lists store
every candidate, every full group, every shell table and the F2 data of every
decision, so they are large: `task-b-det201.json` is about 121 MB and
`task-e-det168.json` about 459 MB; parse them, do not paste them into a model
context. Both are regenerated by the commands below in under eight minutes.

```powershell
python rank4_parity_facts_check.py                       # control, unchanged
python rank4_parity_search.py checks
python rank4_parity_search.py A --output results/rank4-parity/task-a.json --markdown results/rank4-parity/task-a.md --progress
python rank4_parity_search.py B --dmax 49  --output results/rank4-parity/task-b-det49.json --progress
python rank4_parity_search.py B --dmax 101 --output results/rank4-parity/task-b-det101.json --progress
python rank4_parity_search.py B --dmax 201 --output results/rank4-parity/task-b-det201.json --progress
python rank4_parity_search.py C --dmax 64  --output results/rank4-parity/task-c-det64.json --progress
python rank4_parity_search.py estimate --output results/rank4-parity/task-e-estimate.json
python rank4_parity_search.py E --dmin 25 --dmax 168 --output results/rank4-parity/task-e-det168.json --progress
python rank4_parity_search.py certify --task-b results/rank4-parity/task-b-det201.json --output examples/rank4-parity/c1-counterexamples.json
python rank4_parity_search.py D --task-a results/rank4-parity/task-a.json --task-b results/rank4-parity/task-b-det201.json --task-c results/rank4-parity/task-c-det64.json --output results/rank4-parity/task-d.json
python rank4_parity_search.py e8 --output results/rank4-parity/e8-check.json
python verify_rank4_parity.py classes results/rank4-parity/task-b-det201.json --output results/rank4-parity/replay-task-b-det201.json
python verify_rank4_parity.py classes results/rank4-parity/task-e-det168.json --fast-box --output results/rank4-parity/replay-task-e-det168.json
python verify_rank4_parity.py task-a results/rank4-parity/task-a.json --output results/rank4-parity/replay-task-a.json
python verify_rank4_parity.py certificates examples/rank4-parity/c1-counterexamples.json --output results/rank4-parity/replay-certificates.json
python verify_rank4_parity.py pair-members examples/rank4-parity/pair-members-parity.json --output results/rank4-parity/replay-pair-members.json
python -m unittest -v test_rank4_parity
```

The pair-member file `examples/rank4-parity/pair-members-parity.json` is the
`level4.pair_members` block of `task-a.json` (same documents, hash-bound to
`results/rank4-level4/task-a.json`). Replay verdicts: Task A, Task B (49, 101,
201), Task C, the six certificates and the pair members all `verified`
(box backends: rectangular box of `verify_lattice_census` for bounds <= 2000,
the outer-box enumeration of `verify_rank4_level4` for the large cutoffs and,
with `--fast-box`, for the Task E list). The Task E replay (46,494 candidates
regenerated, 3,487 classes, 69,578 pairwise box isometry searches, 1,144
decided classes, 13 failures at the same determinants) is `verified` in 407 s.
All replays completed without reaching a node limit.
