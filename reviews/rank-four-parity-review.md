# Review of the parity-property search (item 4)

Curated copy (release 1.1.0) of the review note RANK4-PARITY-REVIEW.md, written 2026-10-20 by the Claude Fable 5.1 session that wrote the handoff brief and amended 2026-10-21. File paths are relocated to the release layout.

Reviewed 2026-10-20 by the Claude Fable 5.1 session that wrote the handoff
brief (not distributed). The brief was executed by a second Claude Fable 5.1
session; its deliverables are [rank-four-parity-property.md](../docs/rank-four-parity-property.md),
`code/rank4_parity_search.py`, `code/verify_rank4_parity.py`, `tests/test_rank4_parity.py`,
`data/rank4-parity/` and `reports/rank4-parity/`. Nothing in the reviewed
files was edited. The control `code/rank4_parity_facts_check.py` and its output
`reports/rank4-parity/control-facts.txt` still have the hashes recorded when
the brief was written (da9b8dce..., 2724808e...).

## 1. Verdict

Accepted. Both conjectures of the brief are refuted by complete class lists,
and the structural theorems of the report are correct as stated.

- (C1) "even, indecomposable, rank four, odd determinant implies the parity
  property" fails first at determinant 121:

      G = [[4,-2,-1,0],[-2,4,0,1],[-1,0,4,-2],[0,1,-2,4]],

  the trace form of the Hermitian Z[i]-lattice with Gram [[4, 2+i],[2-i, 4]].
  It has Aut(L) of order 24, no splitting involution, one Eichler component,
  and exactly six vanishing even characteristics, the six ordered pairs of
  distinct nonzero classes of the totally singular plane (1+i)L/2L; each is a
  symmetry zero whose witnesses are the two elements g with g^2 = -I. It is
  the only failure among the 255 indecomposable classes of even lattices with
  odd determinant at most 201 (300 classes, 45 decomposable; universe and
  dedup criterion in section 3).
- (C0) for odd lattices fails first at determinant 49:

      G = [[3,-1,-1,0],[-1,3,0,-1],[-1,0,3,1],[0,-1,1,3]],

  Aut(L) of order 16, indecomposable, seven symmetry zeros with involution
  witnesses. It is the smallest of 13 failing classes among the 1144
  indecomposable classes of odd lattices with 25 <= det <= 168 and 4 not
  dividing det (3487 classes in all); the census covers det <= 24 and has no
  such failure. The brief's bracket [25, 169] is replaced by the exact value.
- Theorems 1 to 4 of the report (lines b = Ga mod 2 never vanish; coset
  lattices of odd-det zeros are not integral; no involution witness in an
  even indecomposable rank-four lattice of odd determinant; a Z[i]-structure
  forces six zeros on its fixed plane) were checked line by line (section 6).
  Theorem 3 is the one genuinely new argument; it uses all three hypotheses
  and is correct. Theorem 4 is the rank-four form of the E8 mechanism.

What is proved and what is not, in the reviewer's words: for even
indecomposable rank-four lattices of odd determinant, every symmetry zero has
a witness whose 2-part has order 4 or 8, a square root of -I always produces
zeros (plus type only), and no other mechanism was seen through det 201.
Minus type (det = 5 mod 8) has no zero in the data and no proof. The parity
property is therefore not a determinant condition in rank four; the surviving
open statement is the one the report names (no zeros in minus type, and no
zeros outside Z[i]-structures in plus type).

## 2. Independent replication of the headline lattices

`code/rank4_parity_review_check.py` (filed at the workspace root; output in
`reports/rank4-parity/review-independent-check.txt`) uses only the accepted
modules `rank4_feasibility` (isometries, modular_data, cutoff_build,
determinant) and `exact_theta` (complete_shells, symmetry_data), with the
reviewer's own stabilizer, decomposability and F2 helpers; it imports nothing
from the new builder or verifier. Results, all exact (17 s in total):

| Claim in [rank-four-parity-property.md](../docs/rank-four-parity-property.md) | Reviewer's result |
|---|---|
| det 121, even, primitive, plus type (det = 1 mod 8) | 121; even; gcd 1; 9 nonzero singular classes |
| Aut(L) of order 24, element orders 1:1, 2:13, 3:2, 4:2, 6:2, 12:4 | identical; 2 elements with g^2 = -I |
| plane (g+I)L/2L = {(0,1,1,0), (1,0,0,1), (1,1,1,1)}, totally singular | identical; Q = 0 on all three, B = 0 pairwise |
| 105 characteristics: 99 nonzero, 6 zeros, all with order-4 witnesses | identical; the six zeros are exactly the six ordered pairs (a, c) of the plane; epsilon(g) = 1 on all six; max first nonzero norm4 12 |
| stabilizer of (a,b): four elements, epsilon 1 on the two of order 4 | identical |
| cutoff d0 484, N0 22, M 352, mu 576, B 96; 21 shells, 260 vectors, all zero | identical (extra: all signed coefficients zero through norm4 2000, 112,272 vectors) |
| G0 = [[4,-4,1,0],[-4,16,-8,2],[1,-8,8,-3],[0,2,-3,4]], Aut(L0) of order 8, 2 mappers, 2 swappers | identical |
| indecomposable by both criteria | 0 splitting involutions; 1 Eichler component (6 indecomposable vectors up to sign) |
| trace form of [[4, 2+i],[2-i, 4]] is isometric to G | Re h on (e1, i e1, e2, i e2) has det 121 and 24 isometries to G; multiplication by i preserves it and squares to -I |
| det-49 lattice: odd, Aut 16, indecomposable, 7 zeros with involution witnesses | identical (the zero a = b = (1,1,1,1) also has an order-8 witness); the quoted witness has epsilon 1 |
| family G_k, k = 4, 5, 7, 8: det 70, 91, 133, 154, Aut 4, fail | identical; 4 zeros each, involution witnesses |
| 18 pair members: 17 hold, CS remark-(x) G1 (det 961) decomposable and fails | identical; G1 has 2 splitting involutions and 36 zeros; largest first-nonzero norm4 among holders 74 (Earnest-Nipp) |
| ten even pair members are plus type with 36 totally singular ordered pairs | identical |
| eight level-4 lattices: 104 nonzero, 1 zero (the certified one); only det 4480 even | identical (the certified characteristic is the only one undecided through norm4 400; the zero itself rests on the item-2 certificates) |
| E8: g^t G g = G, g^2 = -I, fixed 4-space totally singular, 120 of 210 ordered pairs witnessed by g | identical |
| the three decomposable classes with g^2 = -I (det 9, 49, 121) fail | identical (36 zeros each) |

## 3. Independent re-enumeration of the class lists

`code/rank4_parity_review_enum.py` (outputs `reports/rank4-parity/review-independent-enum-{B201,C64,E168}.txt`)
re-derives the three universes with a different box: a Minkowski-reduced
basis (G11 <= G22 <= G33 <= G44; Q(x) >= G_kk for every x in {0,+-1}^4 with
x_k = 1 and x_l = 0 for l > k; G11 G22 G33 G44 <= 4 det, the rank-four
Minkowski constant; G12, G13, G14 >= 0 by sign changes), not the project's
successive-minima box. Dedup is exact isometry (complete basis-image search of
`rank4_feasibility.isometries`) between candidates of equal determinant and
equal vector counts through norm 40. Decisions use the reviewer's ladder
(symmetry witness; complete shells through norm4 100, 400, 1600). Control
runs with the constant doubled (`--lam 8`) give the same candidate sets, so
the constant is never binding at these sizes.

| Universe | Report | Reviewer |
|---|---|---|
| B: even, odd det <= 201 | 10,188 candidates, 300 classes, 45 decomposable, 255 indecomposable, 254 hold, 1 fails (det 121) | 5,489 candidates, 300 classes, 45 / 255 / 254 / 1 (det 121, Aut 24, 6 zeros, order-4 witnesses); classes per determinant identical for all 50 determinants; Aut orders 2:1, 4:54, 8:68, 12:46, 16:17, 24:39, 48:25, 72:2, 144:1, 240:2 identical; 124 plus, 131 minus; type = det mod 8 rule on every class; exactly four classes with g^2 = -I (det 9, 49, 121 decomposable, 121 indecomposable), all failing; first-nonzero distribution over the 26,769 nonzero decisions identical, maximum 132 at det 197 (Aut 24); 0 undecided |
| B: det <= 49 and <= 101 | 22 and 79 classes | the per-determinant counts above restrict to 22 and 79 |
| C: even, 4 \| det <= 64 | 1,101 candidates, 69 classes, 42 decomposable, 27 indecomposable, 7 hold, 20 fail | 876 candidates, 69 classes, 42 / 27 / 7 / 20, same failing determinants, per-determinant counts identical; 122 zeros, every one with an involution witness, 18 also with order-8 witnesses |
| E: odd lattices, 25 <= det <= 168, 4 not dividing det | 46,494 candidates, 3,487 classes, 2,343 decomposable, 1,144 indecomposable, 1,131 hold, 13 fail | 20,862 candidates, 3,487 classes, 2,343 / 1,144 / 1,131 / 13 at the same determinants 49, 70, 91, 91, 119, 126, 130, 130, 133, 133, 154, 154, 161; classes per determinant and indecomposable classes per determinant identical for all 108 determinants; 55 zeros, every one with an involution witness; twelve failing groups of order 4 and the det-49 group of order 16; largest first-nonzero norm4 142; 0 undecided (731 s) |
| E control to det 100 with the constant doubled | 1,197 classes (sum of the report's per-determinant counts through 99) | 7,104 candidates (same as with the constant 4), 1,197 classes, per-determinant counts identical |

The first-nonzero distribution of a class is an isometry invariant (an
isometry carries characteristics to characteristics with equal theta
functions), so its agreement across two different representatives per class
is a real check and not a coincidence of bases.

## 4. Replays of the delivered material

Run from the workspace root on 2026-10-20 with fresh output paths in the
reviewer's scratchpad (the delivered replays in `reports/rank4-parity/` were
not overwritten), while the re-enumerations ran on other cores:

| Command | Verdict | Time |
|---|---|---|
| `python scripts/run_tests.py   # runs tests/test_rank4_parity.py with the other rank-four tests` | 11 tests OK | 38 s |
| `rank4_parity_search.py checks` | unrestricted generator equals the frozen 1,510 det-24 candidates; 6 even odd-det classes in bijection with the census | - |
| `verify_rank4_parity.py certificates c1-counterexamples.json` | verified; 6 certificates `verified_c1_counterexample`, det 121, level 1, cutoff 96, witness order 4 | < 1 s |
| `verify_rank4_parity.py pair-members pair-members-parity.json` | verified; 17 hold, det-961 G1 fails | 2.5 s |
| `verify_rank4_parity.py task-a task-a.json` | verified; 169 census classes, 27 indecomposable; Eisenstein smallest level 1 det 169, level 2 det 2401; 17 + 1 pair members; 8 level-4 lattices fail | 20 s |
| `verify_rank4_parity.py classes task-b-det49.json` / `-det101.json` / `-det201.json` | verified; 22 / 79 / 300 classes, 16 / 65 / 255 indecomposable, 0 / 0 / 1 failing (det 121) | 13 s / 40 s / 109 s |
| `verify_rank4_parity.py classes task-c-det64.json` | verified; 69 classes, 27 indecomposable, 20 failing | 25 s |
| `verify_rank4_parity.py classes task-e-det168.json --fast-box` | verified; 46,494 candidates regenerated, 3,487 classes, 69,578 pairwise box searches, 1,144 indecomposable, 1,131 hold, 13 fail at det 49, 70, 91, 91, 119, 126, 130, 130, 133, 133, 154, 154, 161 | 394 s |

The `--fast-box` option only selects the outer-box enumeration backend of the
accepted `verify_rank4_level4` for every bound; it does not loosen any check.

## 5. Reading of the code

- Generator (`candidate_grams`): the project's accepted box with the bound as
  a parameter; even-only restricts the diagonal step to 2 (an even lattice has
  even diagonals in every basis, an odd lattice in none, since x^t G x = sum
  G_ii x_i^2 mod 2), and the determinant filter is applied at the last
  diagonal. At D0 = 24 it reproduces `rank4_feasibility.candidates(24)`
  exactly. The verifier regenerates every candidate list with a loose-box
  scanner and requires equality of the sorted lists.
- Dedup (`class_list`): a differing invariant (determinant, vector counts
  through norm 12) is an exact proof of non-isometry; equal invariants go to
  the complete basis-image search, which raises on interruption. The replay
  rechecks every stored transport and the pairwise non-isometry of
  representatives of equal determinant with the box search.
- Decisions (`decide_characteristic`): the stabilizer and epsilon follow the
  accepted convention of `exact_theta.symmetry_data` (T on a, T^{-t} on b;
  epsilon = a dot (T^{-t} b - b)/2 mod 2), and epsilon is a homomorphism on
  the stabilizer, so an epsilon-one element is a complete proof of vanishing.
  A nonzero complete-shell coefficient is a complete proof of non-vanishing.
  The ladder exhausted leads to the modular cutoff (accepted `modular_data`,
  `cutoff_build`) and otherwise to `unresolved`, never to a verdict. The
  verifier recomputes the whole coefficient histogram of every nonzero
  decision through the stored bound with its own box enumerations, and the
  epsilon of every stored witness with `symmetry_data`.
- Decomposability: the splitting-involution criterion and the Eichler
  component criterion are both correct as argued in the docstrings (a vector
  of norm <= max G_ii is a sum of indecomposable vectors of no larger norm;
  indecomposable vectors of an orthogonal summand lie in it), and the code
  requires them to agree on every lattice. The verifier recomputes both.
- Task E universe: the all-parity box with 4 not dividing det, even Gram
  matrices discarded before dedup. That is the complete list of odd lattices
  in the range, since a lattice is even iff some basis has even diagonals.
- Task A: census verdicts are taken from the stored pair records of the
  accepted census (every zero rechecked with `symmetry_data`); Eisenstein and
  level-4 inputs are hash-bound; the five level-4 zeros whose cutoffs exceed
  the vector budget are marked `zero_certified_externally` with the hash and
  level of the item-2 certificate, and the verifier rechecks that binding.
- Tests: eleven adversarial tests (removed or added classes, dropped
  candidates, wrong bounds, wrong coefficients, tampered groups, flipped
  decomposability, wrong summaries, witness removed, node limits) all fail
  the replay as they should.

## 6. The proofs

Checked by the reviewer in full; the notation is that of the report.

- Theorem 1. For y = a + 2x in an even lattice, y.y = a.a + 4 B(a,x) mod 8,
  and b = G a mod 2 makes the sign (-1)^{x.b} = (-1)^{B(a,x)}, so the signed
  coefficient of every shell is plus or minus its vector count. Correct, any
  rank, no determinant condition.
- Theorem 2. L0 + Z a/2 is integral iff a.x is even for all x in L0 = c^perp,
  i.e. a in {0, c}. Correct; the doubled coset lattices are even of
  determinant 256 det L.
- Theorem 3. The eigenlattice decomposition a = a_+ + a_-, c = c_+ + c_- is
  valid because g fixes a and c mod 2L (odd determinant converts the dual
  condition g^t b = b mod 2 into g c = c mod 2L); epsilon(g) = a_-.c_- mod 2
  and a.c even force a_+.c_+ odd, which excludes rank 0 and rank 1 for both
  eigenlattices (rank 1 gives even products in an even lattice), so both have
  rank two. The glue group is killed by 2, embeds in both discriminant groups,
  has order 2^k with k >= 1 by indecomposability, and det L_+ det L_- =
  4^k det L with det L odd forces the 2-part of each to be exactly 2^k, so
  k <= 2. k = 1 contradicts det = 0 or 3 mod 4 for an even binary form;
  k = 2 makes the 2-part of L_+^*/L_+ equal to (Z/2)^2, so both Smith
  invariants of Gram(L_+) are even and a_+.c_+ is even. Correct; each of
  "even", "odd determinant" and "indecomposable" is used exactly where the
  report says.
- Corollary. epsilon is a homomorphism of the stabilizer (checked above), so
  the 2-part of a witness is a witness, and the finite orders in GL_4(Z) leave
  2-parts of order 4 or 8. Correct.
- Theorem 4. (g+I)(g-I) = -2I gives 2L <= F = (g+I)L <= L of index
  |det(g+I)| = 4; g^t G = -G g gives (g+I)x.(g+I)y = 2 x.y, so F/2L is a
  totally singular plane (plus type); g fixes every class of F mod 2L;
  epsilon(g) = B(alpha, c) for a = (g+I) alpha; the functional is well
  defined on L/F' (B((g-I)x, (g+I)y) = -2 x.gy), nonzero, and vanishes at
  gamma with c = (g+I) gamma because gamma.g gamma = 0; so exactly the two
  classes a not in {0, c} have epsilon 1. Correct. The same argument in rank
  2n gives a totally singular fixed n-space on which g witnesses 2^{n-1} of
  the 2^n - 2 partners of each class, which is the 8 of 14 the report checks
  on E8. Corollaries (i) and (ii) follow; (iii) is a remark about parameters
  (det = D^2 with D = 3 mod 4 for even Z[i]-lattices of odd determinant).
- Section 5.3 quotes two standard facts without proof (all 1575 totally
  singular planes of E8/2E8 form one orbit, and W(E8) maps onto O_8^+(2)); the
  computed statements (35 planes and 210 ordered pairs in the fixed 4-space,
  120 witnessed by g) are replicated.

## 7. Points to state more precisely (none affects the verdict)

- Section 2.3 and 5.2 count "135 plus-type" and "134 minus-type" classes by
  adding the census classes (det 17; det 5, 13, 21) to the Task B list, but
  Task B already contains every class with det <= 201, so the distinct
  counts are 134 and 131. The statements made about them are unaffected.
- Section 4 says that of the 122 Task C zeros "104 have involution witnesses,
  18 have witnesses of order 8", and section 2.1 says the nine D4 zeros have
  witnesses of order 8. All 122 have involution witnesses (the reviewer's
  enumeration lists witness orders 2, 4 and 8 for D4 and 2D4); the split is
  by the stored (first-found) witness, not a partition of the zeros. The
  structural statement that matters, that every degenerate zero and every
  odd-lattice zero has an involution witness while no odd-det even zero can
  (Theorem 3), is unaffected.
- Section 6: the det-49 lattice's zero a = b = (1,1,1,1) has witnesses of
  orders 2 and 8; "all 55 zeros ... with involution witnesses" is right. The
  same section says the minimal shell decides "107,389 of the 120,120 nonzero
  decisions"; 120,120 is the number of decisions (1144 x 105), of which
  120,065 are nonzero and 55 are zeros.
- Observation not in the report (data only, no claim): the det-49 odd lattice
  also has two automorphisms g with g^2 = -I; their fixed plane
  {(0,1,1,0), (1,0,0,1), (1,1,1,1)} is isotropic for the bilinear form mod 2,
  and g has epsilon 1 on four of its six ordered pairs; the other two, and the
  seventh zero a = b = (1,1,1,1), are witnessed by involutions, which Theorem
  3 forbids for even lattices and Theorem 1 forbids for a = c in even
  lattices. So the Z[i] mechanism is present on the odd side too, in a form
  that does not need the quadratic refinement.
- Section 2.2 lists "level 2, det 2401" as the smallest odd Eisenstein
  failure at level 2; this is a statement inside the family F (Task D
  representatives), not about all lattices, and the report says so.
- Section 1 says "none does through det 201" for mechanisms other than
  Z[i]-structures; the universe is the 254 holding indecomposable classes of
  Task B (no zero at all) plus the one failing class (all six zeros on the
  plane). Nothing is said about even lattices of odd determinant above 201
  beyond the ten pair members of section 2.3.

## 8. Consequences recorded by the reviewer

- The handoff brief (not distributed) carries an outcome note at the top
  from this review: (C1) is refuted at det 121, and the bracket
  [25, 169] for the smallest odd failure is replaced by 49.
- [rank-four-eisenstein-review.md](rank-four-eisenstein-review.md) section 7 amended to mark item 4 complete. No
  sentence in [rank-four-eisenstein-family.md](../docs/rank-four-eisenstein-family.md), [rank-four-eisenstein-review.md](rank-four-eisenstein-review.md) or
  [rank-four-isospectrality-literature.md](../docs/rank-four-isospectrality-literature.md) presents det 169 as the smallest odd failure outside
  the family, so nothing else needed correction (report section 7, item 3).
- Item 4 of the workspace to-do list (not distributed) is the user's text
  and is not edited; its
  conjecture is false on both sides (det 49 odd, det 121 even), and its
  remark that the det-23,716 example rules out an involution-glue criterion
  is unaffected.
- The release folder, the manuscript and the deferred-research
  tables are untouched; porting is the user's decision at the next release
  update.
- 2026-10-21: the user decided against item 5 (no arXiv submission; the
  release stays pseudonymous on GitHub and Zenodo) and item 6 was done (a
  status pointer at the top of the workspace `README.md`).
- 2026-10-21: the external assessment [project-assessment-2026-10-21.md](project-assessment-2026-10-21.md) was
  incorporated. Its corollary that the ordinary theta series does not
  determine the parity property in rank four (the det-961 remark-(x) pair:
  equal theta series by a weight-two Sturm bound on Gamma_0(124), not
  isometric, one member fails P and the other holds) was added to report
  section 2.3 with its modular proof; the assessment's script
  `code/check_small_pairs.py` was rerun by this
  session with identical output. The precision points of section 7 above
  (distinct class counts 134 and 131, witness orders by stored witness,
  120,120 decisions with 120,065 nonzero) were applied to the report as
  dated corrections.
- 2026-10-21: ported to the release folder as version 1.1.0
  (`docs/rank-four-parity-property.md`, `reviews/rank-four-parity-review.md`,
  code, tests, certificates and compressed class lists). The new research
  directions of the assessment (minus-type conjecture, plus-type companion
  question, level-3 existence, global minimum, human review) were added to
  the workspace to-do list as items 7 to 11 and to
  [deferred-research.md](../docs/deferred-research.md).

## 9. Files added by this review

- [rank-four-parity-review.md](rank-four-parity-review.md) (this note).
- `code/rank4_parity_review_check.py` and
  `reports/rank4-parity/review-independent-check.txt`.
- `code/rank4_parity_review_enum.py` and
  `reports/rank4-parity/review-independent-enum-B201.txt`,
  `review-independent-enum-C64.txt`, `review-independent-enum-E168.txt`.
- The verifier replays and the two box-constant controls were written to the
  reviewer's scratchpad only; their verdicts equal those of
  `reports/rank4-parity/replay-*.json`.
