# Review of the Eisenstein search (item 3)

Curated copy (release 1.1.0) of the review note RANK4-EISENSTEIN-REVIEW.md, written 2026-10-20 by the Claude Fable 5.1 session that wrote the handoff brief and amended 2026-10-21. File paths are relocated to the release layout.

Reviewed 2026-10-20 by the Claude Fable 5.1 session that wrote the handoff
brief (not distributed). The brief was executed by a second Claude Fable 5.1
session; its deliverables are [rank-four-eisenstein-family.md](../docs/rank-four-eisenstein-family.md),
`code/rank4_eisenstein_search.py`, `code/verify_rank4_eisenstein.py`,
`tests/test_rank4_eisenstein.py`, `data/rank4-eisenstein/det1156-level2.json`
and `reports/rank4-eisenstein/`. Nothing in the reviewed files was edited.

## 1. Verdict

Accepted. The headline result is correct and complete as stated:

    G = [[4,-2,1,-2],[-2,8,2,-2],[1,2,6,-4],[-2,-2,-4,12]],  a = (0,1,0,1),  b = (1,0,1,0)

is an even primitive positive definite lattice of determinant 1156 = 34^2
with Aut(L) = {I,-I}; Theta_G[a,b] (a dot b = 0, even) vanishes identically;
the sign kernel L0 is the trace form of the binary Hermitian Z[w]-lattice with
(p, q, delta) = (8, 6, 4 + 10w); an order-three isometry of L0 carries one
coset lattice onto the other and no automorphism of L0 interchanges them
(level 2). It is a counterexample to the rank-four symmetry converse of the
same kind as the determinant-23,716 example, with ambient determinant 17^2
instead of 77^2. It is the only level-2 example of the family F with
determinant below 1729, and the smallest in F up to D' = 154.

The bracket for the minimum counterexample determinant in rank four is now
24 < det(L) <= 1156 (census lower bound; this example). No minimality claim
is made between 25 and 1155. The determinant-1729 level-4 example of
[rank-four-level-four-zero.md](../docs/rank-four-level-four-zero.md) keeps its own record: it remains the smallest known
zero with non-isometric coset lattices.

**Correction, 2026-10-21 (determinant 280).** The determinant-4480 lattice of
the (2,8,14,20) configuration has every Gram entry even. Halving it gives the
odd primitive lattice G = [[3,1,1,0],[1,5,-1,2],[1,-1,4,-1],[0,2,-1,7]] of
determinant 280 with the same a = (1,0,0,1), b = (1,1,1,1). Scaling a Gram
matrix changes neither the vanishing nor the isometry type of the coset
lattices, so this is a level-4 zero; Aut(L) = {I,-I} with epsilon 0 on both
elements, so it is also a counterexample to the rank-four symmetry converse.
It is certified directly (zero by the modular cutoff, B = 3072, 383 complete
signed shells; group and non-isometry by complete searches;
`data/rank4-level4/det280-level4.json`, replayed by
`code/verify_rank4_det280.py`). The smallest known level-4 determinant is
therefore 280, not 1729, and the bracket for the minimum counterexample
determinant becomes 24 < det <= 280. No minimality claim is made.

## 2. Independent replication of the example

`code/rank4_eisenstein_review_check.py` (filed at the workspace root; output in
`reports/rank4-eisenstein/review-independent-check.txt`) uses only the
previously accepted modules `rank4_feasibility`, `exact_theta` and
`rank4_level4_search.analyse_kernel`, nothing from the new builder or
verifier. Results, all exact:

| Claim in [rank-four-eisenstein-family.md](../docs/rank-four-eisenstein-family.md) | Reviewer's result |
|---|---|
| det 1156, a dot b even, primitive, even | 1156; 0; gcd 1; even |
| modular data d0, N0, M, mu, B | 4624, 34, 544, 864, 144 |
| complete signed shells through B: 16 shells, 192 vectors, all zero | identical |
| (extra) signed shells through norm4 2000 | 36,512 vectors, all signed coefficients zero |
| Aut(L) = {I,-I} | 2 elements |
| G0 (project kernel basis) det 4624; Aut(L0) order 12 | identical |
| elements of Aut(L0) mapping v1 to v2: 2; swapping v1, v2 while fixing v3: 0 | 2 (orders 3 and 6, R0 and -R0); 0 |
| O_v1, O_v2, O_v3 of determinant 1156; O_v1 isometric to O_v2; O_v3 isometric to G | yes; 2 isometries; 2 isometries |
| trace form of (8, 6, 4+10w) is the sign kernel; D' = 34 | det 4624 = (3pq - N)^2, R-invariant, 12 isometries to G0; its line overlattice is isometric to G |
| `analyse_kernel` on that kernel (prefix 400) | 24 isospectral pairs: 6 at level 2 (two of the five lines), 18 at level 1, none at level 3 or 4 |

The direct call `rank4_counterexample.replay` on the embedded manuscript-schema
document returns `verified_rank4_symmetry_converse_counterexample`,
determinant 1156, full group order 2, all stabilizer epsilons zero.

## 3. Independent re-enumeration of the family

`code/rank4_eisenstein_review_enum.py` (output in
`reports/rank4-eisenstein/review-independent-enum.txt`) re-implements Task B
from the brief's formulas: reduced triples with p0 <= sqrt(D'),
N(delta0) <= p0^2, q0 = (2D' + N)/(3p0) >= p0; five line lifts; integrality
of the three overlattices by the direct test (G0 u even, u^t G0 u divisible
by 4); the line level by the complete `rank4_feasibility.isometries` group of
G0 (swap present or not); the primitive determinant of L = O_v3. It agrees
with the report on every count:

| Quantity | Report | Reviewer |
|---|---|---|
| reduced triples, 5 <= D' <= 41 | 498 | 498 |
| (triple, line) candidates / excluded by parity | 2490 / 1576 | 2490 / 1576 |
| members | 914 | 914 |
| members per D' | table in section 3 of the report | identical for all 24 values of D' |
| level-2 members | 12, all D' = 34 | 12, all D' = 34, all primitive of determinant 1156, each isometric to G (2 isometries), Aut(L0) order 12 |

So the completeness statement rests on the covering argument (which the
reviewer derived independently in the brief and the executing session
re-derived) and on two independent implementations of the scan.

## 4. Replays of the delivered material

Run from the workspace root with fresh output paths in the reviewer's
scratchpad:

| Command | Verdict | Time |
|---|---|---|
| `python test_rank4_eisenstein.py` | 13 tests OK | 5.7 s |
| `verify_rank4_eisenstein.py --enumeration task-b.json --certificates det1156-level2.json --task-e task-e.json` | verified: 914 members, 123 examples, 1 level 2 (1156); certificate `verified_level_2`, group 2, kernel group 12, three identity routes, cutoff 144 by box backend; Task E eight cases all level 1 | 22.2 s |
| `verify_rank4_eisenstein.py --enumeration task-d.json` | verified: 5539 triples, 9659 members, 1092 examples, 133 level 2, smallest 1156 | 224 s |

## 5. Reading of the code

- Enumeration (`reduced_triples`): the box |m0|, |n0| <= 2p0 contains every
  delta0 with N(delta0) <= p0^2 (N <= p0^2 forces |m0|, |n0| <= 2p0/sqrt 3);
  duplicates are harmless. D' = 2 mod 3 is impossible because N takes only
  the residues 0 and 1 mod 3, which matches the empty rows.
- Membership (`analyse_member`): the line is moved to Z[w]e1' by a basis
  change commuting with R, the invariant 3pq - N(delta) is rechecked, and the
  parity test is applied to the transformed parameters. The level comes from
  `analyse_kernel` (complete Aut(L0)) and is cross-checked against a second
  group search; the three pairs of a line must share one level in {1, 2}.
- Deduplication (`same_example`): exact isometry of the reduced primitive
  L's carrying (a, b) to (a', b') by primal and dual transport. The verifier
  recomputes every member, every isometry to its representative, and the
  pairwise non-equivalence of representatives of equal primitive determinant
  with `box_isometries`.
- Certificate: three identity routes. (i) R0 with R0^t G0 R0 = G0, R0^3 = I,
  R0 v1 = v2 is an isometry of L0 carrying O_v1 onto O_v2, which is a complete
  all-level proof. (ii) The cutoff is the project's `modular_data` and
  `cutoff_build`, hypotheses met (even G, b nonzero, a dot b even). (iii) The
  manuscript schema is built with M = G0/2 and V = 2K^{-1}, which is exactly
  the ambient (1/2)L0 of the manuscript; its replay is the unchanged
  `rank4_counterexample.replay`.
- Decomposition test: L is decomposable iff some involution T other than
  +-I in the full group has (I+T)/2 integral; both directions are correct
  (the projection onto the fixed sublattice splits L orthogonally because T
  is an isometry). With Aut(L) = {I,-I} the test is vacuous and the
  automorphism-count argument already applies.
- The prefix-60 incident of Task D is handled correctly. A differing coset
  prefix is an exact proof of non-vanishing at any prefix; equal prefixes are
  only a candidate filter, so the flaw could create spurious level-4
  candidates but never hide a zero or change a line level. The rerun at
  prefix 400 leaves no candidate, and the per-line results are identical in
  both runs. Keeping the flawed file as evidence is the right call.

## 6. Points to state more precisely (none affects the verdict)

- Section 3 of the report says the det-1156 class "contains 12 members (four
  reduced triples with p0 = 4, q0 = 6 and three lines each)". The stored
  table and the reviewer's enumeration both give six reduced triples
  (delta0 with N(delta0) = 4, that is the six units times 2) and two lines
  each. The count 12 is right; the breakdown is not.
- "The 154 non-primitive members have primitive determinants 1, 9, 16, 36,
  49, 81, 100 and are all level 1, in agreement with the census for det <= 24"
  is accurate only for the determinants 1, 9 and 16; the others are
  consistent with the level-1 finding of this scan itself, not with the
  census.
- "No level-3 or level-4 candidate appeared anywhere" has as its universe
  the 105 pairs over each of the 914 kernels of Task B (23,040 off-line
  isospectral pairs) and the 9659 kernels of Task D (150,480 pairs at prefix
  400); the report gives these numbers and should tie the sentence to them.
- The family F excludes level-2 configurations whose mapping automorphism
  moves the line; the report says so (section 6). The bracket 24 < det <= 1156
  is therefore a statement about all rank-four lattices, while "smallest in
  F" is the complete statement.
- The report's section 5 count "27,336 conjugate-linear and 23,388
  non-normalising swapping elements" is over all level-1 members with
  multiplicity, not over distinct examples; the structural conclusion (every
  level-1 line has a conjugate-linear swap) is the one worth quoting.

## 7. Consequences recorded by the reviewer

- [rank-four-level-four-review.md](rank-four-level-four-review.md) sections 1, 5 and 6 and [rank-four-isospectrality-literature.md](../docs/rank-four-isospectrality-literature.md)
  section 7 amended on 2026-10-20: the determinant record is 1156, not 1729;
  the level-4 statement about 1729 is unchanged.
- The release folder, the manuscript and the deferred-research
  tables are untouched. Whether the det-1156 lattice replaces the det-23,716
  lattice as the manuscript's headline example, with the Hermitian
  description added, is the user's decision at the next release update; the
  certificate already satisfies the manuscript's own schema, so the swap is
  mechanical.
- Items 1, 2 and 3 of the workspace to-do list (not distributed) are complete. Item 4 (parity
  conjecture from the census) was run next from a further
  handoff brief and is reviewed in
  [rank-four-parity-review.md](rank-four-parity-review.md) (the conjecture fails on both sides: det 49 for odd
  lattices, det 121 for even lattices). Item 5 (arXiv and authorship) was
  decided by the user on 2026-10-21: no arXiv submission, the release stays
  pseudonymous. Item 6 (root README pointer) was done the same day, and the
  release folder was updated to version 1.1.0 with the October work
  (`docs/rank-four-eisenstein-family.md`, `reviews/rank-four-eisenstein-review.md`).

## 8. Files added by this review

- [rank-four-eisenstein-review.md](rank-four-eisenstein-review.md) (this note).
- `code/rank4_eisenstein_review_check.py`, `code/rank4_eisenstein_review_enum.py`, and
  their outputs `reports/rank4-eisenstein/review-independent-check.txt`,
  `reports/rank4-eisenstein/review-independent-enum.txt`.
- The verifier replays were regenerated in the reviewer's scratchpad only;
  their verdicts equal those of `reports/rank4-eisenstein/replay-*.json`.
