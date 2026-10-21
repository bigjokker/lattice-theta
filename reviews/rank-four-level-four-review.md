# Review of the level-4 search (item 2)

Curated copy (release 1.1.0) of the review note RANK4-LEVEL4-REVIEW.md, written 2026-10-20 by the Claude Fable 5.1 session that wrote the handoff brief and amended 2026-10-21. File paths are relocated to the release layout.

Reviewed 2026-10-20 by the Claude Fable 5.1 session that wrote the handoff
brief (not distributed). The brief was executed by a second Claude Fable 5.1
session; its deliverables are [rank-four-level-four-zero.md](../docs/rank-four-level-four-zero.md),
`code/rank4_level4_search.py`, `code/verify_rank4_level4.py`, `tests/test_rank4_level4.py`,
`data/rank4-level4/` and `reports/rank4-level4/`. This note records what
the reviewer checked, with which tools, and the verdict. Nothing in the
reviewed files was edited.

## 1. Verdict

Accepted. Question Q of the brief has a positive answer, and the central
example is correct as stated:

    G = [[3,1,1,0],[1,9,-2,4],[1,-2,7,-2],[0,4,-2,13]],  a = (1,0,0,1),  b = (1,1,1,1)

is a primitive odd positive definite lattice of determinant 1729 = 7*13*19
with Aut(L) = Aut(L0) = {I,-I}; Theta_G[a,b] (an even characteristic, a dot b
= 2) vanishes identically; its two coset lattices L0 + Z(a/2) and
L0 + Z((a+2e1)/2) are isometric to Conway-Sloane's L^-(1,7,13,19) and
L^+(1,7,13,19) (Schiemann's pair), which are not isometric. That is level 4 of
the hierarchy in [rank-four-isospectrality-literature.md](../docs/rank-four-isospectrality-literature.md) section 1: the vanishing has no
geometric explanation, inside L0 or outside it. Since Aut(L) = {I,-I} and
epsilon vanishes on both elements, the lattice is also a second counterexample
to the rank-four symmetry converse, of smaller determinant than the
determinant-23,716 example (1729 against 23,716), and it is orthogonally
indecomposable (a decomposable lattice has at least four automorphisms).
(Later the same day the Eisenstein sweep of item 3, [rank-four-eisenstein-family.md](../docs/rank-four-eisenstein-family.md),
reviewed in [rank-four-eisenstein-review.md](rank-four-eisenstein-review.md), gave a level-2 counterexample of
determinant 1156; the determinant record is now 1156 and the level-4
statement about 1729 is unchanged.)

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
`code/verify_rank4_det280.py`). The smallest level-4 determinant is therefore
280, not 1729, and the bracket for the minimum counterexample determinant
becomes 24 < det <= 280. Later the same day the determinant census through 279
(`docs/rank-four-minimum-determinant.md`, Theorem D of the follow-up note)
showed that every rank-four lattice with det <= 279 satisfies the converse, so
280 is the minimum.

The one new computational fact behind the result is that the index-two
sublattices {all e-coordinates even} of L^+ and of L^- are isometric although
L^+ and L^- are not; everything else is Conway-Sloane 1992 plus the project's
own bridge ([rank-four-isospectrality-literature.md](../docs/rank-four-isospectrality-literature.md) section 1). The reviewer confirmed that fact
independently (section 2, line "isometries ker(L+) vs ker(L-): 2").

## 2. Independent replication

`code/rank4_level4_review_check.py` (filed at the workspace root; output in
`reports/rank4-level4/review-independent-check.txt`) uses only the previously
accepted modules `rank4_feasibility` (isometries, modular_data, determinant)
and `exact_theta` (complete_shells) plus local helpers for Gram matrices and
overlattice bases. It imports nothing from the new builder or verifier and
constructs the Conway-Sloane lattices from the reviewer's own transcription
(a workspace check script, not distributed). Results, all exact:

| Claim in [rank-four-level-four-zero.md](../docs/rank-four-level-four-zero.md) | Reviewer's result |
|---|---|
| det G = 1729, a dot b even | 1729, a dot b = 2 |
| modular data d0, N0, M, mu, B | 6916, 3458, 55328, 107520, 17920 |
| complete signed shells through B: 2238 shells, 2,382,392 vectors, all signed coefficients zero | identical (4,908,667 LDL nodes, about 15 s) |
| Aut(L) = {I,-I} | 2 elements, I and -I |
| G0 = [[12,-4,-4,-6],[-4,10,-1,6],[-4,-1,8,0],[-6,6,0,16]], det 6916 | identical |
| Aut(L0) = {I,-I} | 2 elements |
| O_v1 isometric to L^-(1,7,13,19), O_v2 to L^+(1,7,13,19) | 2 isometries each; 0 to the other member |
| O_v1, O_v2 not isometric; L^+, L^- not isometric | 0 isometries (complete LDL basis-image search) |
| ker b=(1,1,1,1) of L^+ and of L^- isometric, det 6916 | 2 isometries; both also isometric to G0 |
| Earnest-Nipp g1, g2 isometric to L^-(7,13,19,49), L^+(7,13,19,49) | 2 isometries each |
| det-3844 pair S1 (decomposable B + 2B), S2: non-isometric, isospectral to 2000 | 0 isometries; theta equal through norm 400 (reviewer's bound); S1, S2 are the b=(1,1,1,1) kernels of the det-961 pair (4 and 2 isometries; 0 crosswise) |

The cutoff proof rests on `rank4_feasibility.modular_data`, whose hypotheses
([rank-four-census-domain-and-cutoff.md](../docs/rank-four-census-domain-and-cutoff.md) section 2) are any positive definite integral G, b
nonzero, a dot b even; odd G is allowed, so the odd lattice here is in scope.
The builder calls that function and `cutoff_build` rather than reimplementing
them; the verifier recomputes the cutoff with `modular_check` from the
accepted feasibility verifier.

## 3. Replays of the delivered material

Run from the workspace root on 2026-10-20 with fresh output paths in the
reviewer's scratchpad (the delivered replays in `reports/rank4-level4/` were
not overwritten):

| Command | Verdict | Time |
|---|---|---|
| `python test_rank4_level4.py` | 12 tests OK | 20.6 s |
| `verify_rank4_level4.py --task-a task-a.json --bijection cs-bijection.json --certificates det1729-level4.json` | verified; 9 Task A rows (8 `klein_configuration_found`, 1 `no_index_two_neighbour_relation`), bijection verified, certificate `verified_level_4` | 10.7 s |
| `verify_rank4_level4.py --bijection cs-bijection.json --certificates other-level4-certificates.json` | verified; 7 certificates `verified_level_4` | 67.0 s |
| `verify_rank4_level4.py --task-d task-d-det24.json --census data/census-rank4/det24-bound32.json` | verified; 169 classes, 2235 level 1, 15,510 nonzero, 0 at levels 2, 3, 4 | 1.9 s |

Note for reruns: `python -I` cannot import the project modules (it drops the
script directory from the path); run without `-I` from the workspace root.

## 4. Reading of the code

- Task A (`analyse_pair`, `cross_table`): class histograms of L/2L through
  norm 2000 give the fifteen kernel thetas and fifteen overlattice thetas
  exactly; every full-prefix agreement goes to the complete basis-image
  search of `rank4_feasibility.isometries`, which raises on interruption and
  never returns a partial list. The verifier recomputes the histograms with
  its own box enumeration (exact Cauchy-Schwarz radii) and the isometries with
  `box_isometries` from the accepted feasibility verifier, and checks all 450
  cells per pair.
- Klein construction (`klein_from_sublattice_match`): L0, O1 = U(L1), O2 = L2,
  L and L' are built as doubled HNF bases in L2 coordinates. The verifier's
  `check_klein_construction` recomputes the indices (32, 16, 16, 16, 8), the
  inclusions, the distinctness of the three intermediate lattices, O2 = L2,
  O1 isometric to L1 by box search, the Gram of L, and that the characteristic
  cosets a/2 + L0 and a/2 + t + L0 are O1 - L0 and O2 - L0.
- Certificate (`certify_configuration`): the LLL transport of (a,b) is the
  correct one (a' = U^{-1} a, b' = U^t b); the coset lattices are rebuilt from
  the kernel basis; the links to L1, L2, the non-isometry (complete search
  plus the multiset of index-two sublattice theta prefixes), the full group,
  the stabilizer signs, the kernel group and its action on the two classes are
  all recomputed by the verifier from the Gram matrices, not copied.
- Bijection certificate: the piece the verifier does not test explicitly is
  that a single coordinate sign change maps M = 2D4 + Z(1,1,1,1) onto
  M' = 2D4 + Z(1,1,1,-1). That holds: an odd vector of M has coordinate sum
  0 mod 4 and the flip changes the sum by twice an odd number, giving 2 mod 4,
  the defining residue of M'; an even vector of 2D4 changes its sum by a
  multiple of four and stays in 2D4. Hence flip(r + 3M) = flip(r) + 3M', and
  with the nine images pairwise distinct modulo 3M' the pieces assemble to a
  norm-preserving bijection L^+ -> L^- for every diagonal form. The verifier
  checks the parity-and-sum structure of M and M', the containments
  3M < L^+ < M and 3M' < L^- < M', the index nine, the nine representatives,
  the flip rule and the distinctness. This is Conway-Sloane's Theorem 1 proof
  made explicit; the reviewer regards it as a valid all-level identity proof
  independent of the cutoff.
- Task D (`analyse_kernel`): a differing coset prefix is an exact proof of
  non-vanishing; a swap in Aut(L0) is an exact proof of vanishing at all
  levels; levels 2 and 3 are decided by exact isometries; level 4 would be a
  candidate only. The verifier recomputes every class from the census pack
  (hash checked) with box backends.

## 5. Points to state more precisely (none affects the verdict)

- "No level-3 configuration was found anywhere" has two universes. In Task A
  every Klein configuration comes from a known non-isometric pair, so by
  construction it is level 4; level 3 could only have appeared in Task D,
  whose universe is the 169 isometry classes of integral rank-four lattices
  with det(L0) <= 24 (exact-isometry dedup from the accepted census) and
  their 105 configurations each, 17,745 in all: 0 at level 3.
- Task B (symmetric L0 constructions) was skipped because Gate 1 was
  positive. That is what the brief asked; it means no statement is made about
  level-4 examples outside the Conway-Sloane family.
- At the time of this review the smallest known counterexample determinant
  was 1729 (odd lattice); the Eisenstein sweep of item 3 then lowered it to
  1156 ([rank-four-eisenstein-review.md](rank-four-eisenstein-review.md)). The census lower bound is det(L) > 24.
  No minimality claim is made in between.
  (Superseded 2026-10-21: both records are now 280, and 280 is the minimum by
  the determinant census of release 1.2.0; see the note in section 1.)
- The det-3844 pair (S1 decomposable, S2 indecomposable) is only
  isospectral-to-a-bound; it is not certified and should not be cited as a new
  isospectral pair without an all-level proof.
- The report lists Earnest-Nipp and the (7,13,19,49) member as separate rows;
  they are one configuration in two bases, as the report itself says.

## 6. Consequences recorded by the reviewer

- [rank-four-isospectrality-literature.md](../docs/rank-four-isospectrality-literature.md) sections 1, 3, 5 and 7 amended on 2026-10-20: the
  statement that the Conway-Sloane configuration "has index three and does
  not qualify directly" is withdrawn; the re-embedding and the determinant-1729
  example are described in section 3.
- The release folder, the manuscript and the deferred-research
  table are untouched; porting (the deferred row on a smaller counterexample
  is now partly answered by a different construction) is the user's decision
  when the release is next updated.
- Item 2 of the workspace to-do list (not distributed) is complete. Item 3 (smaller counterexample
  from the Eisenstein structure) was run next from a second
  handoff brief and is reviewed in
  [rank-four-eisenstein-review.md](rank-four-eisenstein-review.md) (result: determinant 1156).
- 2026-10-21: the external assessment [project-assessment-2026-10-21.md](project-assessment-2026-10-21.md) proved
  the all-level identity for the determinant-3844 pair S1, S2 of report
  section 2 (even Grams, level 62, weight-two Sturm bound 16, coefficients
  agree, no isometry) and refined the wording of the verdict ("no geometric
  explanation" means no single Euclidean isometry; the sign-change bijection
  is a piecewise norm-preserving bijection). Both are recorded as dated
  addenda in [rank-four-level-four-zero.md](../docs/rank-four-level-four-zero.md); the assessment's script was rerun by
  this session with identical output.
- 2026-10-21: ported to the release folder as version 1.1.0
  (`docs/rank-four-level-four-zero.md`, `reviews/rank-four-level-four-review.md`,
  code, tests, certificates and reports).

## 7. Files added by this review

- [rank-four-level-four-review.md](rank-four-level-four-review.md) (this note).
- `code/rank4_level4_review_check.py` and
  `reports/rank4-level4/review-independent-check.txt`.
- The three replay JSON files were regenerated in the reviewer's scratchpad
  only; their verdicts equal those of `reports/rank4-level4/replay-*.json`.
