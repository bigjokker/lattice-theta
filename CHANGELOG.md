# Changelog

## 1.4.0 (2026-10-23)

- The minus-type parity conjecture is proved at every determinant
  (Theorem S; `docs/rank-four-theorem-s.md`, section 7 of the follow-up
  note). An even, orthogonally indecomposable, positive definite rank-four
  lattice with det = 5 mod 8 has no identically vanishing nontrivial even
  characteristic. More generally, in an even rank-four lattice of odd
  determinant every such zero is a symmetry zero; with Theorem E, in an
  indecomposable lattice it comes from a square root of -I on a totally
  singular plane (plus type only). This answers the plus-type companion
  question of the deferred-research table for nontrivial even
  characteristics.
- Theorem S (computer-assisted): if a marking tau = I + N (N^2 = 0,
  rank 2) of Z^4/2Z^4 pairs classes with equal minimal norms and equal
  minimal multiplicities, then tau lifts to an automorphism of the form.
  It is proved by an exact classification over the faces of Voronoi's three
  primitive cones (matched to the published Delta, Delta', Delta''): 2,394
  positive definite faces x 210 markings = 502,740 cases, every one with an
  exact certificate (477,492 multiplicity mismatches, 19,992 Stiemke
  exclusions, 5,256 symmetry certificates verified on the whole equality
  subspace).
- Final audit by a Codex session (`reviews/rank-four-theorem-s-final-audit.md`),
  independent standard-library replay by Grok Bot and literature and proof
  audit by Grok; the audited files are preserved byte for byte in
  `reports/rank4-theorem-s/package/` (manifest with sha256; files over one
  million bytes stored gzipped, listed in its `COMPRESSED.json`). Release
  copies of the audit's checkers, `code/rank4_theorem_s_check_bundle.py` and
  `code/rank4_theorem_s_check_published.py` (standard library only), run in
  CI. The generating classifier in the package needs numpy, scipy and
  pycddlib; nothing else in the repository does.
- Status notes added to the minus-type report, the stage-two reductions and
  the deferred-research table; README, citation metadata, indexes,
  reproduction guide, export provenance and export validation updated;
  checksums regenerated.

## 1.3.0 (2026-10-21)

- Even rank-four lattices of odd determinant (section 6 of the follow-up
  note). Theorem E: in an even, orthogonally indecomposable rank-four lattice
  of odd determinant, every symmetry zero comes from an automorphism h with
  h^2 = -I on the plane (h+I)L/2L, so minus-type lattices (det = 5 mod 8)
  have no symmetry zeros at any determinant. Theorem F: vanishing depends
  only on the plane span(a, c); at odd nonsquare determinant no zero has
  level 2; in indecomposable minus type the two coset lattices are never
  isometric, so any zero would have level 4. The minus-type parity
  conjecture is thereby reduced to the absence of isospectral non-isometric
  coset pairs in that setting; it remains open.
- Proposition G: complete isometry-class lists of all even rank-four
  lattices with odd det <= 501 (1,890 classes, 1,696 indecomposable; new
  lists for det = 5 mod 8, 202..501, and det = 1 mod 8, 281..497, both
  independently replayed). Every zero of an indecomposable class lies on the
  plane of a square root of -I; the only such classes are det 121 and a new
  Gaussian class of det 361 = 19^2.
- Reports `docs/rank-four-minus-type.md` and
  `docs/rank-four-stage-two-reductions.md` (the latter by a Codex session);
  review `reviews/rank-four-stage-two-review.md` (line-by-line check of
  Lemma A, Theorem E and Theorem F, an independent computation over 31,500
  configurations). Scripts `code/rank4_minus_type.py`,
  `code/rank4_minus_type_constructions.py`, `code/rank4_stage2_checks.py`,
  `code/rank4_stage2_review_check.py`, `code/rank4_plus_type_census.py`;
  preserved outputs under `reports/rank4-minus-type/`, `reports/rank4-stage2/`
  and `reports/rank4-plus-type/`; seven tests (default suite 82, complete
  suite 166); CI steps for the stage-two checks and the plus-type replay and
  plane comparison. README, citation metadata, indexes, deferred-research
  table and reproduction guide updated.

## 1.2.0 (2026-10-21)

- The minimum determinant of a rank-four counterexample to the symmetry
  converse is 280 (Theorem D of the follow-up note). Every positive definite
  integral rank-four lattice with det <= 279 satisfies the converse: complete
  isometry-class lists of the four regions not covered by the determinant-24
  census and the 1.1.0 class lists (odd lattices with 4 !| det, 169..279;
  odd lattices with 4 | det, 25..279; even lattices of odd determinant,
  202..279; even lattices with 4 | det, 65..279) give 14,987 classes, 6,239
  of them indecomposable, with all 655,095 even characteristics decided and
  every one of the 5,330 zeros witnessed by an automorphism; every chunk is
  independently replayed. At determinant 280 exactly one class fails the
  converse (243 odd and 27 even classes censused). Script
  `code/rank4_gap_census.py`; preserved per-determinant chunks, region
  summaries, replays and `closure-25-279.json` under
  `reports/rank4-gap-census/`; report `docs/rank-four-minimum-determinant.md`;
  four tests (default suite 75, complete suite 159); a CI step replaying the
  determinant-280 chunks and checking the closure summary.
- The follow-up note, README and citation metadata now state the
  determinant-280 level-four zero (Theorem A'') and the minimum (Theorem D)
  as results. The dated correction notes of 1.1.1 in the archival reports and
  reviews are completed with the minimality statement; the deferred-research
  table marks the global-minimum question answered.

## 1.1.1 (2026-10-21)

- Correction: the determinant-4480 level-four lattice has every Gram entry
  even, and halving it gives an odd lattice of determinant 280 with
  automorphism group {I, -I} and a level-four zero. It is the smallest known
  level-four zero and the smallest known symmetry counterexample; the
  bracket for the minimum counterexample determinant is now
  24 < det <= 280 (was 24 < det <= 1156). Certificate
  `data/rank4-level4/det280-level4.json`; independent replay
  `code/verify_rank4_det280.py` (checks the determinant-4480 certificate
  with the unchanged level-four verifier, then re-proves the group, the zero
  by the modular cutoff, and non-isometry of the coset lattices directly);
  four tests (default suite 71, complete suite 155); CI step. Dated
  correction notes added to the level-four, Eisenstein and literature
  reports, their reviews, the assessment, the follow-up paper and the
  deferred-research table.

## 1.1.0 (2026-10-21)

Rank-four follow-up results, stated in
[papers/rank-four-follow-up-results.md](papers/rank-four-follow-up-results.md).

- Level four exists in rank four: an odd lattice of determinant 1729 has an
  identically vanishing even characteristic whose two coset lattices are the
  Conway-Sloane pair L^-(1,7,13,19), L^+(1,7,13,19), isospectral and not
  isometric. Eight such configurations are certified
  (`data/rank4-level4/`, report `docs/rank-four-level-four-zero.md`).
- A smaller counterexample to the rank-four symmetry converse: an even
  lattice of determinant 1156 with automorphism group {I, -I}, level 2, the
  same Eisenstein construction as the determinant-23,716 example. The
  bracket for the minimum counterexample determinant is 24 < det <= 1156
  (`data/rank4-eisenstein/det1156-level2.json`, report
  `docs/rank-four-eisenstein-family.md`).
- The rank-four parity property: both census conjectures are false, with
  complete class lists. Odd lattices with determinant not divisible by four
  fail first at determinant 49; even lattices of odd determinant fail first
  at determinant 121 (a Gaussian-integer structure with witnesses squaring
  to -I). Four structural theorems, including the absence of involution
  witnesses for even indecomposable lattices of odd determinant
  (`data/rank4-parity/`, report `docs/rank-four-parity-property.md`).
- Two deductions from the independent assessment of 2026-10-21: the ordinary
  theta series does not determine the parity property in rank four (the
  determinant-961 Conway-Sloane pair), and the determinant-3844 pair found in
  the level-four search is isospectral at every level (weight-two Sturm
  bounds; `code/check_small_pairs.py`).
- A literature note on the isospectrality sources (Schiemann, Conway-Sloane,
  Earnest-Nipp, Cervino-Hein) with the Klein-configuration bridge
  (`docs/rank-four-isospectrality-literature.md`).
- New code (`code/rank4_*_search.py`, `code/verify_rank4_*.py`, independent
  review checks), 36 new tests (default suite 67, complete suite 151),
  preserved outputs and replays under `reports/`, four review notes and the
  assessment under `reviews/`.
- Large preserved reports are gzip files listed in
  `reports/compressed-reports.json`; `python scripts/prepare_data.py --reports`
  expands them. The deferred-research table is updated; the two parity
  conjectures are recorded as false and the assessment's new directions are
  listed.

The determinant-23,716 manuscript is unchanged apart from a status note.

## 1.0.0 (2026-09-27)

Initial release: the symmetry converse in ranks one through three, the
rank-four counterexample of determinant 23,716, the low-rank parity
classification, the root-lattice classification and the determinant-24
censuses. Zenodo version DOI 10.5281/zenodo.23002465.
