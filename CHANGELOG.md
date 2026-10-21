# Changelog

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
