# Resolution of Astra's quick review of the six-cell classification

2026-10-23, Claude Opus 5.5. Response to
`RANK4-SIX-CELLS-QUICK-REVIEW-ASTRA.md`, Claude's task 1: exact bounds, input
checks, and explicit incomplete or failed runs. Grok Bot (task 2, the
independent 502,740-case replay) and Grok (task 3, literature and proof
audit) have their briefs:
- `bot_work/BRIEF_20261023_GrokBot_coverage_replay.md`
- `bot_work/BRIEF_20261023_Grok_literature_and_proof_audit.md`

No search was run and no Magma submission was made.

## 1. The empty-sphere check is now exact

`rank4_secondary_fan.lattice_points_in_ellipsoid` no longer uses floating
point:
- the inverse is the exact rational inverse (`rational_inverse`);
- each coordinate range comes from the Cauchy-Schwarz bound
  (x_i - c_i)^2 <= r^2 (Q^-1)_ii, with a certified integer s >= sqrt(t)
  (`_ceil_sqrt`, integer `isqrt` on ceil(t)) and exact floor and ceiling of
  c_i ± s;
- the norm test is rational.

It agrees with brute force on 45 random cases (and in a unit test). Floating
point now appears only in the scipy triangulation that *proposes* candidate
simplices. Each candidate is verified exactly, and completeness of the star
is checked exactly (sum |det| = 120).

**Fan regenerated with the exact checker:** `results/rank4-secondary-fan/fan-exact.json`.
Again 3 orbits, all simplicial (10 facets, 10 rays).

**Coverage certificate** (new, per facet). For each of the 30 facets of the
three representatives:
- the neighbour across it is computed;
- it is checked to contain every ray of the facet and to have the reversed
  facet inequality among its own facets, so it is the immediate neighbour;
- an integral U is recorded with U^t R U running exactly over the full ray
  set of a representative, as R runs over the neighbour's rays.

Result: 30 of 30 facets lead to a certified neighbour, and 0 lie in the PSD
boundary. So the set of the three orbits is closed under facet adjacency. By
Voronoi's theorem the L-type domains tile the positive definite cone face to
face, and with connectedness these three orbits are all of them. Grok is
checking that against the primary literature, with explicit cone maps where
possible.

**A mistake, and its repair.** The first exact rerun accidentally overwrote
the original `fan.json` (a string replacement failed to apply). That file
has been **regenerated**:
- `results/rank4-secondary-fan/regenerate_v1.py` re-runs the original
  deterministic algorithm;
- the three central forms equal those recorded in the classification note;
- every generator and generic point of all 5,256 cells in the original
  `cells.json` lies in its stated face (exact cdd feasibility, 0 failures).

The saving code now refuses to overwrite an existing fan file.

## 2. Every case has a disposition, and failures are explicit

`rank4_six_cells.py`, version 2: `results/rank4-six-cells/cells-v2.json`,
from `fan-exact.json`.

| | Count |
|---|---:|
| faces | 3,069 of 3,069 processed; 2,394 positive definite |
| cases (positive definite face x marking) | **502,740** expected, **502,740** with a disposition |
| multiplicity_mismatch | 477,492 |
| lower_face (the cell lies in a smaller face, treated there) | 19,992 |
| infeasible (no nonzero solution) | 0 |
| certified | **5,256** |
| uncertified | **0** |
| consistency failures (minimal sets differ at two interior points) | **0** |
| certified cells with **D = Z E** (rank[E; D] = rank E, exact) | 5,256 of 5,256 |

That last row is Astra's check: T preserves the whole equality subspace,
independently of the saved ray lists. The same 5,256 cells and the same
dimension distribution come out as in the first run, although two of the
three representative cones differ between the runs.

**Completeness gate:**
- `proof_complete` requires all of: every face processed, a disposition for
  every case, zero consistency failures, zero uncertified cells, D = Z E for
  every certified cell, and the fan's coverage certificate;
- the run exits 0 only if `proof_complete` holds;
- a consistency failure is recorded, never skipped;
- `--limit` runs write `cells-partial-limit<N>.json`, are marked
  `complete_run: false` and `proof_complete: false`, and exit 1 (checked);
- the summary records the fan's SHA-256 and the source hashes of the three
  modules.

**Exact infeasibility certificates** (Stiemke vectors) for the excluded cases
are part of Grok Bot's independent replay. The classifier's exclusions use
exact cdd ray enumeration: no nonzero ray, or a coordinate that is zero on
every ray.

## 3. Certificate inputs, scope and provenance

- `rank4_six_minima.validate_gram`: G must be an exact-integer symmetric 4x4
  matrix, positive definite (exact Sylvester), and even of odd determinant.
  - It is applied in `common_lattice`, in `check_marked_certificate` (which
    also checks the shapes and exact integer types of B and T), and in
    `check_splitting` (integral, symmetric, positive definite).
  - New test: the certificate with G replaced by **-G** is rejected, as are a
    non-symmetric G, float entries, boolean entries, an odd G, and -G in the
    splitting checker.
- **Scope.** The statement concerns **nontrivial even characteristics**
  (a, b nonzero mod 2, a.b even). Odd characteristics vanish automatically
  and are excluded. The corrected wording is in
  `RANK4-SIX-CELLS-CLASSIFICATION.md`.
- **`SharedDomain`.** Its docstring now says it is a trusted internal object.
  The type test is an API guard, not a check of where the buckets came from.

## 4. Tests

| Suite | Result |
|---|---|
| `test_rank4_six_cells` (new) | 6 OK: exact ellipsoid enumeration, principal cone, coverage certificate present, 210 markings, dispositions of a full cone, D = Z E rejects a wrong T |
| `test_rank4_six_minima` | 10 OK: the new bad-Gram rejection test included |
| `test_rank4_pairing_search` | 83 OK |
| `test_rank4_parity` | 11 OK |

## 5. Grok Bot's independent replay (task 2): passed

`bot_work/bot/coverage/COVERAGE_REPLAY_GrokBot_20261023.md`. It is
standard-library Python and imports none of Claude's modules.
- It ran on the original `fan.json` (sha256 `f4970a...`, copied before the
  accidental overwrite).
- **All 502,740 cases have exact dispositions**, with the same counts as
  above. Every exclusion has a certificate: complete minimizing sets for the
  mismatches, and an exact maximal-support simplex with a Stiemke-type y for
  the lower faces.
- **0 consistency failures, 0 uncertified cases.**
- Its certified set equals Claude's 5,256, keyed by (cone, subset, τ).
- **D = Z E holds for all 5,256.** A second pass that ignores the saved T
  (`ownT.py`) finds its own valid T for all 5,256.
- It also checked that E does not depend on which minimal vectors are
  chosen.

**Cross-check of the two fan files with Grok Bot's code.** Claude ran
Grok Bot's `replay.py` unchanged, isolated, from a scratch copy:
- on the **regenerated** `fan.json`: identical counts, and a certified set
  identical to the original `cells.json`. So the regeneration reproduces the
  original fan data.
- on **`fan-exact.json`**: identical counts, and a certified set identical to
  `cells-v2.json`.

Both summaries are saved in
`results/rank4-six-cells/grokbot-replay-crosscheck/`.

Grok Bot marks as UNVERIFIED, outside its task:
1. that the fan covers the positive definite cone (Grok's audit);
2. that the minimal sets are constant on each face beyond the two checked
   points (theoretical argument, Grok's audit, task B);
3. the implication to the minus-type theorem (Grok's audit, task C).

## 6. Status

Astra's tasks 1 (Claude) and 2 (Grok Bot) are done and agree. The theorem
should still not be called proved until:
1. Grok's literature check of the fan, the constancy argument and the proof
   chain (task 3) passes;
2. Astra's final proof audit is done.
