# Coverage replay of Theorem S: all 502,740 face/marking cases (Grok Bot, 2026-10-23)

## Result
Every case has an exact disposition. **0 consistency failures, 0 uncertified cases, 0 mismatches with Claude's `cells.json`.**

| disposition | count |
|---|---|
| multiplicity_mismatch | 477,492 |
| infeasible (only mu = 0 is a nonnegative solution) | 0 |
| lower_face (nonzero nonnegative solutions exist, but every one has some mu_j = 0) | 19,992 |
| certified | 5,256 |
| uncertified | 0 |
| **total** | **502,740 = 2,394 x 210** |

## Universe and counts (all recomputed)
- Input is the 3 cones of `fan.json`, each with 10 rays, used as input data. That gives 3,069 nonempty ray subsets.
- **2,394** of those subsets are positive definite. This is checked by an exact LDL^T over Q of the sum of the rays, requiring every pivot to be > 0.
- **210** markings N with N^2 = 0 and rank N = 2. These were enumerated over all 2^16 matrices N.
- Minimal sets M_k for k = 1..15 at the sum-of-rays point:
  - The upper bound per class is the minimum of H(bits(k) + 2w) over w in {0, ±e_i, ±e_i±e_j}.
  - The enumeration is complete Fincke–Pohst on the exact rational LDL^T. Coordinate bounds are integers, computed via an exact `isqrt` of the rational radius, and every candidate is tested exactly. No floating point is used anywhere.
- Second interior point: random positive integer weights 1..9 per face. Its minimal sets are identical on all 2,394 faces, so there are 0 consistency failures.
- Further check: every v in M_k has the same vector (R_j(v))_j over j in S, so E does not depend on the choice of representatives. It held for all classes on all faces.

## Certificates (each verified exactly before it is written)
- **multiplicity_mismatch:** records the first unequal pair and both complete minimizing sets.
- **Feasibility:** an exact two-phase simplex (Fractions, Bland's rule) computes the maximal support J of {mu >= 0, E mu = 0}. It returns:
  - mu with E mu = 0, mu_j >= 1 on J, and mu_j = 0 off J;
  - when J != S, a Stiemke-type y with E^t y >= 0 and (E^t y)_j > 0 for every j not in J.
  Together these prove J is exactly the set of weights that can be positive.
  - J = S means the case is feasible.
  - J empty would mean infeasible; this never happens.
  - Otherwise the case is lower_face. These are certified by y; Claude's zero-coordinate rule was not used.
- **certified:**
  - mu > 0 exact, with E mu = 0.
  - An integral T with det ±1, mod-2 permutation equal to tau, and T^t H(mu) T = H(mu).
  - D (10 x |S|, columns sym(T^t R_j T − R_j)) satisfies **D = Z E** with Z rational, found by exact RREF and recorded. So T preserves H on the whole space {E mu = 0}.
- **T provenance:**
  - Main run (`replay.py`): it tried the saved T from `cells.json` first, and all 5,256 passed my own D = Z E check.
  - Independent second pass (`ownT.py`): ignores the saved T. For each of the 5,256 cases it runs my own complete isometry search at a random generic rational point of the cell, restricted to columns in the classes tau(e_i). **All 5,256 have such a T satisfying D = Z E.**
- **Comparison with `cells.json`, keyed by (cone, subset, tau):**
  - Claude certified-but-not-mine: 0.
  - Mine-but-not-Claude: 0.
  - Both sets have 5,256 cases, and they are identical.

## Dependencies trusted
- `fan.json` as input data: that its 3 cones and their GL4(Z) images cover the PD cone, with the stated ray sets. This is **UNVERIFIED here**; Grok is auditing the fan separately.
- That minimal sets are constant on relint(S), justified by the constant Delaunay subdivision on a face of the secondary fan. I only checked two points per face, so constancy is **UNVERIFIED** beyond those points.
- Python 3 standard library (`fractions`, `math.isqrt`). I did not import `rank4_six_cells.py`, `rank4_secondary_fan.py` or `rank4_six_minima.py`. I read `rank4_six_cells.py` only to confirm the conventions. No Magma, no determinant search.
- Scope: this replays the finite case analysis only. The implication from Theorem S to the minus-type conclusion was not audited here.

## Hashes (sha256)
- replay.py `6d62342f629082e114f7ef9a41f16138c8276f7ef67834df5c98cd5f0729902a`
- ownT.py `122c71f8fbe796e2fd7af0eca17ce3e701abbde1cf2cbefd7128008c0f997a82`
- fan.json `f4970a10fa7732e0c15f951ebf53071b8b5b7528bfed9df338eecc65e8e51617`
- cells.json `687b91803a90e7cd729ad20a35742a0666acb9ccfd374ecf21e179183fcdda5c`
- dispositions.jsonl.gz `d630fed22ff084c4e145bc0b8b6503391178319bb1d25104f8a8578199b4dd2c`

## Run
- Main replay ran on Grok Bot's box: single thread, nice 19, plain Python 3. Wall 396.8 s, CPU 396.7 s.
- The ownT pass took 37.4 s.
- Command: `python3 replay.py fan.json cells.json OUTDIR`, then `python3 ownT.py fan.json OUTDIR/dispositions.jsonl.gz`.
- Outputs: `summary.json` and `ownT.json`.
