# Theorem S and the minus-type corollary: audited package (2026-10-23)

This is a frozen snapshot of everything the computer-assisted proof rests on,
**accepted in Astra's final audit** (`final_audit_astra/RANK4-THEOREM-S-FINAL-AUDIT-ASTRA.md`).
`MANIFEST.json` lists every file with its sha256 and original path. Nothing
here is edited after packaging; corrections are listed under Errata.

## The results

**Theorem S.** Let H be a positive definite real quadratic form on Z^4 and
τ = I + N an involution of F2^4 with N^2 = 0 and rank N = 2. Suppose that for
each of the six pairs {k, τk} of classes of Z^4/2Z^4 moved by τ, the two
classes have equal minimal norm and equal numbers of minimal vectors. Then
some T in GL4(Z) satisfies T^t H T = H and T ≡ τ (mod 2).

**Corollary (minus type).** An even, orthogonally indecomposable, positive
definite rank-four lattice with determinant ≡ 5 (mod 8) has no identically
vanishing theta series of a nontrivial even characteristic (a, b nonzero
mod 2, a.b even), **at any determinant**. More generally, in an even positive
definite rank-four lattice of odd determinant, every such zero is a symmetry
zero.

Scope:
- Odd characteristics are excluded.
- Decomposable minus-type lattices can have zeros (the det-21 control).
- This does not classify all isospectral quaternary pairs, and says nothing
  about odd lattices or even determinants.

## The proof chain

The full argument is in `documents/RANK4-SIX-CELLS-CLASSIFICATION.md`
section 1 and in the final audit.

1. **Plane invariance** (`documents/RANK4-STAGE2-STRATEGY.md` section 2; the
   multiplier argument as corrected in the final audit). A zero of one
   characteristic gives zeros of all three shifts of its plane P.
2. **The common lattice** (`documents/RANK4-SIX-MINIMA-COMMON-LATTICE.md`).
   Each of the six halves is θ(v + 2K) - θ(w + 2K), so a zero makes all six
   minimal shells balanced. That is Theorem S's hypothesis for
   H = B^t G B / 2 with the marking τ.
3. **Theorem S** (computed; section "Computation" below).
4. **Extension.** T ≡ τ preserves W = 2L/2K, extends to h in Aut(L), fixes P
   mod 2L, and reverses the character. So h is a symmetry witness.
5. **Theorem 5** (`documents/RANK4-MINUS-TYPE.md`): an indecomposable,
   minus-type lattice has no symmetry witness.

## Computation

- **The fan.** Voronoi's secondary fan of quaternary forms has three
  GL4(Z)-orbits of maximal cones, all simplicial with 10 rays.
  `data/fan-exact.json` was computed with exact arithmetic and a per-facet
  coverage record. Both saved fans map bijectively onto Δ, Δ′, Δ″ of
  Vallentin, Weißbach and Zimmermann (arXiv:2407.03513 section 2.1;
  Indag. Math. 36 (2025), p. 992), by explicit unimodular maps checked on all
  ten rays (`checks_claude/match-published.json`,
  `final_audit_astra/published-cones-and-hashes.json`). Coverage of the
  positive definite cone by these cones is classical: Voronoi; Delone 1929
  p. 161; arXiv:math/0601084 Theorem 2.3.
- **The classification** (`sources/rank4_six_cells.py`). For every positive
  definite face (2,394) and every marking (210), there is an explicit
  disposition:

  | Disposition | Count |
  |---|---:|
  | multiplicity mismatch | 477,492 |
  | lower face (Stiemke-type certificate) | 19,992 |
  | certified (μ > 0 with E μ = 0, an integral T ≡ τ, and **D = Z E**) | 5,256 |
  | uncertified | **0** |
  | total | 502,740 |

  D = Z E means T preserves H(μ) for every real μ with E μ = 0.

**Three independent confirmations:**
1. Grok Bot's standard-library replay (`independent_replay_grokbot/`):
   identical dispositions, certificates for every exclusion, and his own T
   search for all 5,256.
2. Grok's literature and proof audit (`literature_audit_grok/`).
3. Astra's direct bundle check (`final_audit_astra/check_bundle.py`,
   `bundle-audit.json`): all 502,740 cases, recomputing E, the dispositions
   and D = Z E, without using Claude's classifier or cdd.

## Verify

The package's own code, run from a copy of `sources/` (checked
2026-10-23, 8-10 s on 12 workers):

```text
python -B rank4_six_cells.py run --workers 12 --fan <package>/data/fan-exact.json --out check.json
  -> 502,740 cases; 477,492 / 19,992 / 5,256 certified; 0 uncertified; proof_complete true;
     source hashes equal MANIFEST "audited_hashes"
```

The independent replay:

```text
python3 independent_replay_grokbot/replay.py data/fan.json data/cells.json OUTDIR
```

**`proof_complete: true` alone is not the proof artifact** (see Erratum 4).
The proof rests on the audited fan hashes, the published-ray identification
and the certificate bundle checks above.

## Provenance and the fan files

| File | Origin |
|---|---|
| `data/fan.json` | The first fan run. **Overwritten by mistake** on 2026-10-23 at 00:51 by the exact rerun, then **regenerated** by the original deterministic algorithm (`checks_claude/regenerate_v1.py`, `regenerate_v1.log`). It is not the untouched original file. Its consistency is established by: the three central forms equal the recorded ones; every generator of the 5,256 cells of `data/cells.json` lies in its stated face; Grok Bot's replay (made on a pre-overwrite copy, sha256 `f4970a...`) gives identical results on it; its rays map onto the published cones. |
| `data/fan-exact.json` | The rerun with exact empty-sphere bounds and per-facet coverage records. |
| `data/cells.json` | The original classification, on the original fan. |
| `data/cells-v2.json` | The rerun on `fan-exact.json`, before two small edits to the code. |
| `data/cells-v3.json` | The rerun on the audited sources. Its recorded hashes equal `MANIFEST.audited_hashes`. |
| `data/cells-astra-final-20261023.json` | Astra's fresh rerun in the final audit. |

## Errata from the final audit

These do not affect the conclusion.

1. **Grok's audit C.4** (`literature_audit_grok/PROOF_AUDIT_Grok_20261023.md`):
   an integral matrix of order 8 has minimal polynomial dividing x^8 - 1, and
   **Φ8 = x^4 + 1** is required. It is not the fourth cyclotomic polynomial.
   The original Theorem 5 note uses the correct argument.
2. **Plane invariance.** The automorphy multiplier depends on the
   characteristic. Use the final audit's argument: at a fixed positive
   definite period, the theta function in its elliptic variable is not
   identically zero, so the invertible change of that variable rules out a
   zero multiplier. Do not rely on the zero-characteristic value.
3. **References.** The explicit cones are in **section 2.1 of
   arXiv:2407.03513**. The numbered coverage theorem ("Theorem 2.3") is in
   **arXiv:math/0601084**. Some notes attributed Theorem 2.3 to
   arXiv:2407.03513; keep the two references distinct.
4. **`coverage_ok`** in the packaged `rank4_six_cells.py` checks only the
   number of facet records. It is not a verifier of their content, or of an
   arbitrary fan. A distribution must pin the audited fan hashes or run a
   content-level coverage check. The working copy in the repository now does
   both (after this package; not part of the audited sources).
5. **Fan provenance:** keep the original, exact and regenerated fan
   provenance, and the overwrite explanation, as recorded above.

## Related but separate

- **Theorem R / R′ (marked vonorm rigidity)** is a separate computed
  candidate, not needed for Theorem S. It awaits Grok Bot's independent
  replay. See `RANK4-MARKED-RIGIDITY.md` in the repository.
- No further minus-type determinant census is needed for this pairing
  question.
