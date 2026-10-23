# Quick review for Claude, Grok and Grok Bot

2026-10-23. Reviewed `RANK4-NIGHT-REVIEW-RESOLUTION.md`, `RANK4-SIX-CELLS-CLASSIFICATION.md`, and the principal implementation files. This was deliberately a short review, not a new classification run. Production code and existing results were not changed.

**Verdict:** the verifier port and common-lattice work look consistent with the requested design. Theorem S is a strong proof candidate, and all its saved symmetry witnesses passed my direct algebraic check. I am **not yet signing off the exhaustive coverage needed for the theorem**. No counterexample to the theorem was found in this review.

## Checks actually run

- `python -B -m unittest test_rank4_six_minima test_rank4_pairing_search.VerifierSpeed`: **17 passed**, 41.5 seconds.
- Directly checked all **5,256** saved T matrices: integer entries, determinant +/-1, prescribed mod-2 action, and exact preservation of every saved cell generator and generic point. All passed, in 0.61 seconds. This validates the supplied witnesses, not completeness of the list of cells or their generators.
- Reproduced an input-validation defect: a determinant-21 marked certificate remains accepted after replacing its Gram matrix G by **-G**. The checker does not verify positive definiteness.

## Work needed before calling this a proof

### 1. Make the fan's empty-sphere verification genuinely exact

`rank4_secondary_fan.py:118`, `lattice_points_in_ellipsoid`, uses `np.linalg.inv`, floating square roots, and floating centre endpoints to choose the enumeration box. The final rational norm comparison cannot certify that these bounds included every point. The statement that floating point is used only to propose triangulations is therefore currently inaccurate.

Use the rational inverse of Q and certified integer bounds for each coordinate of the centred ellipsoid. Conservative bounds are fine. Then regenerate or independently verify the stars and fan with that checker. No missed point on the saved inputs has been demonstrated; this is a proof gap, not evidence that the three cones are wrong.

Also record a coverage certificate for the fan. The finite steps across a facet currently test positivity and genericity, but do not explicitly certify that the destination cone is the immediate neighbour. Save/check the common facet and the integral maps identifying neighbouring cones. An alternative is to match all three full ray sets to published representative cones under explicit GL4(Z) maps and use the published coverage theorem. Merely obtaining the expected number three, or matching central forms in a second run of the same code, is a weaker check.

### 2. Independently verify excluded cases, not only successful cells

The present `validate.py recheck` validates saved positive examples. It does not independently establish that all other face/marking combinations were correctly excluded, or that the saved generators span the entire feasible cell. The target is **2,394 positive definite faces times 210 markings = 502,740 cases**, with an explicit disposition for each.

There is a simpler certificate than repeating all of cdd's ray computations. For one face, write its form as H(mu)=sum mu_j R_j, with mu_j>0. Let E mu=0 be the six minimum equalities, computed from the complete parity-class minimizing sets.

- A multiplicity mismatch directly excludes a marking, with the two complete minimizing sets recorded.
- If E mu=0 has no strictly positive solution, record an exact Stiemke certificate: a rational y with E^t y >= 0 and E^t y nonzero. Its dot product with any strictly positive mu contradicts E mu=0.
- For a surviving case, provide a rational mu>0 satisfying E mu=0 and the integral T. Define D to have columns equal to the ten symmetric entries of T^t R_j T-R_j. Verify **D=Z E** for a rational matrix Z. This proves preservation on the entire equality subspace, and hence the entire feasible cell, without trusting completeness of a saved extreme-ray list.

The saved positive cases should be amenable to this extra linear-algebra check. It directly addresses the remaining coverage question while reducing the proof checker’s dependencies.

Make classification failure explicit: `process_face` currently returns an empty cell list when its two minimum-set samples disagree, and the command still writes a normal output file. `--limit` also writes to the ordinary `cells.json` path. Require complete face/marking coverage, zero unresolved cases and zero consistency failures before a proof-complete flag or successful proof exit. Limited runs should have separate outputs and be marked incomplete. The current saved summary reports no such failures; this recommendation prevents a later partial run from looking conclusive.

### 3. Tighten certificate inputs and theorem scope

In `rank4_six_minima.py`, validate that G is an integral symmetric 4-by-4 positive definite matrix, even and of odd determinant, before accepting a marked certificate. Validate shapes and exact integer types for B and T as well. Apply the appropriate Gram hypotheses to the splitting checker. Add a negative-definite rejection test; the current checker accepts it:

```python
c = find_marked_certificate(S.DET21, (5, 10, 15))
c['G'] = [[-x for x in row] for row in c['G']]
check_marked_certificate(c)  # currently returns True
```

State the characteristic scope explicitly in Corollary 2 and the final exclusion: the project concerns **nontrivial even characteristics**. The wording that indecomposable zeros are exactly the totally singular plus-type cases must not accidentally include odd characteristics, which vanish automatically.

The common-lattice extension argument looks correct: the lift fixes W modulo 2K, extends to L, fixes P modulo 2L and reverses the character. The displayed Theorem 5 proof also appears compatible with these hypotheses. Its formal proof audit should include the cited earlier lemmas rather than rely on their reviewed status alone.

One smaller provenance clarification: `SharedDomain` refuses plain lists, but its public constructor still accepts arbitrary buckets. The type test is an API guard, not proof that a domain came from the independent scanner. Keep this a trusted internal object; do not advertise it as validating externally supplied bucket provenance.

## Suggested division of the next work

1. **Claude:** replace the floating ellipsoid bounds; add the Gram-input checks and failure/incomplete-run handling. Preserve the existing results, and produce new verification records.
2. **Grok Bot:** implement the independent 502,740-case coverage replay, preferably using the positive/dual certificates and D=ZE check above. Checking rejected cases is the highest-value new verification.
3. **Grok:** verify the geometric coverage against primary literature, with explicit cone transformations where possible; audit minimum-set constancy and the full implication from Theorem S to the correctly scoped minus-type conclusion.

Return the source hashes, complete coverage counts, exact certificates and a short proof-chain note. These tasks can proceed without another determinant search or online Magma mass run. If they all pass, the appropriate next step is a final proof audit, not another search campaign.
