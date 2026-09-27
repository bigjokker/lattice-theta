# Astra review: exact rank-four counterexample

Reviewed 2026-09-27 against `Astra-rank4-counterexample-review.md`.

**Verdict: ACCEPT.** The proposition in `papers/rank-four-counterexample.md` is proved.
The displayed primitive positive definite integral rank-four lattice, of
determinant 23,716, has an identically zero even half-characteristic and no
epsilon-one element in its full characteristic stabilizer. The unrestricted
rank-four symmetry converse is therefore false. This achieves this phase's
scientific rank-four goal. No mathematical correction is required.

The decisive checks are as follows.

1. **Lattice and characteristic.** Independently computed leading principal
   minors are 12, 48, 1232, 23716. They prove positive definiteness and give
   the stated LDL pivots 12, 4, 77/3, 77/4. The entries have gcd one (in
   particular, an entry is -1), so the integral Gram matrix is primitive.
   The given binary vectors satisfy a dot b = 0 exactly and hence define an
   even characteristic with phase one.

2. **Complete cosets and all-level vanishing.** Exact multiplication verifies
   V^t M V = 2G, det(V) = -8, R^t M R = M, det(R) = 1, and R^3 = I.
   The index-two sign kernel has basis K = diag(2,1,1,1), and
   2VK = 4W with det(W) = -1. Consequently the images of the two sign sets
   are the entire cosets 2e2 + 4Z^4 and 2e1 + 4Z^4, not merely subsets of
   them. The offsets follow from Va = 2e2 and
   V(a + 2e1) = -2e1 + 4e2. Since Re1 = e2 and R(4Z^4) = 4Z^4, R maps
   the entire negative coset bijectively onto the entire positive coset,
   preserving w^t M w / 2. Their unsigned coefficients agree at every
   norm. Absolute convergence, supplied by positive definiteness, justifies
   subtracting their theta series and proves F identically zero. No finite
   coefficient bound or modular-form cutoff enters this implication.

3. **Ambient equality versus an original lattice symmetry.** The next
   image under R is 2e1 + 2e2 + 4Z^4, distinct from both sign cosets.
   For an additional direct check, in the original basis the rotation is

       V^(-1) R V = [[1,    2,0, 0],
                    [-3/2,-2,0, 0],
                    [0,    0,0,-1],
                    [0,    0,1,-1]].

   Its nonintegral entry shows directly that this ambient rotation is not
   an automorphism of the original integral lattice. The complete group
   calculation below also excludes every other possible symmetry witness.

4. **Full automorphism group.** Independently inverting G gives diagonal
   (4/11, 3/11, 4/77, 4/77). The coordinate bound
   v_i^2 <= Q(v)(G^(-1))_ii follows from Cauchy-Schwarz in the G inner
   product, and confines every possible automorphism column to the stated
   315-point box. Direct enumeration of all those points gives exactly
   2, 6, and 12 vectors at norms 12, 16, and 28, respectively, matching
   every representative and negative in the proof. Exhaustively imposing
   the Gram products on all candidate column tuples gives precisely I
   and -I. The paper's shorter elimination is also correct: after fixing
   v1 = e1, the allowed columns are the two listed choices for v2 and v3
   and only e4 for v4; their products with e4 exclude the alternative
   columns by -4 != 2 and -11 != -14. Thus the group computation is
   complete, not a subgroup certificate. The stated orthogonal
   indecomposability corollary follows as well.

5. **Stabilizer, epsilon, and original normalization.** Both group elements
   fix a modulo 2 and b modulo 2 under inverse transpose. For I the dual
   shift is zero; for -I it is -b. Hence epsilon is respectively zero and
   -a dot b = 0. To check the original definition, take the primal shift
   xi = a/2 and the dual shift with basis coordinates delta = G^(-1)b/2.
   Then 2xi is integral, 2delta is dual-integral, and the phase is
   exp(pi*i*x dot b + pi*i*(a dot b)/2) = (-1)^x0.
   Also |x + xi|^2 = Q(2x+a)/4, so the original theta function is exactly
   F(tau/8). Thus the claimed identity and stabilizer sign use the
   project's original normalization.

Verification performed:

- `python scripts/run_tests.py --all`: all five tests passed,
  including acceptance with zero corroboration bound, rejection of altered
  data, and failure on interrupted enumeration.
- Fresh replay saved to `reports/rank4-counterexample/astra-box.json` using
  the supplied certificate. Its parsed contents agree exactly with
  `independent-box.json`, including input SHA-256
  `a1effafe9871bbb833a39a187d98f718879986d24fd23444e167b8635b55e36c`.
- A separate standard-library calculation, importing no project modules,
  checked the principal minors, inverse diagonal, all matrix and offset
  identities above, the complete 315-point enumeration, all compatible
  column tuples, and both stabilizer signs using exact integer/rational
  arithmetic. All checks passed.

The only presentation suggestion is optional: include V^(-1)RV above when
explaining why the ambient map supplies no original lattice symmetry.
The existing proof is already sufficient without this addition.

Proceed with the bounded closing write-up: integrate the proposition and
proof into the rank-four conclusion, record this accepted review and its
replay in the review snapshot, and distinguish the unrestricted negative
result from the previously accepted rank-one through rank-three results
and determinant-at-most-24 census. Those earlier results were taken as
context and were not re-reviewed here. Existing reports, manuscripts,
inputs, inventories, and ZIPs were preserved. No new unrelated research
lead arose; no side-quest entry is needed. No larger census, rank-five
continuation, minimum-determinant search, or priority claim is part of this
acceptance.
