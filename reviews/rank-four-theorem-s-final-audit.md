# Final audit: Theorem S and the minus-type corollary

*Curated copy for release 1.4.0 of `RANK4-THEOREM-S-FINAL-AUDIT-ASTRA.md`, the final audit by a Codex session (Astra). Only the five workspace links are relocated to the release layout; the original file is preserved in [`reports/rank4-theorem-s/package/final_audit_astra/`](../reports/rank4-theorem-s/package/final_audit_astra/). The release copies of its two checkers are `code/rank4_theorem_s_check_bundle.py` and `code/rank4_theorem_s_check_published.py`. This is an AI-assisted audit with independent computation, not peer review; the theorem and its proof are stated in [rank-four-theorem-s.md](../docs/rank-four-theorem-s.md).*

2026-10-23. This audit resolves the outstanding questions in `Astra-handoff-20261023-morning.md` and supersedes the candidate-only verdict of my quick review.

**Verdict: I accept Theorem S and its correctly scoped corollary as a computer-assisted proof.** The finite computation is exhaustive, its positive and negative certificates have been checked, and the geometric and lattice-theoretic implications close. The result applies at every determinant, not just through 1997.

Precisely: an even, orthogonally indecomposable, positive definite rank-four lattice with determinant congruent to 5 modulo 8 has no identically vanishing theta series of a nontrivial even characteristic. More generally, in an even positive definite rank-four lattice of odd determinant, every such zero is a symmetry zero. Odd characteristics are excluded. Decomposable minus-type lattices can have these zeros, as the determinant-21 control shows.

This closes the particular minus-type search. It does not classify all isospectral quaternary pairs or assert absence of level-4 zeros for odd lattices or even determinants.

## Checks performed for this audit

| Check | Result |
|---|---|
| Current classifier, fresh output, 12 workers | 502,740 dispositions; 5,256 certified; 0 unresolved or consistency failures; 8.2 seconds |
| Tests of the classifier and common-lattice certificates | 16 passed, including the negative-Gram rejection; 45.4 seconds |
| Fresh exact parity-class minima | All 2,394 positive definite faces, using Grok Bot's separate standard-library enumeration |
| Direct certificate-bundle audit | All 502,740 cases checked, including every exclusion; 22.1 seconds total, with 16 workers for the minima stage |
| Uniform positive certificates | Every integral T has determinant +/-1, the correct mod-2 action, and explicitly satisfies D=ZE |
| Coverage accounting | Every positive definite face occurs with each of the 210 distinct markings exactly once |
| Published cone identification | Both saved fan files map bijectively to the three published cones, checking all ten rays under each integral unimodular map |
| Provenance | The three current source hashes and `fan-exact.json` match the hashes recorded in `cells-v3.json` and the handoff |

The independent bundle has exactly 477,492 multiplicity mismatches, 19,992 lower-face exclusions and 5,256 certified cases. The audit does not trust the stored disposition labels, E, or the asserted D=ZE: it recomputes these quantities and checks the equations. The certified set matches the original `cells.json` exactly. This also ties the regenerated original fan to the actual independent certificates, resolving the overwrite incident for the purpose of the proof.

Artifacts:

- [Fresh classifier output](../reports/rank4-theorem-s/package/data/cells-astra-final-20261023.json.gz)
- [Direct bundle checker](../reports/rank4-theorem-s/package/final_audit_astra/check_bundle.py)
- [Bundle audit result and input hashes](../reports/rank4-theorem-s/package/final_audit_astra/bundle-audit.json)
- [Published-cone and source-hash checks](../reports/rank4-theorem-s/package/final_audit_astra/published-cones-and-hashes.json)

The direct checker uses Grok Bot's exact short-vector enumerator, whose bounds and recursion I reviewed. It does not call Claude's classifier, cdd, Grok Bot's simplex solver, or his `solve_left` routine. Thus it checks the certificate equations without repeating the algorithms that found the feasibility and row-space certificates. It does not need to search for a new T: a valid supplied T proves existence.

## Why the computation proves Theorem S

The three simplicial secondary cones and their explicit rays are given in section 2.1, page 5 of [Vallentin, Weissbach and Zimmermann](https://arxiv.org/pdf/2407.03513). Their faces cover positive definite forms up to integral change of basis. I read the displayed matrices and independently transcribed them for the cone-match check; R11 is replaced by its primitive half, which does not change its ray. Recomputing the 52-orbit quotient is unnecessary because the computation treats all faces with repetitions.

Each ray is positive semidefinite. The kernel of a positive linear combination of such rays is their common kernel, so a face has positive definite relative interior exactly when the sum of its rays is positive definite. This justifies the 2,394-face selection, including boundary faces of the maximal cones.

For any lattice vector v, minimality in v+2Z^4 is equivalent to v/2 lying in the Voronoi cell of zero: the two defining norm differences differ by a factor of four. By central symmetry and convexity, this is equivalent to the Voronoi cells of zero and v meeting, hence to their sharing a Delaunay cell. The Delaunay subdivision is constant on the relative interior of a secondary-cone face. Therefore complete parity-class minimizing sets, including their multiplicities, are constant there. The two sampled interior points are a check, not the justification for constancy.

All markings are covered: there are 35 choices of W=ker(N)=im(N), and six isomorphisms V/W to W, giving 210 rank-two square-zero N. The enumeration independently finds these 210. Changing an integral basis conjugates the marking, so treating every marking on the three representative cones covers every marked form.

On a fixed face, write H(mu)=sum mu_j R_j with all mu_j positive. After the multiplicity test, the required minimum equalities are E mu=0. Each excluded case carries an exact y with E^t y nonnegative and nonzero. If mu were strictly positive, its scalar product with E^t y would be positive, contradicting E mu=0. Such points may lie on a proper face, which is treated separately.

Every surviving case has a strictly positive rational solution mu, an integral unimodular T with the requested reduction, and a checked identity D=ZE, where column j of D represents T^t R_j T-R_j. Thus E mu=0 implies T^t H(mu)T=H(mu) for **every real** solution, including irrational weights. The proof consequently covers all positive definite real forms, not just rational samples. The certificates cover all cases without omissions or duplicates.

## From Theorem S to the lattice statement

The common-lattice construction and extension argument in [RANK4-SIX-MINIMA-COMMON-LATTICE.md](../reports/rank4-theorem-s/package/documents/RANK4-SIX-MINIMA-COMMON-LATTICE.md) are valid. For P a B-isotropic plane, K is its inverse image in L. Each of the six signed halves is a difference of two ordinary 2K-coset theta series. Their pair-swapping marking has rank(tau-I)=2 and square zero. Theorem S supplies an integral lift, which preserves W=2L/2K and hence extends to an automorphism h of L.

The extension fixes P modulo 2L. For distinct nonzero a,c in P, the class j(a) annihilates a but not c under the perfect W-by-P pairing. Consequently B((ha-a)/2,c)=1 modulo 2. This is exactly the sign-reversing symmetry witness required in Theorem 5. The same h works for every ordered basis of P.

Plane invariance also holds. For an even rank-four Gram G of odd determinant d, expansion modulo four gives

    d = (G12 G34 - G13 G24 + G14 G23)^2 modulo 4,

so d=1 modulo 4. Put A=adj(G) and delta=(3d+1)/2. The integral symplectic matrix with blocks (2I,3G; A,delta I) sends zG to ((2z+3)/(dz+delta))G. The characteristic formula sends (a,b) to (c,Ga) modulo two, where c=G^{-1}b. Together with z->z+1 this generates the ordered bases of P. This application follows directly from [DLMF 21.5.9](https://dlmf.nist.gov/21.5.E9); the scalar fractional-linear transformation preserves the upper half-plane because its determinant is one.

For completeness, the multiplier is nonzero for **each** characteristic: at any fixed positive definite imaginary period, the theta function as a function of its elliptic variable is not identically zero, by its Fourier expansion. The invertible change of that variable in the transformation formula therefore rules out a zero multiplier. At elliptic variable zero the determinant factor is also nonzero, so identically vanishing one-variable theta constants transfer. This avoids relying on the zero-characteristic value to justify other characteristics' multipliers.

The a=c case cannot vanish: evenness makes the sign constant on each norm shell. Otherwise a,c span P, so plane invariance and the six-half decomposition apply.

I checked the preceding symmetry results as well. An involution witness would give two rank-two eigenlattices. Their 2-primary glue has exponent two and order 2^k, with k=1 or 2. The first case is incompatible with the determinant modulo four of an even binary Gram; the second forces all inner products in an eigenlattice to be even, contradicting the witness sign. Thus the 2-part of a witness has order four or eight. In rank four an order-eight integral matrix has characteristic polynomial Phi_8(x)=x^4+1. An order-four element of Sp_4(F2) fixes only one nonzero vector, so it cannot stabilize the independent a,c. The remaining order-four case either supplies h^2=-I or an integral splitting involution, excluded by indecomposability. Finally h^2=-I gives the totally singular plane (I+h)L/2L, forcing plus type. This proves the stated minus-type exclusion.

## Corrections and packaging notes

These do not reopen the mathematical conclusion above:

1. Grok's audit C.4 has a transcription error: an order-eight matrix has minimal polynomial dividing **x^8-1**, with **Phi_8=x^4+1** required. It is not the fourth cyclotomic polynomial and does not divide x^4-1. The original Theorem 5 note uses the correct argument.
2. Replace the zero-characteristic multiplier explanation with the argument above, or explicitly apply the transformation and its inverse. The multiplier depends on the characteristic.
3. The handoff attributes a general fan theorem to “Theorem 2.3” of arXiv:2407.03513. The explicit cones there are in **section 2.1**; Grok's numbered Theorem 2.3 citation is to **arXiv:math/0601084**. Keep the references distinct.
4. `coverage_ok` currently checks only the number of facet records. It is not a standalone verifier of their content or of an arbitrary user-supplied fan. The present theorem is supported by the independently checked **published-ray identification** above. A distribution should either pin those audited fan hashes or run the actual ray/coverage certificate checker before accepting another fan. A bare `proof_complete: true` is not the proof artifact.
5. Retain the original, exact and regenerated fan provenance and the accidental-overwrite explanation. The independent bundle and exact ray matches now establish the needed consistency; silently describing the regenerated file as the untouched original would be inaccurate.

## What to do next

**Theorem R/R-prime:** keep it as a separate computed candidate until its full independent replay is complete. Theorem S does not require it. A lack of a located published proof is not yet a novelty claim, and symmetry lifting within one integral isometry orbit should not be confused with a global marked-reconstruction theorem for arbitrary pairs of forms.

**For finding a new pairing:** the present even, indecomposable, odd-determinant minus-type route is closed. The useful discovery work is now in odd lattices or even determinant, using the known determinant-280 odd and determinant-4480 even level-4 examples as controls. The mod-2 setup needs redesign: for even determinant, G modulo two is singular and c=G^{-1}b may not exist. The six-minimum implication proved here must not be carried over without its hypotheses. Start from the sign-kernel and its two relevant cosets, then retain exact theta/isometry checks and the known-family exclusions.

My recommended split is to have one effort package Theorem S and its certificates, a second independently replay R-prime, and the discovery effort build the odd/even-determinant variant. No further minus-type determinant census is needed for this pairing question.
