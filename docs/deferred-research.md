# Deferred research alternatives

Prepared 2026-09-27; updated 2026-10-21 for releases 1.1.0 and 1.2.0. These are possible
directions, not active assignments. The original scientific goal (the
symmetry converse) is complete; the rank-four follow-up of October 2026 is
stated in [the follow-up note](../papers/rank-four-follow-up-results.md).

## Questions recorded in release 1.0.0, with their 1.1.0 status

| Question | Status after 1.1.0 |
|---|---|
| What replaces the indecomposability criterion for P in rank four? | Sharpened, not answered. Indecomposability plus odd determinant does not suffice (the determinant-121 Gaussian-integer lattice), and neither does indecomposability plus determinant not divisible by four for odd lattices (determinant 49). The surviving candidate is the minus-type conjecture below. |
| Does a restricted rank-four geometric family admit a general signed-shell argument? | Open. The minimal shell of a + 2L decides 90.4% of the 26,775 decisions for even lattices of odd determinant through 201, but the first nonzero index reaches norm4 = 132, so no uniform minimal-shell argument is available there. |
| What are the minimum orders of epsilon-one witnesses in the higher-order examples? | Partly answered. For even indecomposable rank-four lattices of odd determinant no involution is a witness (Theorem 3 of the parity report), so witnesses have 2-part of order 4 or 8; the determinant-121 witnesses have order 4 with square -I. Odd lattices and the degenerate even case always had involution witnesses in the lists. |
| Can the historical glue/eigenkernel claims be repaired? | Open in general. Theorem 3's eigenlattice and glue argument is a valid instance of that kind of reasoning, for involutions only. |
| Can stored CM matrices be given exact arithmetic provenance? | Answered for the examples at hand. The determinant-23,716 ambient lattice and the determinant-1156 sign kernel are trace forms of binary Hermitian lattices over the Eisenstein integers; the determinant-121 lattice is the trace form of a binary Hermitian lattice over the Gaussian integers. All three identifications are exact isometries. |
| Can the ambient order-three construction give a smaller rank-four counterexample? | Answered within the family: determinant 1156 is the minimum of the Eisenstein family F through D' = 154, and no level-2 configuration has sign-kernel determinant at most 24. The global minimum is 280 (release 1.2.0, [rank-four-minimum-determinant.md](rank-four-minimum-determinant.md)): the determinant-280 level-four lattice is the smallest counterexample, and every lattice of determinant at most 279 satisfies the converse. |

## Directions proposed by the assessment of 2026-10-21

| Question | Potential value | Smallest useful check | Promotion criterion | Status |
|---|---|---|---|---|
| Minus-type parity conjecture: every even, orthogonally indecomposable, positive definite rank-four lattice with det = 5 mod 8 has the parity property. | The sharpest surviving form of the parity question; a classification of symmetry witnesses would be a theorem even if identity zeros stayed open. | Stage 1: classify the remaining integral representations (order-4 actions faithful of order 4 on L/2L, order-8 actions) and their signs on characteristic stabilizers; involutions and square roots of -I are already excluded. Stage 2: decide whether every zero in the even, indecomposable, odd-determinant setting forces span(a, G^{-1} b mod 2) to be a totally singular plane; minus type has no such plane. | A witness-classification theorem or a certified counterexample to it; extend the census only to test an explicit obstruction found in stage 1. | Deferred, not started; evidence is all 131 distinct minus-type indecomposable classes through det 201 |
| Plus-type companion: under the same hypotheses with det = 1 mod 8, does the square-root-of-minus-I mechanism exhaust all zeros? | Would complete the even odd-determinant picture. | Examine the one failing indecomposable class (det 121) against the mechanism; search for a plus-type zero without a square root of -I. | A proof or a certified counterexample. | Deferred; second stage after the minus-type question |
| Level 3: a zero whose coset lattices are isometric while no isometry between them preserves the sign kernel. | Fills the one unexhibited level of the hierarchy. | Start from a lattice M, compare its fifteen index-two sublattices, look for isometric kernels in different Aut(M)-orbits, transport M along a kernel isometry, then check an integral third lattice and the full kernel group to exclude levels 1 and 2. | A certified level-3 configuration; a bounded failure does not prove impossibility. | Deferred; ranked below the minus-type question |
| Global determinant minimum for the rank-four symmetry converse. | Makes the bracket sharp. | A complete census of the four regions through 279 not covered by the earlier lists. | A complete census to 279 or a structural lower bound. | Answered in release 1.2.0: the minimum is 280 ([rank-four-minimum-determinant.md](rank-four-minimum-determinant.md)) |
| Independent human review of the core proofs and the attribution. | Required before any claim beyond "AI-assisted local review with computational replication". | Theorems 1 to 4 of the parity report, the level-four and determinant-1156 certificates, and the Conway-Sloane, Cervino-Hein, Schiemann and Earnest-Nipp attributions. | Acceptance by a human reviewer. | Not a computational task |

Rank five is not needed to establish the rank-four threshold and is not
proposed. Historical code and assertions are leads only. Add new entries when
an actual question appears; do not expand the list for its own sake.
