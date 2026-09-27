# Adversarial review of the rank-four census and its proof target

Prepared 2026-09-27 from `Grok-rank4-census-and-proof-target.md`. The feasibility review is the accepted prerequisite in `RANK4-REVIEW-RESOLUTION.md`. This note audits the determinant-`≤24` census and the Walsh reduction in `RANK4-PROOF-TARGET.md`. It does not reopen the rank-three classification, and it does not enlarge the determinant bound.

Earlier reports, manuscripts, census inputs, inventories, and ZIPs were left unchanged.

**Checked** means a certificate field recomputed in this pass, or a source read in this pass. **Deduction** means a step proved from those computations and from `docs/foundations-and-exact-enumeration.md` §§2 and 7, `docs/rank-two-classification.md`, and the accepted rank-three symmetry converse. The unrestricted rank-four symmetry converse remains open. Astra is not indicated by the obstruction below.

## Decisions

**(a) Census: accept.** Every positive definite integral rank-four lattice of determinant at most 24 occurs exactly once, up to integral isometry, among the 169 representatives. All 43,264 binary characteristics are resolved. Every zero has an epsilon-one witness, and every nonzero has a complete shell at norm4 at most 27. The search bound is 32. There is no converse counterexample in this domain.

**(b) Reduction: accept.** For each fixed Gram matrix, the symmetry converse is equivalent to the pair-specific claim that `Theta_G(P_A f_{a,b})` identically zero implies `P_A f_{a,b}=0`. The group-average criterion is unconditional and holds in every rank. It does not prove that claim in rank four.

**(c) Next task.** Prove the following statement, which is still conjectural. If a positive definite integral rank-four lattice is orthogonally indecomposable, and the stabilizer of a binary pair `(a,b)` has trivial epsilon, then `Theta_G(f_{a,b})` is not identically zero. The determinant-`≤24` census is a regression for that statement, not a proof. Any argument has to survive the determinant-15 cancellation recorded below: the minimal shell in the support of the average can vanish.

## 1. Census

Checked by `python code/verify_rank4_census.py data/census-rank4/det24-bound32.json` and the same command with `--backend ldl`. Both returned `verified` in this pass, in 35.379 and 35.900 seconds. Each report has

```
candidate_bases 1510
lattice_classes 169
even_pairs 22984
odd_pairs 20280
even_verdicts proved_nonzero 16608, proved_zero 6376
unresolved_pairs 0
symmetry_converse_in_domain true
```

The input SHA-256 of both runs is

```
f3411b2ae36cda48b1d93b92c8edd16eddf659e7a3031e47f613b484507ef982
```

That is the hash of `data/census-rank4/det24-bound32.json` (91,064,021 bytes), and it is the source hash recorded in `reports/census-rank4/det24-profile.json`. The profile’s 169 Gram matrices and summary match the certificate. `python scripts/verify_integrity.py

The arithmetic of the summary is `169 * 256 = 43264` characteristics, `136 * 169 = 22984` even pairs, and `120 * 169 = 20280` odd pairs. A direct pass over the certificate found exactly two verdicts: 26,656 zeros and 16,608 nonzeros. The 26,656 zeros are the 6,376 even symmetry zeros plus the 20,280 odd symmetry zeros. No pair has verdict `proved_zero_without_symmetry`, `proved_nonzero_by_modular_replay`, or `unresolved`. Odd pairs are symmetry zeros. Every stored zero has epsilon 1 in its stabilizer list, and the replay recomputes that sign from the full group.

### 1.1. Domain, transports, and uniqueness

The verifier regenerates the candidate list with `independent_candidates(24)` and requires the stored generation list to equal it. That is the accepted feasibility domain: 1,510 bases, before isometry deduplication. For each stored transport `T` it requires `det T = ±1` and `T^t G_rep T = G`. The representative Grams occur in order, each equal to its own first transported basis. Two representatives of equal determinant are rejected if the independent box search finds an isometry. Different determinants cannot be isometric. A finished run of that search is exhaustive: the integer box contains every vector of squared length at most the largest basis diagonal, by the same radius identity used in the feasibility review.

The determinant histogram in the certificate is

| Determinant | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Classes | 1 | 1 | 2 | 4 | 3 | 3 | 3 | 7 | 6 | 4 | 4 | 12 |

| Determinant | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Classes | 5 | 5 | 8 | 15 | 6 | 9 | 6 | 16 | 11 | 7 | 7 | 24 |

These 24 rows sum to 169. The determinant-`≤4` subtotal is 8, which is the pilot replayed by the tests.

### 1.2. Groups, stabilizers, and shells

For each representative the stored automorphism list is required to equal the box join. Stabilizer indices and epsilon signs are recomputed from that group by transporting `(a,b)`. A symmetry zero must name a matrix `T` in the list with `symmetry_data` epsilon 1. A nonzero must name the first shell whose signed coefficient is nonzero, and the stored shell histogram must equal a fresh enumeration through the run bound 32. The box enumeration and the LDL enumeration were both run. They use different vector searches and the same exact arithmetic. The box search does not use the producer’s LDL sphere.

The largest first nonzero norm4 in the certificate is 27. It occurs once, on

\[
G=\mathrm{diag}(1,1,1,24),\qquad a=(1,1,1,1),\qquad b=(0,0,0,0),
\]

with coefficient `+16`. Every coordinate is odd, so the squared length is at least `1+1+1+24=27`, and the 16 sign patterns of `(±1,±1,±1,±1)` all have `b=0` and contribute `+1`. This is an unsigned coset minimum inside the bound 32. The modular kernels, fifteen nonzero `b` on each representative, were recomputed by the replay. Their cutoffs run from 8 to 512. None of those cutoffs was used to decide a census verdict.

### 1.3. Interruptions and adversarial tests

Checked by `python scripts/run_tests.py --all`, eight tests, all passing in 6.696 seconds. A one-node generation raises before a certificate is written. A one-node replay returns `unresolved` and does not assert the domain converse. Building determinant 1 at shell bound 0 leaves unresolved pairs and does not assert the converse; relabeling such a pair as a symmetry-free zero is rejected. Deleting a candidate, a class, a group element, a stabilizer index, or flipping an epsilon sign is rejected. An invalid transport and a second representative of the same class are rejected. A synthetic modular fallback is accepted only when its cutoff certificate replays the same characteristic, and a tampered kernel is rejected. The completed determinant-24 certificate contains none of those fallback verdicts.

## 2. The Walsh reduction

### 2.1. The action sign

The functions `f_{a,b}` are defined on `(Z/4Z)^4` as in `RANK4-PROOF-TARGET.md`. For fixed `a` there are 16 residues `y ≡ a (mod 2)`, and the 16 choices of `b` are the Walsh characters of that coset. Distinct cosets have disjoint supports, so the 256 functions are an orthogonal basis of `R[(Z/4Z)^4]`. The value depends on `y` only through `(y-a)/2 (mod 2)`. Replacing `y` by `y+4k` changes that vector by `2k` and does not change the sign.

The action is `(T.f)(y)=f(T^{-1}y)`. If `T^t G T=G`, the change of variables `y=Tx` preserves `y^t G y`, so `Theta_G(T.f)=Theta_G(f)` and therefore `Theta_G(P_A f)=Theta_G(f)`.

Deduction of the sign. Suppose `T` preserves the labels: `Ta=a+2α` and `T^{-t}b=b+2β`. On the support,

\[
(T.f_{a,b})(y)=(-1)^{((y-Ta)/2)\cdot T^{-t}b}.
\]

The exponent is `((y-a)/2-α)·(b+2β)`. Modulo 2 this differs from `((y-a)/2)·b` by `α·b`. The pairing identity `(Ta-a)·b=a·(T^t b-b)` gives

\[
α·b=a·(T^t b-b)/2.
\]

From `T^{-t}b=b+2β`, applying `T^t` yields `(T^t b-b)/2=-T^t β`, so

\[
α·b=-(Ta)·β=-(a+2α)·β\equiv a·β\pmod 2.
\]

The right-hand side is the project’s epsilon, `a·(T^{-t}b-b)/2 (mod 2)`. Thus `T.f_{a,b}=(-1)^{epsilon(T)} f_{a,b}`.

Checked on every label-stabilizing automorphism of `I_4` (13,440 pairs) and of the order-1,152 determinant-4 class (13,824 pairs), and on all 1,200 label-stabilizing pairs in the order-24 group of the determinant-15 Gram in §3. No sign disagreed.

Outside the label stabilizer, `T.f_{a,b}` is `±f_{a',b'}` for the transported labels. Distinct basis functions are not scalar multiples of one another, so only the label stabilizer contributes to the coefficient of `f_{a,b}`.

### 2.2. Vanishing of the average

Deduction. Let `S` be the label stabilizer and let `v=sum_{g in A} g.f`. If some `σ` in `S` has `σ.f=-f`, then as `g` runs through `A` so does `gσ`, and

\[
v=\sum_g (gσ).f=\sum_g g.(σ.f)=-v.
\]

Hence `v=0` and `P_A f=0`. The identity `(gσ).f=g.(σ.f)` is the group action. If instead every element of `S` has epsilon 0, each contributes `+f`, and the coefficient of `f_{a,b}` in `P_A f` is `|S|/|A|`, which is positive because `S` contains the identity. Therefore `P_A f_{a,b}=0` if and only if `S` contains an epsilon-one element.

`docs/foundations-and-exact-enumeration.md` §2 already shows that epsilon is a homomorphism on the characteristic stabilizer and that an epsilon-one element forces the theta series to vanish. The calculation above identifies that sign with the Walsh action and gives both directions of the average criterion. Odd `a·b` is the case `epsilon(-I)≡1`. For `b=0` every sign is `+1`, and for `a=0` the constant term is `+1`, so epsilon is identically 0 and the series does not vanish. The pairs that remain are even, with `a` and `b` both nonzero.

The census coefficient of `q^{y^t G y}` in `Theta_G(f_{a,b})` is the signed norm4 coefficient. The symmetry converse for this Gram matrix is exactly the assertion that `Theta_G(P_A f_{a,b})` identically zero implies `P_A f_{a,b}=0`. The note’s separation is correct: this equivalence does not assert the nonvanishing, and it does not require `Theta_G` to be injective on every invariant function.

## 3. One proved reduction, and one failed shortcut

### 3.1. Orthogonal sums

Deduction. Let `L=L_1 ⊥ L_2` be an orthogonal sum of positive definite integral lattices of positive ranks adding to 4. In a splitting basis, `G=diag(G_1,G_2)` and `f_{a,b}=f_{a_1,b_1}⊗f_{a_2,b_2}`. Absolute convergence gives

\[
\Theta_G(f_{a,b})=\Theta_{G_1}(f_{a_1,b_1})\,\Theta_{G_2}(f_{a_2,b_2}).
\]

`docs/foundations-and-exact-enumeration.md` §7 records the same product for the project’s theta series of component characteristics. The ring of holomorphic functions on the upper half-plane is an integral domain, so a vanishing product has a vanishing factor. Ranks 1 and 2 have the symmetry converse at every determinant (`docs/rank-two-classification.md`; rank 1 is the constant-term and `-I` analysis in `docs/foundations-and-exact-enumeration.md` §§2 and 8). Rank 3 has the symmetry converse at every determinant (the accepted uniform theorem). The vanishing factor therefore has an epsilon-one automorphism `T_i`. Extending it by the identity on the complementary summand preserves the block-diagonal Gram matrix. The dual shift of `b` is supported on that factor, so the epsilon on `(a,b)` equals the factor’s epsilon and is 1.

A lattice that becomes block diagonal only after a change of basis `U∈GL(4,Z)` is the same case after the audited transport `a ↦ U^{-1}a (mod 2)`, `b ↦ U^t b (mod 2)`. The witness on the splitting basis conjugates back by `U`. A sum of more than two factors is the same argument applied to one vanishing factor.

This proves the symmetry converse for every orthogonally decomposable positive definite integral rank-four lattice, at every determinant. The determinant-24 census is not an input.

Checked decomposition of the 169 representatives, by searching vectors of squared length at most 24. A primitive vector `v` with `Q(v)` dividing every entry of `Gv` spans a line summand. A saturated pair whose `2×2` Gram determinant is positive and divides the pairings with the ambient basis spans a plane summand. The bound 24 is large enough for a reduced basis of either summand: the summand determinant divides the lattice determinant, and a reduced binary diagonal is at most that determinant. The search reports 142 decomposable lattices (138 with a line summand and 4 with a plane summand) and 27 indecomposable lattices. It finds a line in `I_4` and in `diag(1,1,1,24)`, a plane in `A_2⊕A_2`, and no summand in the standard Grams of `D_4` and `A_4`. The 27 indecomposable classes are exactly the classes whose full groups have orders

```
8, 12, 12, 12, 16, 16, 16, 24, 24, 24, 24, 24, 32, 48, 48, 48, 48, 64, 96, 96, 96, 96, 96, 96, 96, 240, 1152.
```

None has order 2. They carry 73 of the even symmetry zeros. Those 73 zeros have epsilon-one witnesses, so their averages vanish. `D_4` and `A_4` are in this set. The product argument does not apply to them.

### 3.2. The minimal shell of the average can cancel

The natural shortcut is: if `P_A f≠0`, the minimal norm attained on its support has a nonzero shell sum, because one orbit contributes a nonzero weight. That shortcut is false.

Checked example. Representative 131,

\[
G=\begin{pmatrix}2&-1&-1&-1\\ -1&2&0&0\\ -1&0&3&0\\ -1&0&0&3\end{pmatrix},
\]

has determinant 15 and automorphism group of order 24. The pair `a=(0,0,1,1)`, `b=(0,1,1,1)` is even. Its stabilizer has order 8 and every epsilon is 0, so the coefficient of `f_{a,b}` in its average is `8/24=1/3`. The census witness is the norm4-14 coefficient `-2`. The same average was recomputed on all 256 residues. Sixteen residues have nonzero average. They form four orbits, with constant values `1/3`, `-1`, `-1`, and `1/3` after dividing the group sum by 24.

At norm 6 the support consists of eight vectors:

| Vectors | Residue orbit | Value of `P_A f` | Contribution |
|---|---|---:|---:|
| `(0,0,1,1)`, `(0,0,-1,-1)`, `(2,0,1,1)`, `(-2,0,-1,-1)`, `(2,2,1,1)`, `(-2,-2,-1,-1)` | six residues | `1/3` | `+2` |
| `(0,0,1,-1)`, `(0,0,-1,1)` | two residues | `-1` | `-2` |

The norm4-6 coefficient is `2+(-2)=0`. Two further orbits meet at norm 14. Their multiplicity-weighted sum is `-2`, in agreement with the stored shell. The first nonzero coefficient is produced by unequal multiplicities of two orbits at the same norm, after an earlier exact cancellation.

The same scan found six epsilon-trivial pairs on indecomposable representatives for which the minimal positive norm in the support of `P_A f` has shell coefficient 0. One fully expanded pair is enough to retire the shortcut. The census coefficient at the later norm is still nonzero, so these six pairs are not converse counterexamples.

## 4. What remains

The unconditional reduction and the orthogonal-sum lemma leave one conjecture.

**Conjecture.** Let `G` be the Gram matrix of an orthogonally indecomposable positive definite integral lattice of rank four, and let `a,b` be binary. If every automorphism that stabilizes `(a,b)` has epsilon 0, then `Theta_G(f_{a,b})` is not identically zero.

The gap is a reason for some norm4 coefficient to be nonzero after several orbits of `P_A f` are allowed to share a norm. The determinant-15 Gram is the smallest explicit regression: a proposed identity has to return the coefficient 0 at norm 6 and a nonzero coefficient by norm 14 on `a=(0,0,1,1)`, `b=(0,1,1,1)`. The determinant-`≤24` census already confirms the conjecture for the 27 indecomposable classes, including the six pairs whose minimal support shell cancels. That confirmation uses the shell search through bound 32. It does not supply the identity.

No consultation is indicated. The cancelled minimal shell is a failed shortcut inside the main argument, and the conjecture above is the main rank-four claim restricted to indecomposable lattices.

## 5. What this pass ran and read

Ran: both coefficient backends of `code/verify_rank4_census.py` on `data/census-rank4/det24-bound32.json`; `python scripts/run_tests.py --all`; `verify_manuscript_manifest.py` on `reports/rank4-census-artifacts.json`; a direct pass over all 43,264 stored verdicts, shell indices, and modular cutoffs; the Walsh sign comparison on `I_4`, the order-1,152 class, and the determinant-15 group; the orthogonal-summand search on all 169 representatives and on `I_4`, `D_4`, `A_4`, `A_2⊕A_2`, and `diag(1,1,1,24)`; the support-orbit calculation for the determinant-15 pair.

Read: `Grok-rank4-census-and-proof-target.md`; `RANK4-CENSUS-STATUS.md`; `RANK4-PROOF-TARGET.md`; `RANK4-REVIEW-RESOLUTION.md`; `code/rank4_census.py`; `code/verify_rank4_census.py`; `tests/test_rank4_census.py`; the two saved final replay summaries; `docs/foundations-and-exact-enumeration.md` §§2, 7, and 8; `docs/rank-two-classification.md`’s statement of the rank-one and rank-two symmetry converse; the orthogonal-sum argument in `reviews/rank-three-proof-review.md` §6. The 91 MB certificate was parsed programmatically and was not copied into this note.
