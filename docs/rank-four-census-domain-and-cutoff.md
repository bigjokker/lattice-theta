# Rank-four feasibility: finite generation and an identity cutoff

Complete finite candidate generation and inclusive rank-four modular coefficient cutoff. The completed determinant-24 census and counterexample are described in the main manuscript.

## 1. A complete finite candidate domain

Given any positive definite integral rank-four lattice, choose a shortest
primitive vector, project along it, choose a shortest primitive vector in
the projected lattice, and repeat. Each projected lattice is discrete: the
projections of a completion basis are linearly independent and generate it.
Each primitive vector extends to an integral basis. Lifting those successive
bases produces a lattice basis v1,...,v4. Subtract earlier basis vectors in
reverse order to make all Gram-Schmidt coefficients mu_ij have absolute
value at most 1/2. This preserves the projected vectors and basis saturation.

Write d_i=||v_i^*||^2 for the positive Gram-Schmidt squared lengths and
Delta=product d_i=det G. In the projected lattice at stage i, the projection
of v_(i+1) has norm d_(i+1)+mu_(i+1,i)^2 d_i and cannot be shorter than
the selected shortest vector of norm d_i. Therefore

    d_(i+1) >= (3/4)d_i.

In particular, A=d_1 is a positive integer and

    Delta >= (3/4)^6 A^4 = (729/4096)A^4.

For Delta<=H this bounds A; A<=3 when H=24. For a chosen prefix of length
k-1, put P=product_(j<k)d_j and r=5-k. The remaining determinant satisfies

    Delta >= P d_k^r (3/4)^(r(r-1)/2).

This gives a finite upper bound for each next pivot and a lower bound
d_k>0, d_k>=(3/4)d_(k-1). Every original basis vector also has norm >=A.

To enumerate integral Gram entries, choose each cross entry G_jk successively
for j=1,...,k-1. Its Gram-Schmidt relation is

    mu_kj = (G_jk - sum_(l<j) mu_jl mu_kl d_l)/d_j.

Hence it lies in the finite integral interval centered at that sum with
radius d_j/2. After the cross entries are fixed, put
q=sum_(j<k) mu_kj^2 d_j and d_k=G_kk-q. Enumerate integral G_kk>=A
subject to the positive-pivot, adjacent-pivot and determinant inequalities.
For the fourth diagonal, solve P(G_44-q)=Delta for Delta=1,...,H and
retain only integral G_44 with the same bounds. Positive pivots prove
definiteness. Everything uses rational arithmetic, with no floating roots.

These necessary inequalities may admit extra bases whose first projected
vectors are not shortest. Extra bases are allowed. The proof asserts that
every lattice has at least one admitted basis, not that every admitted basis
is reduced or that candidates are canonical isometry representatives.

For the independent diagonal scanner, d_i>=(3/4)^(i-1)A>=(3/4)^(i-1).
Bounding the other three pivots below gives d_i<=H/(3/4)^(7-i).
Consequently G_ii<=U_i+(1/4)sum_(j<i)U_j with U_i=H/(3/4)^(7-i).
Cross-entry boxes follow from
|G_jk|<=sum_(l<=j)|mu_jl|d_l/2. The scanner uses these loose boxes,
Schur projections and direct diagonal scanning, rather than the builder's
centered cross intervals and last-diagonal determinant solving.

A generator node limit must raise and discard its partial result. Candidate
generation alone gives neither a count of isometry classes nor a converse
result. Complete basis-image isometry tests and independent reconstruction
of the declared candidate domain are still required for a future census.

## 2. Inclusive rank-four identity cutoff

Let G be any positive definite integral 4x4 Gram matrix, a,b binary columns,
b nonzero and a^t b even. Set L0={x:x^t b even}, let G0 be its Gram matrix,
d0=det G0, and let N0 be the least positive inverse-denominator-clearing
integer for G0. Set

    M=16N0,
    mu=[SL2(Z):Gamma0(M)] = product_(p^e || M) p^(e-1)(p+1),
    B=mu/6.

For complete signed norm4 coefficients c_N, cancellation at every index
0<=N<=B implies an identity zero. The constant term and endpoint are included.
If an enumeration is interrupted, none of its partial cancellation is used.
The statement has no determinant cap and does not assert optimal level.

Proof. The checked common-space translation in `docs/coset-modularity.md` gives
F(z)=sum c_N exp(2 pi i Nz), Theta(tau)=i^(a^t b) F(tau/8). In rank four,
[Kane-Kim, Proposition 2.3](https://arxiv.org/html/2211.03987v2) gives weight
two and quadratic character chi_(4d0) on Gamma0(M). Both doubled cosets
belong to that same space, including a trivial coset based on 2L0.
Thus F^2 is holomorphic of weight four with trivial character, at all cusps
as well as in the upper half-plane. No cuspidality is needed.

On Gamma0(M), infinity has width one and -I belongs to the group.
[Brunault, Theorem 1](https://perso.ens-lyon.fr/francois.brunault/recherche/Sturm-bound-general.pdf)
implies ord_infinity(F^2)<=mu/3 for a nonzero square. If F has first nonzero
coefficient at m, its square starts at 2m with nonzero coefficient c_m^2.
Therefore m<=mu/6. Since M is divisible by 16, its two-primary index factor
is 3 times a multiple of eight, so mu is divisible by six and B is integral.
Vanishing through B would force m>=B+1, a contradiction. This is a complex
identity test, not a finite-field congruence theorem. In this proof the weight,
the power and the bound are derived afresh; the ternary constant is not reused.

## 4. Completed local feasibility checks

The new implementation is `code/rank4_feasibility.py`; the independent checker is
`code/verify_rank4_feasibility.py`. The latter uses inverse-Gram integer boxes
and compatible-column-pair joining for full groups, and direct transported
characteristics for every stabilizer. It does not use the builder's action
tables or LDL sphere enumeration. The candidate scanner uses the distinct
loops described in section 1. Both algorithms agree for every determinant
bound from 1 through 24, including all 1,510 candidate bases at 24.

| Determinant bound | Candidate bases, before isometry deduplication |
|---|---:|
| 1 | 1 |
| 2 | 2 |
| 4 | 47 |
| 8 | 166 |
| 12 | 406 |
| 16 | 688 |
| 24 | 1,510 |

The saved candidate domain and sizing are in
`data/rank4-feasibility/candidates-det24.json` and `candidate-sizing.json`.
These counts do not count lattice classes. The later census is documented
in the rank-four manuscript and `reports/census-rank4/det24-profile.json`.

All five named controls replay independently, including all 1,280 binary
characteristic pairs and complete even and odd stabilizers:

| Control | Full group order | Even symmetry zeros | Even nonzeros | Odd zeros |
|---|---:|---:|---:|---:|
| I4 | 384 | 55 | 81 | 120 |
| 2I4 | 384 | 55 | 81 | 120 |
| A4 | 240 | 0 | 136 | 120 |
| D4 | 1,152 | 9 | 127 | 120 |
| A3 perpendicular a norm-one line | 96 | 28 | 108 | 120 |

The input is `data/rank4-feasibility/controls.json`; the independent
saved report is `reports/rank4-feasibility/independent-box.json`. These are
five controls, not a representative list of all rank-four lattices.

Ten separate cutoff controls are in `data/rank4-feasibility/cutoff-controls.json`.
LDL construction and independent box replay agree on their complete data.
Three limited runs remain unresolved: a missing endpoint, an interruption,
and a cancelled first shell below the cutoff. The last is a ternary regression
extended by a norm-one line; its norm4=6 coefficient cancels and its complete
norm4=10 coefficient is two. Its inclusive cutoff is 384, so cancellation
at six alone is not an identity test. The default kernel b=(1,0,0,0) gives
cutoffs 16,32,48,16 for I4,2I4,A4,D4 respectively.

All eight tests in `tests/test_rank4_feasibility.py` pass. They check domain
agreement, D4 coverage, noncanonical bases, saturated kernel changes,
independent index arithmetic, first-shell cancellation, endpoint omissions,
group/coefficient tampering and interruption policy. The complete suite is run by `python scripts/run_tests.py --all`. 

Print-only independent replay and tests:

```powershell
python code/verify_rank4_feasibility.py data/rank4-feasibility
python scripts/run_tests.py --all
```

Builders and verifier output options write files; use deliberate fresh
destinations for a new run. The saved report binds the exact candidate,
characteristic-control and cutoff-control pack hashes.

