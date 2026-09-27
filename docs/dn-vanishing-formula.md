# Uniform Dn vanishing criterion and counting formula

Completed 2026-09-26. The theorem below holds for every n >= 2. The saved finite audit covers D2 through D10; it checks every characteristic class, with explicit proofs for all even zeros and complete coefficients for all even nonzeros. This note uses the theta and sign-character conventions of [docs/foundations-and-exact-enumeration.md](foundations-and-exact-enumeration.md).

## The theorem

For the standard Euclidean lattice

\[
D_n=\{x\in\mathbb Z^n:\sum_i x_i\equiv0\pmod2\},
\qquad D_n^*=\mathbb Z^n\cup(\mathbb Z^n+\tfrac12\mathbf1),
\]

the number of even-p characteristic pairs with identically zero theta is

\[
\boxed{Z_n=4^{n-1}-3^n+3\,2^{n-1}-1
       +\mathbf1_{2\mid n}\frac12\binom n{n/2}.}
\]

Pairs mean classes of xi modulo Dn and delta modulo Dn*, with 2xi in Dn and 2delta in Dn*. There are 4^n pairs in total and 2^(2n-1)+2^(n-1) even pairs. Every even zero has an order-two signed-permutation witness with epsilon=1. Every nonzero characteristic has a nonzero minimal-shell coefficient, at squared norm at most max(1, floor(n/2)/2).

Consequently D2 has one even zero, D3 has none, and every Dn with n >= 4 fails the parity converse. This result supplies the complete Dn counts, not a classification of all integral lattices.

## Complete characteristic parameters

Write A=2xi and B=2delta in ambient coordinates. Every primal class has a unique representative

\[
A=\mathbf1_S+2e e_1,\qquad |S|=s=2k,\qquad e\in\{0,1\}.
\]

Indeed S is the set of odd coordinates of A. After subtracting its indicator, the resulting even vector divided by two has a well-defined sum modulo two; call this e. Two vectors have the same S and e exactly when their difference belongs to 2Dn. There are 2^(n-1) even subsets S, hence 2^n primal classes.

Every dual class has a representative

\[
B=w+\frac h2\mathbf1,\qquad h\in\{0,1\},\qquad w\in\{0,1\}^n.
\]

For fixed h, w is taken modulo simultaneous complementation. This follows because 2Dn* consists of integral vectors whose coordinates are either all even or all odd. Choose w_n=0 to represent each class once. There are 2^n dual classes. Set T=supp(w), t=|T| and r=|S intersect T|. Then

\[
p=A\cdot B=r+2e w_1+h(k+e)
\equiv r+h(k+e)\pmod2.
\]

Complementing w changes neither parity nor any of the vanishing tests below: s is even, r changes to s-r, and t changes to n-t.

## Parity projection, with its phase

Put q=exp(pi i tau) and define the one-dimensional function

\[
\vartheta[\alpha,\beta]
=\sum_{m\in\mathbb Z}q^{(m+\alpha)^2}
                       e^{2\pi i(m+\alpha)\beta}.
\]

Absolute convergence for Im(tau)>0 permits the parity projection onto Dn:

\[
\Theta_{D_n}[\xi,\delta]
=\frac12\left(P+(-1)^{k+e}Q\right),\qquad
P=\prod_i\vartheta[\xi_i,\delta_i],\quad
Q=\prod_i\vartheta[\xi_i,\delta_i+\tfrac12].
\]

The sign comes from exp(-pi i sum(xi_i))=(-1)^(k+e) when replacing (-1)^sum(m_i) by a shift of delta. Integer shifts of a coordinate xi_i leave its individual factor unchanged, but e must remain in this relative sign. Erasing e by reducing xi coordinatewise modulo integers would lose the balanced cancellations.

## Integer B: product zeros and balanced cancellations

First let h=0. In one dimension, vartheta[1/2,1/2] is identically zero by pairing m and -m-1; the other three binary factors are nonzero by their leading coefficients.

If s>0, P has a zero factor exactly when r>0, and Q has a zero factor exactly when r<s. Thus both products vanish exactly when 0<r<s. If r=0 or r=s, exactly one survives; its leading coefficient has absolute value 2^s, so Theta has a nonzero coefficient of absolute value 2^(s-1) at squared norm s/4. These endpoint classes have even p. Restricting to even p in the zero case requires r even.

If s=0 and e=0, the coset contains zero, with constant coefficient 1. If s=0 and e=1, set

\[
\vartheta_3=\vartheta[0,0]=1+2q+O(q^4),\qquad
\vartheta_4=\vartheta[0,\tfrac12]=1-2q+O(q^4).
\]

Then

\[
\Theta=\frac12\left(\vartheta_3^{n-t}\vartheta_4^t
                     -\vartheta_4^{n-t}\vartheta_3^t\right).
\]

This is identically zero when t=n/2. Otherwise its coefficient at q is 2(n-2t), which is nonzero. Squared norm 1 is the minimum of this coset: it consists of integer vectors of odd coordinate sum. Hence these additional zeros occur exactly when n is even and t=n/2.

## Half-integer B: no even zeros

Now let h=1. Put g=vartheta[0,1/4] and f=vartheta[1/2,1/4]. Reindexing m or pairing opposite terms gives

\[
\begin{aligned}
\vartheta[0,\tfrac34]&=\vartheta[0,\tfrac54]=g,
&g&=1+O(q),\\
\vartheta[\tfrac12,\tfrac34]&=\vartheta[\tfrac12,\tfrac54]=-f,
&f&=\sqrt2\,q^{1/4}+\text{higher terms}.
\end{aligned}
\]

In passing from P to Q, a coordinate on S supplies a minus sign precisely when w_i=0; a coordinate outside S supplies none. Thus Q=(-1)^(s-r)P=(-1)^r P and

\[
\Theta=\frac12\bigl(1+(-1)^{k+e+r}\bigr)P
       =\frac12\bigl(1+(-1)^p\bigr)P.
\]

For even p this equals P, whose leading coefficient has absolute value (sqrt(2))^s=2^k at squared norm s/4. It cannot vanish. For s=0 even parity forces e=0, and the coefficient is 1. This also checks the phase in the quarter-shift sector directly.

Combining the two sectors, the full criterion is: Theta vanishes if p is odd, or if p is even and h=0 with either

1. s>0 and 0<r<s; or
2. s=0, e=1, n even and t=n/2.

All other cases have the nonzero leading coefficients given above. The bounds and minimal-shell assertion in the theorem follow. A minimal-shell test is valid here because of this proof; it is not a valid universal stopping rule for arbitrary lattices.

## Exact counting

Fix a nonempty even set S of size s. There are 2^(s-1)-2 choices of an even intersection with T that is neither empty nor all of S. There are independently 2^(n-s) choices outside S. Divide by two for complementation, and multiply by two for e. Therefore the number of product-zero even classes is

\[
\sum_{\substack{2\le s\le n\\s\text{ even}}}
\binom ns(2^{s-1}-2)2^{n-s}
=4^{n-1}-3^n+3\,2^{n-1}-1.
\]

For the simplification, use sum_even C(n,s)=2^(n-1) and
sum_even C(n,s)2^(n-s)=(3^n+1)/2, subtracting the s=0 term in each sum. The balanced classes have S empty and e=1, with C(n,n/2)/2 choices of T modulo complementation when n is even. They are disjoint from the product zeros. Adding them proves the formula for every n >= 2.

## Uniform symmetry certificates

For a product zero choose i in S intersect T and j in S outside T, and flip those two coordinates. This is an involutive isometry of Dn. It shifts xi by -A_i e_i-A_j e_j, an integral vector of even sum; it shifts delta by -B_i e_i-B_j e_j, an integral vector in Dn*. Its sign pairing is

\[
\epsilon=-A_iB_i-A_jB_j\equiv1\pmod2,
\]

since A_i,A_j are odd, B_i is odd and B_j is even.

For a balanced zero, pair T with its complement and swap each pair. This involution preserves Dn. Here xi is an integer vector of odd sum. Its shift has sum zero and belongs to Dn; the shift of delta has every coordinate congruent to 1/2 modulo integers and belongs to Dn*. Pairing it with 2xi is congruent to sum(xi_i), hence odd. The sign-character theorem proves vanishing. These constructions explain every even zero without computing the full automorphism group.

## Exact bases and independent finite audit

The code uses the column basis M=(e1-e2, ..., e_(n-1)-e_n, e_(n-1)+e_n), with G=M^t M and det(G)=4. It generates Dn: for any integral even-sum vector v, its coordinates y satisfy

\[
y_i=\sum_{j=1}^i v_j\ (1\le i\le n-2),\qquad
y_n=\tfrac12\sum_jv_j,\qquad y_{n-1}=y_n-v_n.
\]

The displayed prefix is empty for n=2. For n=2, G=diag(2,2), identifying D2 with A1 perpendicular A1. For n=3,

\[
G=\begin{pmatrix}2&-1&-1\\-1&2&0\\-1&0&2\end{pmatrix},\qquad
U=\begin{pmatrix}0&1&0\\1&0&0\\0&0&1\end{pmatrix},\qquad
U^tGU=\begin{pmatrix}2&-1&0\\-1&2&-1\\0&-1&2\end{pmatrix},
\]

which identifies D3 with A3 exactly.

The dual description used above also follows directly from this basis: pairing with e_i-e_j forces all dual coordinates to have the same fractional part, and pairing with e_(n-1)+e_n forces twice that part to be integral. Conversely, every integral vector and every integral vector plus half the all-ones vector pairs integrally with Dn.

[code/dn_exact.py](../code/dn_exact.py) implements the criterion and integral witness checks. [code/dn_ambient.py](../code/dn_ambient.py) independently enumerates the complete sphere of ambient integer vectors v with even sum and norm <= norm4_bound. They represent 2(x+xi), so the squared norm of v is four times the theta exponent. Recover y=M^-1 v, set a=y mod 2 and x=(y-a)/2, and histogram the parity of x. A second backend uses rational LDL shell recursion in the lattice basis. Both apply the same exact Walsh transform to their independently computed histograms; the coefficient summation is shared, the sphere enumeration is independent.

The sufficient four-times-norm bound is max(4, 2 floor(n/2)). For each even pair in strict binary-mask order the saved stream contains either a complete nonzero coefficient or an involution recipe. Missing, repeated, extra and corrupted records are rejected. A zero record must pass the actual stabilizer and epsilon equations, rather than merely satisfy the predicted criterion. With A=2xi and C=4delta, these checks are:

- PA-A has even entries and coordinate sum divisible by four;
- PC-C has all entries congruent to zero modulo four, or all congruent to two;
- A dot (PC-C) is divisible by four and its quotient is odd.

Signed permutations preserve Dn and the Euclidean form by construction. The recipes are involutions. [expand_record](../code/dn_exact.py) converts them to the original verifier's matrix T=M^-1 P M, with exact integrality and order checked in the sample tests. No vanishing verdict rests solely on finitely many zero coefficients.

| Rank | Even pairs | Product zeros | Balanced zeros | Total even zeros | Even nonzeros |
|---|---:|---:|---:|---:|---:|
| D2 | 10 | 0 | 1 | 1 | 9 |
| D3 | 36 | 0 | 0 | 0 | 36 |
| D4 | 136 | 6 | 3 | 9 | 127 |
| D5 | 528 | 60 | 0 | 60 | 468 |
| D6 | 2080 | 390 | 10 | 400 | 1680 |
| D7 | 8256 | 2100 | 0 | 2100 | 6156 |
| D8 | 32896 | 10206 | 35 | 10241 | 22655 |
| D9 | 131328 | 46620 | 0 | 46620 | 84708 |
| D10 | 524800 | 204630 | 126 | 204756 | 320044 |

All 700,070 even pairs pass both backends, with identical record and complete-histogram hashes. The remaining 698,026 odd pairs vanish by parity; their computed coefficients also vanish within the audited bounds. There are no unresolved pairs. The primary D10 build took about 21.1 seconds, its ambient replay about 4.3 seconds, and its compressed stream occupies 1,340,336 bytes on this run. These timings describe this machine and run, not an asymptotic guarantee.

Certificate streams and standalone examples are under [data/dn](../data/dn). [dn-audit-summary.json](../reports/dn-audit-summary.json) records the primary results; each `reports/dn-dN-ambient.json` records the independent replay, including coverage and hashes.

## Replay and tests

Python 3.11+, standard library only:

```powershell
python code/build_dn.py 2 3 4 5 6 7 8 9 10
foreach ($dnRank in 2..10) {
    python code/verify_dn.py "data/dn/d$dnRank.jsonl.gz" --backend ambient --output "reports/dn-d$dnRank-ambient.json"
}
python scripts/run_tests.py --all
```

The combined suite passes all 42 tests. Dn tests replay every saved stream, compare enumeration hashes, verify quotient coverage, check the formula against its unsimplified sum through n=60, check sector coefficients, export witnesses to the original verifier, and reject incomplete or corrupted coverage. A limited LDL enumeration returns unresolved, without publishing a partial classification. The finite certificate reader deliberately accepts ranks 2 through 10; the mathematical proof is uniform in n.

The D4 test also enumerates all 384 signed permutations and finds exactly 64 epsilon-one witnesses for xi=e1, delta=(0,0,1/2,1/2). It rejects the alternative flip of coordinates 1 and 3 for that class because it does not stabilize delta. This verifies that subgroup count, not the separate full Aut(D4) order claim.

## Relation to the completed results

The universal implication "even p implies nonzero theta" is false; D4 is an indecomposable counterexample. The original Q1 also asks which lattices satisfy the implication and about its literature history. Those broader requests are not answered simply by saying "no." Within irreducible ADE lattices the answer is exactly An and E6, with all ADE even-zero counts now justified. The unrestricted classification and literature audit remain separate tasks, as does the symmetry converse outside the established families.

The subsequent Bn,Cn,F4,G2 milestone is complete in [docs/root-lattice-classification.md](root-lattice-classification.md), with integral scaling and lattice conventions. The completed manuscripts include the lower-rank classification and the rank-four counterexample. Minimum witness-order questions remain deferred.
