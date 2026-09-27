# Exact foundations for lattice theta functions

Foundational proofs and exact enumeration conventions. Subsequent uniform rank-three results and the rank-four counterexample are in the companion manuscripts.

## 1. Setup and characteristic classes

Let L be a positive definite integral lattice of positive rank in its real span, with dual L*. For Im(tau)>0 define

\[
\Theta_L[\xi,\delta](\tau)=\sum_{x\in L}
 e^{\pi i\tau\|x+\xi\|^2}e^{2\pi i\langle x+\xi,\delta\rangle},
\quad a=2\xi\in L,\quad b=2\delta\in L^*,\quad p=\langle a,b\rangle.
\]

Positive definiteness implies absolute convergence, locally uniform in the upper half plane, and finitely many vectors in each bounded norm range. Thus reindexing the sum is justified. Shifting xi by ell in L leaves Theta unchanged; shifting delta by nu in L* multiplies it by exp(2 pi i <xi,nu>). Consequently vanishing depends only on the two characteristic classes. Replacing a,b by a+2ell,b+2nu changes p by an even integer, so its parity also depends only on the classes.

## 2. The stabilizer sign-character theorem

Let

\[
H=\{\sigma\in\operatorname{Aut}(L):\sigma\xi-\xi\in L,
\ \sigma\delta-\delta\in L^*\},\qquad
\varepsilon(\sigma)=\langle a,\sigma\delta-\delta\rangle\pmod2.
\]

**Theorem.** H is a subgroup, epsilon is a well-defined homomorphism H -> F2 depending only on the characteristic classes, and

\[
\varepsilon\not\equiv0\implies\Theta_L[\xi,\delta]\equiv0.
\]

**Proof.** H is the stabilizer of the ordered pair of classes under the induced group action, hence a subgroup. Put mu_sigma = sigma delta - delta in L*. The pairing <a,mu_sigma> is an integer.

Changing xi to xi+ell changes epsilon by <2ell,mu_sigma>, which is even. Changing delta to delta+nu changes it by

\[
\langle a,\sigma\nu-\nu\rangle
=\langle\sigma^{-1}a-a,\nu\rangle\in2\mathbb Z,
\]

since sigma stabilizes xi and therefore sigma^{-1}a-a is in 2L. These changes also leave H unchanged. This proves independence of representatives, including simultaneous changes.

For sigma,rho in H,

\[
\mu_{\sigma\rho}=\sigma\mu_\rho+\mu_\sigma,
\qquad
\langle a,\sigma\mu_\rho\rangle
=\langle\sigma^{-1}a,\mu_\rho\rangle
\equiv\langle a,\mu_\rho\rangle\pmod2.
\]

Hence epsilon(sigma rho)=epsilon(sigma)+epsilon(rho). Finally orthogonal change of variables and the two shift laws give

\[
\Theta_L[\xi,\delta]
=\Theta_L[\sigma\xi,\sigma\delta]
=\Theta_L[\xi,\delta+\mu_\sigma]
=e^{2\pi i\langle\xi,\mu_\sigma\rangle}\Theta_L[\xi,\delta]
=(-1)^{\varepsilon(\sigma)}\Theta_L[\xi,\delta].
\]

An element with epsilon=1 therefore forces the entire function to be zero. This proves the theorem.

**Consequences.** -I always belongs to H and epsilon(-I)=-p modulo 2, recovering classical odd-parity vanishing. Epsilon vanishes on every odd-order element, because 0=epsilon(sigma^k)=k epsilon(sigma). Epsilon(sigma^{-1})=epsilon(sigma), and epsilon vanishes on every square in H. None of these assertions proves the converse of the theorem.

**Symmetry converse and its subsequent resolution.** For even p, does Theta identically zero imply that epsilon is nontrivial on H? An example with Aut(L)={+I,-I} would disprove this implication, since epsilon is then zero for even p. This small-group condition is a useful search target, not a requirement: a larger automorphism group with trivial epsilon on the relevant full stabilizer would also give a counterexample. The rank-four companion now gives precisely such a small-group counterexample. The classification of lattices with no even-p zeros is a separate question.

Section 15 proves the symmetry converse for every characteristic in ranks one
and two at every determinant; the remaining question starts at rank three.

## 3. Coordinates and the inverse-transpose convention

Choose a lattice basis e_i and its dual basis e_i*. From this point a,b are the integer coordinate columns of 2xi,2delta in those respective bases. Let G=(<e_i,e_j>) and let T act on primal columns. Then

\[
T^tGT=G,\quad p=a^tb,\quad
(T-I)a\in2\mathbb Z^n,\quad(T^{-t}-I)b\in2\mathbb Z^n,
\quad\varepsilon(T)=a^t\frac{(T^{-t}-I)b}{2}\pmod2.
\]

The dual action is T^{-t}: invariance of the pairing requires (Ta)^t b'=a^t b for all a. The division by two is permitted only after verifying the dual stabilizer condition. An integral T preserving G has determinant +/-1 (take determinants), so its rational inverse must have integral entries; the verifier checks this directly too.

The alternative expression a^t(T^t-I)b/2 is epsilon(T^{-1}). If T is in H, so is T^{-1}, and the homomorphism proof gives equality modulo two. Thus the legacy expression is valid on H, after both stabilizer conditions have been verified. Outside H it supplies no certificate.

[docs/theta-transformations.md](theta-transformations.md) identifies this sign with
the reduction of characteristics under the standard GL(n,Z) theta action,
using the checked DLMF characteristic-shift and basis-change formulas. The
unreduced basis-change identity has no extra prefactor, even for determinant -1.

## 4. Signed coefficients and complete shell enumeration

In primal coordinates x in Z^n, put y=2x+a. The phase and norm are

\[
e^{2\pi i\langle x+\xi,\delta\rangle}=i^p(-1)^{x^tb},
\qquad 4\|x+\xi\|^2=y^tGy.
\]

Thus, for the integer N>=0, the signed coefficient is

\[
c_N=\sum_{\substack{y\in\mathbb Z^n,\ y\equiv a\ (2)\\y^tGy=N}}
(-1)^{((y-a)/2)^tb},\qquad
\Theta=i^p\sum_{N\ge0}c_N e^{\pi i\tau N/4}.
\]

If one complete c_N is nonzero, Theta is not identically zero. Indeed, if any coefficient is nonzero, there is a least such N; along tau=it, multiply by exp(pi t N/4) and let t tend to infinity. Absolute convergence at any fixed positive t justifies domination of the remaining tail, and the limit is i^p c_N. A finite run of zero coefficients gives no conclusion beyond its range.

**Enumeration completeness.** Rational LDL decomposition gives G=L D L^t with L unit lower triangular and positive diagonal D=(d_i). The positive pivots certify positive definiteness by this congruence, and

\[
y^tGy=\sum_i d_i\left(y_i+\sum_{j>i}L_{ji}y_j\right)^2.
\]

Choose coordinates from n down to 1. After fixing larger indices, let R be the remaining budget and h_i=sum_{j>i} L_{ji}y_j. Every possible continuation satisfies

\[
|y_i+h_i|\le\sqrt{R/d_i},\qquad y_i\equiv a_i\pmod2.
\]

The verifier takes the ceiling of the square root using integer arithmetic, enumerates every integer of that parity in the resulting rational-endpoint interval, and filters by the exact inequality d_i(y_i+h_i)^2<=R. The larger interval cannot exclude a valid value. Conversely, the recursive filter admits a leaf exactly when its full norm is at most the budget. Each vector occurs once. This proves completeness through the prescribed bound, without a heuristic coordinate box or floating-point comparison. If the node limit interrupts enumeration, all partial coefficients are discarded and the verdict is unresolved.

## 5. D4 example and a uniform Dn obstruction

For the D4 example, take basis columns e1-e2, e2-e3, e3-e4, e3+e4. Then

\[
G=\begin{pmatrix}2&-1&0&0\\-1&2&-1&-1\\0&-1&2&0\\0&-1&0&2\end{pmatrix},
\quad a=(2,2,1,1)^t,\quad b=(0,-1,0,2)^t.
\]

These correspond to xi=e1 and delta=(0,0,1/2,1/2), with p=0. The double transposition (13)(24) preserves D4, shifts xi by e3-e1 in D4, and shifts delta by (1/2,1/2,-1/2,-1/2) in D4*. Pairing the latter with 2xi gives 1. The sign theorem proves vanishing. `data/d4.json` records the same map in the displayed lattice basis.

There is a uniform failure of the parity converse for Dn with n>=4. Take

\[
\xi=(1/2,1/2,1/2,1/2,0,\ldots,0),\quad
\delta=(1/2,1/2,0,0,0,\ldots,0).
\]

Here 2xi is in Dn, 2delta is an integral vector and hence in Dn*, and p=2. Flip coordinates 1 and 3. This preserves Dn, shifts xi by (-1,0,-1,0,...) in Dn, and shifts delta by (-1,0,0,0,...) in Dn*. Its sign pairing is -1. Thus every Dn with n>=4 fails the parity converse. This proves existence, not the formula counting all vanishing classes.

## 6. The parity converse holds for An

Realize An as the integer vectors of sum zero in R^(n+1). Write 2xi=u in An and 2delta=P(w), where w is integral and P projects onto the sum-zero space. This describes every element of An*: pairing with ei-ej says that all coordinate differences are integers; subtract a common real coordinate, then project to obtain w.

The residues of u modulo 2 give an even-size set S, say |S|=2k. Conversely every such S occurs. Two u with the same residues differ by 2An. For dual classes, P(w) and P(w') differ by 2An* exactly when w-w' modulo 2 is constant: if their difference is 2P(v), w-w'-2v is an integral constant vector; the reverse implication is immediate. Thus w modulo 2, up to simultaneous complementation, represents every dual class once. There are 2^n choices for each characteristic class.

In xi+An, coordinates on S are half integers and the others integers, with sum zero. The least possible norm is k/2, attained exactly by vectors with +/-1/2 on S, k of each sign, and zero elsewhere. They belong to this coset because subtracting u/2 gives an integer vector of sum zero.

Let s count the odd w_i on S. Since u_i is odd precisely there,
p=<u,P(w)>=<u,w> is congruent to s modulo 2. For a minimal vector v, choose the set U of its k negative coordinates in S. Its phase is

\[
e^{\pi i\langle v,w\rangle}
=i^{\sum_{i\in S}w_i}(-1)^{\sum_{i\in U}w_i}.
\]

Summing gives a nonzero scalar times

\[
K_k(s;2k)=[t^k](1-t)^s(1+t)^{2k-s}.
\]

For completeness, extend this coefficient to a polynomial in s using generalized binomial coefficients. Its degree is at most k. Its leading coefficient is (-2)^k/k!: in the convolution sum, the terms of degree k sum to (-1)^k sum_{j=0}^k 1/(j!(k-j)!). For odd integral s in {1,3,...,2k-1}, the generating polynomial is anti-palindromic, since t^(2k)F(1/t)=(-1)^s F(t), and hence its middle coefficient is zero. Its k distinct roots and leading coefficient therefore give

\[
K_k(s;2k)=\frac{(-2)^k}{k!}\prod_{j=1}^k(s-2j+1).
\]

For even p, s is even, so the minimal coefficient is nonzero. If k=0 the unique minimal vector is zero and its coefficient is 1; the same conclusion holds. This proves the parity converse for all n>=1. The polynomial identity is included for completeness and is not new: Chihara–Stanton's generating function (2.2) and sum (2.3) use this normalization, and the remark after Theorem 3.4 states its odd roots (author PDF pp. 2 and 5). See [docs/literature-and-source-scope.md](literature-and-source-scope.md) for the checked citation. The An shell interpretation is a separate step; no priority claim is made for it.

## 7. Orthogonal sums

Every positive-rank lattice admits an odd-p characteristic: choose a basis vector e1 and its dual e1*, and set xi=e1/2, delta=e1*/2. Their pairing p is 1. No unimodularity or even-lattice hypothesis is needed.

If L=L1 perpendicular L2 with both summands of positive rank, absolute convergence gives Theta_L=Theta_L1 Theta_L2 for the corresponding component characteristics. Choose odd characteristics on both. The product vanishes and total p=2 is even. Additional summands may be given zero characteristics. Thus every nontrivial orthogonal decomposition obstructs the parity converse. The zero lattice is excluded from the summands in this statement.

## 8. Index-two coset reformulation

Here a,b again denote the vectors 2xi,2delta, rather than their coordinates. If b is not in 2L*, the homomorphism x -> <x,b> modulo 2 is nonzero (duality identifies L* with Hom(L,Z)); its kernel L0 has index two. Choose t in L with <t,b> odd. Then L=L0 disjoint-union (t+L0), and the phase calculation in section 4 gives

\[
\Theta_L[\xi,\delta]
=i^p\left(\theta_{\xi+L_0}-\theta_{\xi+t+L_0}\right),
\qquad\theta_C=\sum_{v\in C}e^{\pi i\tau\|v\|^2}.
\]

Changing t to another odd vector changes it by an element of L0, so the second coset is well defined. Vanishing is exactly equality of the two norm-counting series, by coefficient uniqueness proved in section 4. This statement neither assumes nor establishes an ambient isometry between the cosets.

If b is in 2L*, all signs (-1)^<x,b> are 1, and Theta=i^p theta_(xi+L). At tau=it with t>0 this is a nonzero scalar times a strictly positive real sum. Hence it cannot vanish identically. Similarly, xi in L gives a unique norm-zero vector and a nonzero constant coefficient. These degenerate cases need no search.

## 9. Scope of the computational examples

The A3 certificate uses the simple-root basis ei-e(i+1). Its norm-1 shell has signed coefficient 2. The ternary certificate uses the Gram matrix in `data/ternary.json` and proves its supplied complete-shell claims; no claim about its full automorphism group is part of the certificate.

The rank-four Hermitian example uses Re(z* H w) in basis (e1,i e1,e2,i e2), so multiplication by i has the recorded integral order-4 matrix. The order-8 and order-16 certificates copy the stored Gram and generator matrices exactly. Each verifies positive definiteness, integrality, isometry, exact witness order, both stabilizer conditions and epsilon=1. The stored witness matrices can be checked directly with `code/exact_theta.py`; their minimum possible witness orders are not claimed.

These checks prove the theta verdict for the supplied matrices, even if the original algebraic construction is not independently recovered. They do not prove the full automorphism group, minimal certifying order among all its elements, indecomposability, rootlessness, CM ideal arithmetic, or unbounded minimal certifying orders.

## 10. The parity converse holds for E6

The complete result, Gram matrix, characteristic coverage and enumeration proofs are in [docs/e6-classification.md](e6-classification.md). Two exact methods verify all 2080 even-p classes. Every nonzero a class has a minimal shell of either two or ten vectors, comprising one or five antipodal pairs. For even p, antipodal signed phases agree, so an odd number of pairs cannot cancel. The zero a class has a unique norm-zero vector. Therefore every even-p characteristic has a nonzero minimal coefficient, at squared norm at most 1.

The finite E6 shell counts are certified computational evidence, checked both by the rational LDL verifier and by a separate integer Euclidean enumeration. The characteristic coverage and antipodal-pair argument are proved in the milestone note. The 2016 odd classes vanish by section 2, so E6 vanishing is exactly odd parity.

## 11. E7 and E8: complete classifications

[docs/e7-and-e8-classification.md](e7-and-e8-classification.md) supplies complete even-class coverage with two shell enumeration algorithms. All 8256 even E7 pairs are certified: 1260 vanish by symmetry and 6996 have complete nonzero coefficients. For E8, all 32896 even pairs are certified: 9450 vanish by symmetry and 23446 have complete nonzero coefficients. The odd pairs vanish by section 2.

Zero certificates comprise an exact root witness with epsilon=1 and a checked tree of isometries transporting it to every zero class. Isometry covariance and characteristic shift laws prove vanishing at each vertex; the conjugated witness can also be exported to the original verifier. No complete-group claim or finite-determination theorem is assumed.

The E7 dual-glue restriction, the E8 totally singular-plane criterion, its 1575 planes and its single zero orbit are certified over complete class sets. The milestone gives the quadratic-form conventions and an explicit totally singular four-space certifying plus type.

**ADE parity-converse theorem.** Among irreducible ADE root lattices, the parity converse holds exactly for An and E6. The positive directions are sections 6 and 10; the negative directions are the uniform Dn obstruction in section 5 and the supplied E7/E8 zero certificates. The complete Dn counts are proved independently in section 12.

## 12. Uniform Dn criterion and counts

[docs/dn-vanishing-formula.md](dn-vanishing-formula.md) proves the necessary-and-sufficient vanishing criterion for every n>=2, with complete primal and dual characteristic parametrizations and the phase-preserving parity projection. The number of even-p zeros is

\[
Z_n=4^{n-1}-3^n+3\,2^{n-1}-1
    +\mathbf1_{2\mid n}\frac12\binom n{n/2}.
\]

All even zeros have explicit involution witnesses. All other classes have nonzero minimal coefficients, with squared norm bounded by max(1,floor(n/2)/2). Exact basis data identify D2 with A1 perpendicular A1 and D3 with A3, giving respectively one and zero even zeros.

The saved D2-D10 streams independently check all 700,070 even pairs using rational LDL enumeration and ambient integer spheres. Every zero passes its actual stabilizer and epsilon equations; every nonzero passes a complete shell check. Both enumeration histograms agree, and no finite initial cancellation is used as a proof of identity. This establishes the uniform formula symbolically and reproduces the finite counts computationally.

## 13. Nonsimply-laced types and the root-lattice theorem

[docs/root-lattice-classification.md](root-lattice-classification.md) fixes integral short-root normalization, proves positive-scaling invariance, and identifies the root lattices: Bn is an orthogonal sum of n copies of A1; Cn is Dn for n>=2; F4 is D4; G2 is A2. It supplies exact unimodular basis data and complete F4/G2 certificates. For Bn, the uniform even-zero count is 2^(2n-1)+2^(n-1)-3^n, proved by coordinate factorization and checked through rank eight.

**Root-lattice parity-converse theorem.** For lattices generated by reduced irreducible crystallographic root systems, the positive types are An, E6, C3 and G2; the optional rank-one B1/C1 conventions are A1. Up to normalized lattice isometry, the positive lattices are exactly An and E6. C3 and G2 contribute A3 and A2, rather than additional isometry classes. Root-system irreducibility does not guarantee lattice indecomposability: Bn and C2 illustrate the distinction.

The same note qualifies the E8 two-coset argument: the binary S/complement calculation applies when 2xi is integral. Quarter-coordinate xi requires its own analysis, and the number of products in an expression is not an established invariant. An exported quarter-coordinate E8 zero has an exactly verified symmetry witness.

## 14. Rank-one/two finite-domain classification

**Theorem (determinant 1 through 24).** Among positive definite integral lattices
of rank one or two with Gram determinant at most 24, all 24 rank-one classes
satisfy the parity converse. Among the 70 rank-two classes, precisely the 26
orthogonally indecomposable classes satisfy it; the 44 decomposable classes fail.
For every characteristic in this domain, theta vanishes identically if and only
if epsilon is nontrivial on the full characteristic stabilizer.

**Proof by exact finite certification.** Binary reduction gives
`G=[[A,B],[B,C]]`, `1<=A<=C`, `0<=2B<=A`, and det(G)>=3A^2/4, hence A<=5.
The finite generator includes every such form of determinant 1 through 24,
together with every rank-one matrix [d]. The coverage proof and exhaustive
basis-image isometry test are detailed in [docs/rank-one-and-two-census.md](rank-one-and-two-census.md).
Independent enumeration verifies that the 94 generated matrices are pairwise
nonisometric and checks the full automorphism groups and pair stabilizers.

The separately saved [closure pack](../data/census/rank12-det24-closed.json.gz)
covers all 772 even pairs exactly once. Every one of its 728 nonzero entries has
a complete nonzero shell coefficient with norm4<=25. Each of the other 44 entries
has an exact stabilizing isometry with epsilon one. Independent box and LDL
replay agree, with no unresolved computations. Section 4 proves that these are
nonvanishing certificates; section 2 proves that the symmetry entries vanish.
The counts split as 72 nonzero rank-one pairs and 656 nonzero/44 zero rank-two
pairs. At lattice level, all rank-one classes and 26 rank-two classes are positive.
Each of the 44 diagonal rank-two representatives has one even zero.

Every decomposable integral rank-two lattice has a basis from its two orthogonal
rank-one summands, hence a diagonal Gram matrix. Its determinant bound places that
matrix among the generated diagonal representatives. Exact pairwise nonisometry
therefore proves that none of the non-diagonal representatives is decomposable.
The diagonal classes fail also by section 7.

For even pairs the complete list proves that vanishing always has an epsilon-one
witness. For odd pairs -I supplies that witness. The converse is section 2.
This proves the stated symmetry equivalence throughout the finite domain.
It supplies no all-determinant binary classification or unrestricted symmetry
converse by itself. Section 15 removes the determinant restriction in ranks one
and two. The original norm4=16 run and its 33 unresolved records are preserved
as a separate initial experiment.

## 15. Uniform rank-two classification and low-rank symmetry converse

[docs/rank-two-classification.md](rank-two-classification.md) gives the complete proof, reviewed from
reviews/initial-mathematical-review.md section 2, with no determinant restriction.

**Theorem.** A positive definite integral rank-two lattice has the parity
property if and only if it is orthogonally indecomposable. A decomposable
rank-two lattice has exactly one even zero. The symmetry converse holds for
every characteristic in ranks one and two, at every determinant.

In a reduced Gram matrix `[[A,B],[B,C]]`, `1<=A<=C`, `0<=2B<=A`, decomposition
is equivalent to B=0. Seven of the ten even pairs are settled by a zero vector
or constant signs. The remaining pairs have complete minimal coefficients 2
at norm4=A, 2 at norm4=C, and -2 at norm4=A+C-2B when B>0. The standalone
proof gives inequalities excluding every other vector, including the reduction
boundaries. When B=0, only a=b=(1,1) vanishes, with witness diag(-1,1).
Odd pairs have witness -I, and rank-one even pairs are all nonzero.

This theorem strengthens the classification and symmetry conclusion of section
14 to all determinants in these ranks. The finite census and its original
bounded artifacts retain their stated domains and verdicts. The rank-three companion completes the next rank; the rank-four companion disproves the unrestricted converse.
