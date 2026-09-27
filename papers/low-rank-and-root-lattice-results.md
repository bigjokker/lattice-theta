# Vanishing lattice theta functions with half-integral characteristics

Completed lower-rank and root-lattice results. The rank-four companion supplies the subsequent counterexample. These are research manuscripts; no publication or priority claim is made.

## Abstract

For a positive definite integral lattice, we study identical vanishing of its
one-variable theta functions with half-integral primal and dual characteristics.
The characteristic stabilizer carries a sign character whose nontriviality
forces vanishing and includes classical odd-characteristic vanishing. In ranks
two and three, every even characteristic is nonzero as a function precisely when
the lattice is orthogonally indecomposable, at every determinant. Every
identically zero characteristic in ranks one through three has a stabilizer
element with sign minus one. The ternary proof uses an obtuse superbase and
complete signed-shell evaluations, including zero conorms and cancellation on
the diamond graph. Exact censuses cover all integral lattice classes of ranks
one through three and determinant at most 24. We also give an inclusive ternary
coefficient cutoff from coset modularity. For irreducible crystallographic
root lattices in the stated integral normalization, the positive types are
An, E6, C3 and G2, the last two being A3 and A2. The general classification of the parity property in rank at least four remains open. The symmetry converse fails in rank four, as proved in the companion manuscript [rank-four-counterexample.md](rank-four-counterexample.md).

## 1. The two converse questions

Let L be a positive definite integral lattice of positive rank in its real span,
with dual L*. For 2xi in L and 2delta in L* define

\[
\Theta_L[\xi,\delta](\tau)=\sum_{x\in L}
e^{\pi i\tau\|x+\xi\|^2+2\pi i\langle x+\xi,\delta\rangle},
\qquad \operatorname{Im}\tau>0,
\quad p=\langle2\xi,2\delta\rangle\in\mathbb Z.
\]

Absolute convergence permits the reindexings below. Characteristic classes are
xi modulo L and delta modulo L*. Shifting xi leaves theta unchanged; shifting
delta by nu in L* multiplies it by exp(2 pi i <xi,nu>). Thus identical vanishing
and parity of p are properties of the pair of classes.

Call **P(L)** the property that every even-p characteristic has nonzero theta.
Our first question is to classify lattices with P(L). The separate **symmetry
converse** asks whether every identically vanishing theta with even p has an
epsilon-one element in its full characteristic stabilizer, defined next. A
lattice can fail P(L) while every zero is symmetry-explained.

In a lattice basis and its dual basis, write a,b for the integer coordinates
of 2xi,2delta and G for the Gram matrix. This series is the Riemann theta
constant at Omega=tau G with characteristic (a/2,b/2). Positive definiteness
ensures that Omega belongs to Siegel space. Consequently our question concerns
identical vanishing along this particular one-variable family of period
matrices. It is distinct from the zero of an even theta constant at one period
matrix. The normalization is [DLMF (21.2.5)](https://dlmf.nist.gov/21.2.E5),
and general half-characteristic parity is [DLMF (21.3.6)](https://dlmf.nist.gov/21.3.E6).
At z=0 its sign is (-1)^p, giving odd-p vanishing for every period matrix,
hence along Omega=tau G. Our direct proof uses the primal/dual hypotheses above.

## 2. Stabilizer sign and complete coefficients

Let H be the automorphism stabilizer of the two characteristic classes:

\[
H=\{\sigma\in\operatorname{Aut}(L):\sigma\xi-\xi\in L,
\ \sigma\delta-\delta\in L^*\},\qquad
\varepsilon(\sigma)=\langle2\xi,\sigma\delta-\delta\rangle\pmod2.
\]

**Proposition 1.** Epsilon is independent of representatives and is a
homomorphism H to F2. If epsilon is nontrivial, theta vanishes identically.

**Proof.** H is a pair stabilizer, hence a subgroup. Put mu_sigma=sigma delta-delta
in L*. Shifting xi by ell changes the pairing by <2ell,mu_sigma>, an even
integer. Shifting delta by nu changes it by
<sigma^{-1}(2xi)-2xi,nu>, also even because sigma stabilizes xi. Moreover
mu_(sigma rho)=sigma mu_rho+mu_sigma. Pairing sigma mu_rho with 2xi agrees modulo
2 with pairing mu_rho with 2xi, by the same stabilizer condition. This proves
the homomorphism property. Orthogonal reindexing and the characteristic shift
laws give

\[
\Theta_L[\xi,\delta]=\Theta_L[\sigma\xi,\sigma\delta]
=(-1)^{\varepsilon(\sigma)}\Theta_L[\xi,\delta].
\]

An element with epsilon one therefore forces identical vanishing.

The element -I belongs to H and has epsilon=-p modulo 2. This proves odd-p
vanishing. Epsilon is zero on odd-order elements and squares in H, but these
facts do not bound witness orders or prove the symmetry converse.

For a primal-column matrix T the dual action is T^{-t}, so the coordinate test is

\[
T^tGT=G,\quad Ta\equiv a\ (2),\quad T^{-t}b\equiv b\ (2),
\qquad \varepsilon(T)=a^t\frac{(T^{-t}-I)b}{2}\pmod2.
\]

The division is exact after the dual stabilizer test. Both congruences are
necessary; ambient group generators that individually fix a pair need not
generate its whole stabilizer.

This sign is also obtained from the standard integral basis-change law
[DLMF (21.5.5)](https://dlmf.nist.gov/21.5.E5) and characteristic shift law
[DLMF (21.3.4)](https://dlmf.nist.gov/21.3.E4). For U in GL(n,Z), reindexing the
defining series gives the exact unreduced-characteristic identity

\[
\theta\!\begin{bmatrix}U^{-1}\alpha\\U^t\beta\end{bmatrix}
(U^t z\mid U^t\Omega U)
=\theta\!\begin{bmatrix}\alpha\\\beta\end{bmatrix}(z\mid\Omega).
\]

There is no extra prefactor in this identity, including det U=-1. If U=T^{-1}
and T preserves G, the period matrix is unchanged. Reducing the transformed
characteristic (T alpha,T^{-t} beta) contributes precisely (-1)^epsilon(T).
The multiplier interpretation is standard; the self-contained homomorphism
and vanishing proof above makes its lattice conventions explicit. See
docs/theta-transformations.md for the convention comparison with the symplectic block.

Put y=2x+a. After removing the global phase i^p, the coefficient at norm4=N is

\[
c_N=\sum_{\substack{y\equiv a\ (2)\\y^tGy=N}}
(-1)^{((y-a)/2)^tb},\qquad
\Theta=i^p\sum_N c_N e^{\pi i\tau N/4}.
\]

A complete nonzero coefficient proves nonvanishing. To see this, choose the
least N with c_N nonzero, set tau=it, multiply by exp(pi t N/4), and let t tend
to infinity; absolute convergence controls the tail and the limit is i^p c_N.
Finite initial cancellation alone does not prove identical vanishing.

The supplement enumerates bounded spheres by exact rational LDL recursion.
Each recursive coordinate interval contains every continuation permitted by
the remaining norm budget. An independent backend uses inverse-Gram bounds
|y_i|^2 <= N(G^{-1})_ii to enumerate complete integer boxes. A resource limit
discards partial coefficients and returns unresolved. See docs/foundations-and-exact-enumeration.md section 4
for the full enumeration proof.

## 3. Orthogonal sums and the uniform binary theorem

**Proposition 2.** Every nontrivial orthogonal sum of positive-rank integral
lattices fails P(L).

**Proof.** Every summand has an odd characteristic: take half of a basis vector
and half of its dual vector, giving p=1. Choose these on two summands and zero
characteristics on any others. The theta series factors by absolute convergence,
both chosen factors vanish, and the total p=2 is even.

**Theorem 3.** A positive definite integral rank-two lattice satisfies P(L) if
and only if it is orthogonally indecomposable. In a decomposable rank-two lattice
exactly one even class vanishes. The symmetry converse holds in ranks one and
two at every determinant.

**Proof.** Choose a reduced basis
`G=[[A,B],[B,C]]`, `1<=A<=C`, `0<=2B<=A`. Such a basis exists by integer column
subtraction and swaps, with each swap decreasing the positive integer A.
For N(m,n)=Am^2+2Bmn+Cn^2, opposite signs give, with u=|m| and v=|n|,
N=Au^2-2Buv+Cv^2>=Au^2-Auv+Cv^2=Au(u-v)+Cv^2, using 2B<=A.
When u>=v this is at least A; when u<v, C>=A gives
N>=A((u-v/2)^2+3v^2/4)>=3A. Same signs or a zero coordinate are immediate.
Hence A is the shortest nonzero norm.

If L=Ze perpendicular Zf, with norms d<=h, its shortest vector is on one of
the two axes. The first reduced vector must be such an axis vector; the second
is xe+sf, s=+/-1, after relabeling. The reduced bound |2B|<=A=d forces x=0.
This also covers d=h. Therefore decomposition is equivalent to B=0 in a reduced
basis.

There are ten even binary pairs a,b. Four with a=0 have constant coefficient
1, and three further pairs with b=0 have only positive signs. The other three
are resolved by complete minimal shells:

| a | b | norm4 | Minimal vectors y | Coefficient |
|---|---|---|---|---:|
| (1,0) | (0,1) | A | +/-(1,0) | 2 |
| (0,1) | (1,0) | C | +/-(0,1) | 2 |
| (1,1), B>0 | (1,1) | A+C-2B | (1,-1), (-1,1) | -2 |

For the first coset any nonzero even second coordinate has absolute value
v>=2. The opposite-sign bound is at least 4C if u>=v, or 3A otherwise;
same signs give at least 4C. Thus only +/-(1,0) attains A.
For the second coset, if the even first coordinate is nonzero then u>=2.
Same signs give N>C. For opposite signs,
N-C>=Au(u-v)+C(v^2-1), which is positive if u>=v since their parities differ.
If u<v, then v>=3 and the bound is at least A(3v^2/4-1)>0.
Both axis shells have sign +1 at each vector.

For the odd/odd coset put M=A+C-2B. Same signs give N>=A+C>M when B>0.
For opposite signs u,v>=1, so uv-1>=0, and 2B<=A gives
N-M=A(u^2-1)-2B(uv-1)+C(v^2-1)>=Au(u-v)+C(v^2-1).
This lower bound can equal zero only at u=v=1; for u<v, v>=3
again makes it strictly positive. The two displayed vectors give x=(0,-1),
(-1,0), both with sign -1. These inequalities include A=C and 2B=A.

Thus every even pair is nonzero if B>0. If B=0, only a=b=(1,1) remains.
The matrix T=diag(-1,1) preserves the diagonal G, and
(T-I)a=(T^{-t}-I)b=(-2,0). Its dual shift is (-1,0), so epsilon=-1 modulo 2
is one. The other nine are nonzero,
so there is exactly one even zero. For rank one, a=0 gives constant coefficient
1 and the only other even pair a=1,b=0 has minimal coefficient 2. Odd pairs
in both ranks have witness -I. Proposition 1 proves that an epsilon-one witness
forces vanishing. Conversely, the witnesses -I and diag(-1,1) cover every zero
in these ranks, proving the low-rank symmetry converse. This proves all
assertions. The expanded proof is RANK-TWO.md.

## 4. Uniform ternary classification and symmetry converse

**Theorem 4.** A positive definite integral rank-three lattice satisfies P(L)
if and only if it is orthogonally indecomposable. Every identically zero
half-integral characteristic in rank three has an epsilon-one stabilizer
element. Both statements hold at every determinant.

The obtuse-superbase framework is classical. See
[Conway-Sloane, section 2.0.5, equation (5), and Theorem 8](https://neilsloane.com/doc/fedorov.pdf),
and [Kurlin, Theorem 2.8 and Appendix A](https://kurlin.org/projects/lattice-geometry/lattices3Dmaths.pdf).
Our integral descent and signed-shell argument are given in full below;
neither a finite census nor an isometry-classification table is a premise.

**Proof: decomposable lattices.** A rank-three orthogonal splitting consists
of three lines or a line and an indecomposable binary lattice. An even zero
on three lines has exactly two odd factors: the even factors are nonzero by
the rank-one argument in Theorem 3. On a line plus an indecomposable binary
lattice, Theorem 3 similarly makes every even zero a product of two odd
factors. Products of nonzero factor series are nonzero because their least
nonzero coefficients multiply at the sum of their indices. Negate either
odd factor and extend by identity on the other factors; this gives an
epsilon-one isometry. Odd characteristics have the witness -I in all cases.
Proposition 2 gives failure of P for every nontrivial orthogonal sum.
It remains to prove nonidentity for every even characteristic of an
indecomposable rank-three lattice.

### 4.1. Obtuse superbases by a terminating integer reduction

Start with a lattice basis v1,v2,v3 and v0=-v1-v2-v3. Any three of these
four vectors form a unimodular basis. Suppose vi dot vj>0, with k,l the
other indices. Replace the superbase by

    vi'=-vi, vj'=vj, vk'=vk+vi, vl'=vl+vi.

Its sum is zero. Omitting vj, its three columns are a unimodular change
of the previous basis omitting vj. The sum of all four squared norms
decreases by 2 vi dot vj: use vk+vl=-vi-vj. This is a decrease of at
least two in a positive integer. Hence repeating the step terminates
with vi dot vj<=0 for all distinct i,j. This is an obtuse superbase;
the argument uses no determinant bound or canonical reduction claim.

Set w_ij=-vi dot vj>=0 for i!=j. In the coordinates
x=(x0,x1,x2,x3), modulo addition of a common integer to every coordinate,

    y=sum_i xi vi,
    ||y||^2=sum_(i<j) w_ij (xi-xj)^2.

The identity follows by expanding, since ||vi||^2=sum_(j!=i) w_ij.
Every lattice vector occurs exactly once after fixing any xr=0.

### 4.2. The positive-edge graph

Join i and j when w_ij>0. A nonempty proper subset of a superbase cannot
sum to zero: after subtracting the coefficient of a chosen root, the
other three coefficients would give a nontrivial relation among a basis.
Equivalently, the only relation among the four vectors is their total sum.
If the graph were disconnected, take x to be the indicator of one
component. Its corresponding nonzero proper subset sum would have norm
zero, a contradiction. Thus the graph is connected.

If it has an articulation vertex r, omit vr. The complementary triple
is a unimodular basis, and components of the graph with r deleted have
no positive edges between them. Since conorms are nonnegative, all such
cross inner products are zero. After grouping that basis by components,
its Gram matrix is block diagonal with nonempty blocks. Their basis spans
are saturated orthogonal nontrivial summands. Therefore an indecomposable
lattice's positive-edge graph has no articulation vertex.

A connected simple graph on four vertices without an articulation
vertex has minimum degree two. Its possibilities are the four-cycle,
the complete graph missing one edge (diamond), or K4. For four edges,
degree two everywhere forces a cycle; for five edges it is the diamond;
for six it is K4. Three or fewer edges cannot meet the degree bound.

### 4.3. Characteristics and a complete lower bound

A primal class is a parity vector a on the four vertices modulo adding
the all-one vector. After choosing root r with ar=0, fix xr=0 and require
xi congruent ai mod 2. A dual functional is an integer vector beta with
sum beta_i=0, where beta_i=<vi,2delta>. Its sign is

    (-1)^(sum_i beta_i (xi-ai)/2).

If the vector sum ai vi differs from the original primal representative
by 2ell, the displayed sign differs from the original coefficient sign by
the constant (-1)^<ell,2delta>. Changing the chosen lift of a therefore only
multiplies all coefficients by the same nonzero sign. Root changes and
adding a common integer do not alter y because sum vi=0 and sum beta_i=0.
For a new root r, subtract ar from all parity coordinates and xr from
all integer coordinates, then normalize parity representatives by even
shifts. This preserves the class and changes the sign only by the constant
just described. The parity is p=sum ai beta_i mod 2. Thus a can be
complemented or the vertices relabelled without changing whether a
complete coefficient is zero. For the original signs s(y),

    s(-y)/s(y)=(-1)^(-y^t b)=(-1)^(a^t b)=(-1)^p,

since y congruent a mod 2. Thus antipodal signs agree exactly when p is
even and are opposite when p is odd.

Let S be the odd vertices and T its complement. For every y in the class,

    ||y||^2 >= m(S)=sum_(i in S,j in T) w_ij.

Across the cut the difference is odd, so its square is at least one;
within a parity part it is even, so its square is at least zero.
Equality forces every positive within-part edge to have difference zero
and every positive cross edge to have difference +/-1. These conditions
are necessary and sufficient for equality, and completely enumerate the
minimal shell in the cases below. In these coordinates y^t G y is the
original norm4 index. When S is empty the norm-zero shell contains only
y=0. Its coefficient is +1 for a=0, and +/-1 after an even change of
representative. Hence every zero primal class is nonzero for every b.

### 4.4. All indecomposable cases

Complement S if needed so |S|<=2.

**One odd vertex.** Its complement is connected because the full graph
has no articulation vertex. Fix a root in the complement. Equality
forces its three coordinates to be zero and the odd coordinate to be
+/-1. These are exactly two vectors. Even p makes their signs agree,
so the coefficient is +/-2.

**Two odd vertices, both within-part edges positive.** Equality makes
each part constant. With a root in T, the two coordinates in T are
zero and both coordinates in S are +/-1. Again the complete shell has
one antipodal pair and coefficient +/-2.

**Both within-part edges zero.** The positive-edge graph must be the
four-cycle with all four cross edges positive. Put T={r,t}, xr=0.
Equality allows xt=0, with the two odd coordinates independently
+/-1 (four vectors), or xt=+/-2 with both odd coordinates respectively
+/-1 (two more vectors). No other xt is possible because two length-one
cross differences connect r to t. There are exactly three antipodal
pairs. Their even-p signed contributions are each +/-2, whose sum
cannot be zero.

**Exactly one within-part edge zero.** Complement S so its within edge
is zero; write T={r,t} and h=w_rt>0. The graph is the diamond and all
four cross edges are positive: a missing cross edge would leave one
vertex of S with degree at most one. At the minimum m=m(S), equality forces
xr=xt=0 and the two odd coordinates independently +/-1, giving four
vectors. Even p means the beta values on S have the same parity. If
both are even, all four signs agree and the coefficient is +/-4.

If both are odd, the minimal coefficient cancels. In fact the whole
slice xt=0 cancels at every norm: negating either odd coordinate keeps
its norm, since its only edges meet zero coordinates, and reverses its
sign. This involution has no fixed point on the odd coordinate.

Write xt=2k. Every vector with |k|>=2 has norm at least m+16h, by the
cross-edge lower bound and the contribution h(2k)^2. For k=+/-1 the
least norm is m+4h, attained uniquely by setting both odd coordinates
to +/-1, respectively: every cross difference must have magnitude one.
These two vectors are antipodal and have equal signs for even p.
The k=0 contribution cancels and the |k|>=2 slices lie strictly above
this norm. Consequently the COMPLETE coefficient at m+4h is +/-2.
These two contributing vectors are +/- (vt-vr). The shell can also
contain vectors from k=0; those cancel within that slice. Therefore the
two contributing vectors need not exhaust the shell. This handles the
cancellation boundary without requiring any genericity or unequal-conorm
hypothesis. Independent replay of this example sums the entire shell.

For example, G=[[2,0,-1],[0,2,-1],[-1,-1,6]] has determinant 20 and
occurs in the ternary census. For a=(1,1,0), its norm4=4
shell has coefficient zero when b=(1,1,0) or (1,1,1). At norm4=20,
the full shell has ten vectors; its signed coefficients are respectively
+2 and -2. Eight of those vectors cancel on the k=0 slice, leaving the
two contributing vectors. Here h=4 and m=4.

### 4.5. Conclusion

The four cases prove that every even characteristic of an indecomposable
rank-three integral lattice has a complete nonzero coefficient. Orthogonal
sums fail P by the established product obstruction. Therefore

    P(L) holds in rank three iff L is orthogonally indecomposable.

Combining this with the decomposable-factor argument above proves the
symmetry converse for all even characteristics in rank three, at every
determinant. Every odd characteristic has the witness -I, with
epsilon(-I)=-a^t b=p mod 2, by Proposition 1. Thus the converse
holds for ALL half-integral characteristics in rank three. The graph
reasoning also implies that a positive-edge graph
without articulation cannot represent a decomposable lattice; otherwise
its proved P would contradict the product obstruction. This last
implication is a consequence, not an assumed graph classification.

There is no extrapolation from the 120-class census in this argument.
The shell index m or m+4h is computed in an obtuse superbase and is the
original norm4 index. Degenerate conorms, equal norms, arbitrary positive
integer scalings and boundary symmetries are included explicitly.

## 5. An and the root-lattice classification

Realize An as integer vectors of sum zero in R^(n+1). Write 2xi=u in An and
2delta=P(w), w integral, where P is orthogonal projection onto the sum-zero
space. This represents every dual vector, since integral pairings with ei-ej
are exactly the integral coordinate-difference condition.

Let S be the set of odd coordinates of u; its size is 2k. These residues
represent all primal classes. The binary residues of w, up to simultaneous
complementation, represent all dual classes. The complete minimal vectors in
xi+An have +/-1/2 on S, k of each sign, and zero elsewhere; their norm is k/2.
Let s be the number of odd w coordinates in S. Then p=<u,w> has parity s,
and the phase sum over this shell is a nonzero scalar times

\[
K_k(s;2k)=[t^k](1-t)^s(1+t)^{2k-s}
=\frac{(-2)^k}{k!}\prod_{j=1}^k(s-2j+1).
\]

In the normalization of [Chihara–Stanton](https://www-users.cse.umn.edu/~stant001/PAPERS/krlaur.pdf),
this is k_k(s,2,2k): their (2.2) identifies the polynomial, (2.3) determines the
leading coefficient, and the remark after Theorem 3.4 states the odd roots.
The lattice-shell evaluation above is a separate argument. For even p, s is
even, so the coefficient is nonzero. The k=0 case has its unique norm-zero
vector. This proves P(An) for every n>=1; the supplement includes a direct
proof of the polynomial factorization as well.

For Dn, n>=4, use xi=(1/2,1/2,1/2,1/2,0,...) and
delta=(1/2,1/2,0,...). Here p=2, and flipping coordinates 1 and 3 preserves
Dn and both characteristic classes, with epsilon=-1 modulo 2. Therefore P(Dn)
fails. D2 decomposes and D3 is A3. docs/dn-vanishing-formula.md proves a uniform criterion
and counts all even zeros; this additional counting theorem is part of the
supplement rather than a premise of the obstruction.

The exceptional cases are resolved by exact finite certification:

| Lattice | Even pairs | Symmetry zeros | Nonzero shells | Bound on squared norm of certified nonzero shells |
|---|---:|---:|---:|---:|
| E6 | 2080 | 0 | 2080 | 1 |
| E7 | 8256 | 1260 | 6996 | 3/2 |
| E8 | 32896 | 9450 | 23446 | 1 |

The packs contain the Gram matrices, each binary pair once, and explicit shell
or symmetry evidence. E7/E8 transport trees certify all zero classes from checked
root witnesses; every transport step is an integral isometry with the correct
primal/dual action. Rational LDL and independent ambient Euclidean sphere
enumerations agree on all histograms. Coverage has the expected
2^(n-1)(2^n+1) even pairs. Thus no initial-zero search is used as an identity
proof. docs/e6-classification.md and docs/e7-and-e8-classification.md give the realizations, independent
enumeration proofs and complete replay instructions. These finite checks prove
P(E6) and disprove P(E7), P(E8).

**Theorem 5.** For lattices generated by reduced irreducible crystallographic
root systems, with short roots of squared length 2, the positive types are
An,E6,C3,G2, including the optional rank-one B1/C1 conventions as A1.

**Proof.** The ADE cases follow above. In this integral normalization Bn is
an orthogonal sum of n rank-one lattices, Cn is Dn for n>=2, F4 is isometric
to D4, and G2 to A2. docs/root-lattice-classification.md proves these identifications from the
generators and provides exact unimodular basis changes. Proposition 2 deals
with Bn and C2; C3=D3=A3 is positive and Cn for n>=4 is negative. The F4/G2
packs independently check their complete 136/10 even pairs. This finishes the
case list. Positive rescaling preserves vanishing via xi'=xi, delta'=delta/c
when the metric is multiplied by c, so it does not change these parity verdicts.

## 6. Census as an independent finite check

For rank one enumerate G=[d], 1<=d<=24. In rank two the reduced bounds imply
det G>=3A^2/4, hence A<=5. Enumerating all A,B and
A<=C<=floor((24+B^2)/A), and retaining positive determinants at most 24,
covers the domain. An exhaustive isometry test enumerates all candidate
images of the source basis in the target sphere of squared norm max(G_ii),
then checks the full Gram equation and determinant +/-1. Taking source equal
to target gives the full automorphism group.

There are 94 pairwise nonisometric generated forms: 24 rank one and 70 rank two.
The closed pack resolves all 772 even pairs as 728 nonzero and 44 symmetry
zeros. It gives 50 positive lattice classes and 44 negative classes. Every
minimum, full group, stabilizer sign and shell claim is replayed independently.
The results agree with Theorem 3, including exactly one even zero per diagonal
rank-two class. The initial norm4=16 run is retained with its 33 unresolved
records; a separate norm4=25 run supplies their complete nonzero coefficients.
The follow-up changes neither the lattice domain nor the historical evidence.

The rank-one/two census is not needed to prove Theorem 3. It checks the
implementation and certificates independently of that uniform argument.

### 6.1. Complete ternary domain generation

For determinant Delta<=H=24 choose a shortest primitive vector v1 of norm A.
It extends to an integral basis. Projecting along v1 gives a discrete rank-two
lattice. Gauss-reduce a projected basis and lift it to v2,v3, subtracting
integer multiples of v1 to reduce both pairings. Then

    G=[[A,B,D],[B,C,F],[D,F,E]],
    |2B|,|2D|<=A,
    S=C-B^2/A, Q=F-BD/A, R=E-D^2/A,
    |2Q|<=S<=R.

The projected reduction terminates since its squared norms lie in (1/A)Z
and every swap decreases a positive first norm. Shortestness gives C>=A,
so S>=3A/4 and

    Delta=A(SR-Q^2)>=3AS^2/4>=27A^3/64.

Thus A<=3. Put M=AC-B^2. Enumerate positive A and C>=A, integral B,D with
the reduced pairing bounds, and 3M^2<=4AH. Enumerate integer F with
|2(AF-BD)|<=M. For each Delta=1,...,H solve

    E=(Delta+CD^2+AF^2-2BDF)/M.

Retain integral E with AE-D^2>=M. These finite inequalities cover a basis
for every lattice in the declared domain; they need not give a unique basis.
Positive leading minors prove definiteness. For an independent scanner,
C,E<=2H+A suffice: S<=sqrt(4H/(3A))<=2H, C<=2H+A/4, and
E=Delta/(AS)+Q^2/S+D^2/A<=11H/6+A/4.

Exact isometry testing enumerates all target images of the three source basis
vectors at their prescribed norms, checks every pairing and determinant
plus or minus one. Each possible isometry occurs in that list. Applying the
same procedure to a form itself gives its full automorphism group. The two
generation implementations use determinant-solving and bounded-diagonal
loops, respectively. Both sphere backends use proved complete bounds.

### 6.2. Ternary results and decomposition certificates

The 303 generated matrices reduce to 120 isometry classes, including odd,
even, nonprimitive and decomposable lattices, with scales kept distinct.

| Quantity | Certified value |
|---|---:|
| Integral isometry classes | 120 |
| Even characteristic pairs | 4,320 |
| Even symmetry zeros | 723 |
| Even complete-coefficient nonzeros | 3,597 |
| Odd zeros with witness -I | 3,360 |
| Unresolved census pairs | 0 |
| Indecomposable lattices satisfying P | 25 |
| Decomposable lattices failing P | 95 |

Every minimum, full automorphism group, characteristic stabilizer, sign and
coefficient is independently replayed. Norm4=32 is a search resource choice,
not an identity cutoff. The largest nonzero certificate index is 26; 22
nonzeros have a cancelled first shell. No interrupted search is certified
complete. These checks corroborate Theorem 4 and do not prove its unrestricted
statement.

To certify splitting independently, a rank-one summand Zv with primitive v
and d=v^t G v exists precisely when each entry of Gv is divisible by d.
The integral row lambda=(Gv)^t/d satisfies lambda v=1 and gives the saturated
orthogonal decomposition x=(lambda x)v+(x-(lambda x)v). Conversely a summand
has this integral projection. Every possible generator has d<=max G_ii:
some supplied basis vector has a nonzero integral v-component, whose squared
norm is at least d. A complete sphere through this bound therefore excludes
all splittings when the test fails. Bezout completion and projection of the
remaining basis columns provide an explicit unimodular splitting basis.
Binary reduction of its complement distinguishes the two decomposable types.

Independent replay yields 51 sums of three lines, 44 line/binary sums and
25 indecomposables. A three-line sum has 27 even nonzeros and nine even zeros;
a line plus an indecomposable binary lattice has 30 even nonzeros and six
even zeros. All 7,680 structural characteristic records replay, with 4,083
zeros and 3,597 nonzeros. The 25 indecomposable superbase graphs comprise
12 cycles, 11 diamonds and two complete graphs. All 900 even characteristics
on them have checked symbolic-shell predictions.

## 7. Coset reformulation and modularity

Put a=2xi and b=2delta as actual vectors. If b is not in 2L*, let
L0={x in L:<x,b> even}, an index-two sublattice, and choose t in L with
<t,b> odd. Splitting the sum gives

\[
\Theta_L[\xi,\delta]=i^p
\bigl(\theta_{\xi+L_0}-\theta_{\xi+t+L_0}\bigr).
\]

If b is in 2L*, signs are constant and the series is nonzero at tau=it, t>0.
Equality of the two norm-counting coset series is the identity needed for
vanishing in the other case. It should not be presumed to arise from an
automorphism of the ambient lattice. Work on unsigned ternary lattice-coset
series, such as [Kane–Kim](https://arxiv.org/abs/2211.03987), provides relevant
background but does not directly settle this stabilizer question.
More precisely, for even p, doubling the two cosets gives integral cosets
a+2L0 and a+2t+2L0. In the unsigned variable z, Proposition 2.3 of that
paper places both summands and their difference in weight n/2 on
Gamma0(16N_L0), where N_L0 is the least denominator-clearing integer for
the inverse Gram matrix of L0. The source character is chi_(4 det L0)
for odd n and chi_((-1)^(n/2) 4 det L0) for even n. A trivial coset uses
base 2L0 and conductor one; a nontrivial one uses base L0 and conductor
two. The series substitution is z=tau/8.

**Proposition 6 (rescaling).** Assume p is even and b is not in 2L*. Let chi be the source character just defined.
The signed theta is a holomorphic form of weight n/2 on

\[
H=\Gamma_0(16N_{L_0})\cap\Gamma^0(8),
\]

where Gamma^0(8) requires the upper-right entry to be divisible by eight.
For even n its character is chi; for odd n, in the source's half-integral
slash convention, its character is chi chi_8, with chi_8(v)=(8/v).
The group has cusp width eight at infinity. This uses the modular-form
definition permitting arbitrary cusp widths; the source's convention
requiring translation by one does not apply to H as stated.

**Proof.** Put F=Phi_(a+2L0)-Phi_(a+2t+2L0), so Theta=i^p F(tau/8).
For gamma=[[r,s],[u,v]] in H, gamma'=[[r,s/8],[8u,v]] lies in
Gamma0(16N_L0) and gamma'(tau/8)=gamma(tau)/8. The analytic denominator
is u tau+v in both transformations. Integral weight leaves chi(v) unchanged;
in odd rank the symbol changes by (8u/v)/(u/v)=(8/v). Use the principal
argument in (-pi,pi]; the common denominator makes the branch comparison
valid also for negative v. At -I the half-integral slash factor is
i^n exp(-pi i n/2)=1, matching chi(-1)chi_8(-1)=1; in even rank the
two factors (-1)^(n/2) cancel. For cusps, choose
frames sigma,sigma' so sigma'^-1 diag(1,8) sigma=[[a_c,b_c],[0,d_c]]
with d_c>0 and a_c d_c=8. In those frames the cusp expansion of the
rescaled form is (d_c/8)^(-n/2) times the source cusp expansion at
(a_c tau+b_c)/d_c. The positive dilation preserves nonnegative exponents.
docs/coset-modularity.md supplies the exact branch conventions and normalization.
The resulting expansion at infinity uses exp(2 pi i tau/8). This proves
holomorphy. This proposition asserts neither cuspidality nor optimal level;
the rank-three coefficient cutoff is proved next.

## 8. An inclusive ternary identity cutoff

**Theorem 7.** Under the following rank-three hypotheses, complete cancellation
through the inclusive coefficient index B proves identical vanishing.

Let G be a positive definite integral 3x3 Gram matrix, a,b binary primal
and dual columns, b nonzero and a^t b even. Let K be any integral column
basis for L0={x:x^t b even}, G0=K^t G K, d0=det G0, and N0 the least
positive integer clearing the denominators of G0^-1. Then

    M=16N0,
    mu=[SL2(Z):Gamma0(M)]=product_(ell^e || M) ell^(e-1)(ell+1),
    B=mu/8.

If every signed norm4 coefficient c_N is exactly zero for 0<=N<=B,
then the characteristic theta function is identically zero. The constant
term and the upper endpoint are included. Missing occupied shells are
invalid evidence; unoccupied indices have coefficient zero once the full
sphere through B has been enumerated. No cuspidality is assumed.

Proposition 2.3 of Kane-Kim, with the conventions above, puts
F(z)=sum c_N exp(2 pi i N z) in weight 3/2 on Gamma0(M), character
chi_(4d0), and Theta(tau)=i^(a^t b) F(tau/8). This is the exact signed
doubled-coset difference, with a trivial coset based on 2L0. The Fourier
index in z is precisely our norm4 index, not norm4/8 or norm4/4.

Under gamma=[[r,s],[u,v]], its transformation factor is
chi_(4d0)(v) (u/v) epsilon_v^(-3) (uz+v)^(3/2).
The character and symbol have values +/-1 on the group, epsilon_v^4=1,
and the fourth power of the analytic factor is (uz+v)^6. Thus F^4 has
integral weight six and trivial character on Gamma0(M), including -I.
Its holomorphy in the upper half-plane is immediate.

For completeness at cusps, the weight-six slash by a cusp frame is
periodic with the cusp width and has polynomial growth as Im z grows,
uniformly in a period strip. For a negative Fourier index -m, its
coefficient is bounded by a polynomial times exp(-2 pi m Im z/w),
using the Fourier coefficient integral at height Im z. Letting this
height grow forces the coefficient to be zero. Hence every negative
index vanishes and F^4 is holomorphic at all cusps. Periodicity alone
would not exclude an essential singularity.

The complex identity theorem in [Brunault, Theorem 1](https://perso.ens-lyon.fr/francois.brunault/recherche/Sturm-bound-general.pdf)
applies to holomorphic integral-weight forms. On Gamma0(M), the infinity
width is one and -I belongs to the group, so a nonzero weight-six form
has order at infinity at most mu/2. This is the only sourced inequality
needed; all power and rounding steps here are deductions.

Since M is divisible by 16, its index factor from the prime two is
divisible by eight; therefore B is an integer. If F were nonzero with
all c_N=0 for 0<=N<=B, its first nonzero index m would be at least B+1.
The first coefficient of F^4 is c_m^4 at index 4m, which is nonzero
over C. Its order would satisfy 4m>=4(B+1)=mu/2+4, contradicting the
inequality. Thus F=0. The result is an exact identity test, not a
modular congruence test or a proof of minimal level.

The modular space and slash convention are those of
[Kane--Kim, sections 2.1--2.2 and Proposition 2.3](https://arxiv.org/html/2211.03987v2).
No spinor-genus theorem or unreviewed reduction in level is used.

This theorem does not assert optimal level or cutoff. Its mathematical
statement has no determinant restriction. The supplied cutoff implementation
accepts only rank three with determinant at most 24 and binary even pairs
with b nonzero. A modular identity certificate alone does not certify the
full stabilizer or its sign; those require the separate group evidence.

## 9. Current scope and reproduction

The unrestricted classification of P(L) in rank at least four remains open.
The symmetry converse has now been disproved in rank four: the companion
[manuscript](rank-four-counterexample.md) proves an all-level even zero with complete
original automorphism group {I,-I} and trivial stabilizer sign. Together with
the lower-rank proofs here, this establishes four as the first possible rank
of failure. No minimum-determinant or literature-priority claim is made.
Supplied witnesses of orders 4,8,16 do not establish their minimum possible
orders. Glue and arithmetic-provenance questions remain deferred.

The mathematical appendices are under `docs/`. Exact inputs are under `data/`
and preserved verifier reports under `reports/`. See [reproduction instructions](../docs/reproduction.md)
for the standalone standard-library implementation, both enumeration backends,
and the test suite. Finite tests check implementations; uniform statements rest
on the written proofs. Reviews are AI-assisted local audits, not journal peer review.

## References

1. Laura M. Chihara and Dennis Stanton, *Zeros of generalized Krawtchouk
   polynomials*, Journal of Approximation Theory 60 (1990), 43-57.
   [DOI](https://doi.org/10.1016/0021-9045(90)90072-X),
   [checked author PDF](https://www-users.cse.umn.edu/~stant001/PAPERS/krlaur.pdf).
   The polynomial citation uses equations (2.2), (2.3) and the remark after
   Theorem 3.4; page references in docs/literature-and-source-scope.md refer to that author PDF.
2. NIST, *Digital Library of Mathematical Functions*, Chapter 21,
   [definitions](https://dlmf.nist.gov/21.2),
   [functional properties](https://dlmf.nist.gov/21.3) and
   [modular transformations](https://dlmf.nist.gov/21.5), checked 2026-09-27.
   The formulas used are (21.2.5)-(21.2.6), (21.3.4), (21.3.6),
   (21.5.5) and (21.5.9).
3. Ben Kane and Daejun Kim, *Theta series of ternary quadratic lattice cosets*,
   [arXiv:2211.03987v2](https://arxiv.org/abs/2211.03987v2), 9 May 2024.
   Sections 2.1-2.2 and Proposition 2.3 supply the conventions and modularity
   used here; the ternary spinor-genus theorem is not used.

4. John H. Conway and N. J. A. Sloane, *Low-dimensional lattices. VI. Voronoi
   reduction of three-dimensional lattices*, Proceedings of the Royal Society
   of London A 436 (1992), 55-68. [Author text](https://neilsloane.com/doc/fedorov.pdf).
   Formula (5) in section 2.0.5, Theorem 8 and section 8.0.2 are used for
   classical geometric attribution; the signed classification is proved here.
5. Vitaliy Kurlin, *A complete isometry classification of 3-dimensional lattices*,
   [author PDF](https://kurlin.org/projects/lattice-geometry/lattices3Dmaths.pdf),
   [arXiv:2201.10543](https://arxiv.org/abs/2201.10543).
   Theorem 2.8 and Appendix A, Lemma A.1, supply the obtuse-superbase framework.
6. Francois Brunault, *Sturm bounds for general congruence subgroups*, 28 May 2021,
   [author note](https://perso.ens-lyon.fr/francois.brunault/recherche/Sturm-bound-general.pdf),
   Theorem 1. The fourth-power and inclusive-endpoint deductions are given above.
