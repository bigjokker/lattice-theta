# Rank-four decomposition, regressions and fibre identity

Prepared 2026-09-27 after `reviews/rank-four-census-review.md`. The completed census
and Walsh reduction are accepted. The following corollary is uniform; the
subsequent class counts and regression scan are finite-domain evidence.

## Uniform orthogonal-sum corollary

Every orthogonally decomposable positive definite integral rank-four lattice
satisfies the symmetry converse, at every determinant.

Proof. In an integral splitting basis G=G1 perpendicular G2, the signed
theta series factors as F1*F2. Absolute convergence justifies factorization.
If that holomorphic product is identically zero on the upper half-plane,
one factor is identically zero: if F1 is not identically zero, it is nonzero
on a nonempty open set, so F2 vanishes there and the identity theorem applies.
Both factors have ranks at most three. The already proved converse in ranks
one through three gives an epsilon-one characteristic-stabilizing isometry
on the vanishing factor. Extend it by the identity on the other factor.
Its dual shift is confined to the vanishing factor, so its epsilon is one
for the original pair. This proves the claim, also for more than two factors.

If an integral unimodular basis change U exhibits the splitting, transport
the characteristic as a'=U^{-1}a modulo two and b'=U^t b modulo two. Reducing
representatives changes only a nonzero overall phase. Conjugate the witness
back by U; stabilizing the characteristic and its epsilon are preserved.
The required lower-rank results are in `docs/rank-two-classification.md` and `docs/rank-three-classification.md`;
the factorization is in `docs/foundations-and-exact-enumeration.md`. The finite census is not part of this proof.

Hence the unresolved uniform rank-four converse can be restricted to
orthogonally indecomposable lattices. For a=0 the constant term is one, for
b=0 the coefficients are unsigned, and odd a dot b has -I as a witness.
Only even pairs with a and b nonzero and epsilon identically zero remain.

## A complete finite summand test

In rank four, a decomposition has a summand of rank one or two. If the
lattice determinant is Delta, the determinant d of that integral summand
divides Delta, since an integral splitting basis makes the full Gram block
diagonal. A line has a primitive generator of squared length d<=Delta.

A plane has a Gauss-reduced integral basis with Gram [[a,b],[b,c]],
1<=a<=c, |b|<=a/2 and d=ac-b^2. Its two basis norms are at most d.
Indeed if a=1 then b=0 and c=d. If a>=2, d>=3a^2/4 and
c=d/a+b^2/a<=d/a+a/4<=4d/(3a)<=d. Thus both basis vectors of an actual
plane summand lie in the sphere of squared norm Delta, and that pair is
saturated in the ambient lattice.

For a primitive vector v, Q(v) dividing all coordinates of Gv gives the
integral orthogonal projector v*(Gv)^t/Q(v). For a saturated independent
pair W=[v,w], let B=W^t G W. Its real orthogonal projector is

    P = W*B^{-1}*W^t*G.

It is integral exactly when B^{-1}W^t G is integral: saturation means that
integral vectors in the rational plane have integral coordinates in W.
An integral P yields the direct decomposition Z^4=im(P) plus ker(P).
Self-adjointness for G makes these two groups orthogonal. They have positive
ranks two and two. Conversely an actual orthogonal plane summand has this
integral projector. Saturation is checked by gcd of the six two-by-two
minors of W equalling one.

Therefore a complete sphere search through Delta, followed by all saturated
pairs, proves decomposition or indecomposability. The saved search uses the
larger common bound 24. It may stop once it finds an explicit integral
projector. A negative verdict requires finishing the sphere and every pair.
An interruption aborts the analysis before any new pack is written.

## Checked remainder of the census

`code/rank4_remainder.py` constructs the evidence using LDL spheres. Box replay
uses inverse-Gram coordinate boxes, repeats the complete summand test and
rebuilds the averaged shell data. The projection criterion and label-action
calculation are shared; this is independent vector enumeration, not a wholly
separate mathematical implementation. Every positive projector is additionally
checked by rational matrix multiplication, idempotence and G-self-adjointness.

| Census classes | Count |
|---|---:|
| A line summand found | 138 |
| No line, but a plane summand found | 4 |
| Orthogonally indecomposable | 27 |
| Total | 169 |

The 27 indecomposable classes carry 73 even symmetry zeros and 3,599 even
nonzeros; no pair is unresolved. Their groups all have order greater than
two. This last observation belongs only to this determinant domain.

Evidence: `data/rank4-remainder/det24.json`, bound to the accepted census
SHA-256. Replay: `reports/rank4-remainder/final-box.json`. The CLI requires
the exact accepted census bytes, rather than treating an arbitrary subset
or altered input as the completed domain.

## The six failed-first-shell regressions

For each indecomposable class the scan constructs P_A f exactly as an integer
group sum divided by the full group order, verifies the zero-average criterion
for all 256 pairs, and compares every epsilon-trivial even pair with complete
norm4<=32 shells. All six cancellations are on a=(0,0,1,1):

| Representative | Determinant | b | First supported norm4 | First nonzero norm4 | Coefficient |
|---|---:|---|---:|---:|---:|
| 131 | 15 | (0,1,1,1) | 6 | 14 | -2 |
| 131 | 15 | (1,0,0,0) | 6 | 14 | +2 |
| 131 | 15 | (1,1,1,1) | 6 | 14 | -2 |
| 132 | 22 | (0,1,1,1) | 7 | 15 | -2 |
| 132 | 22 | (1,0,0,0) | 7 | 15 | +2 |
| 132 | 22 | (1,1,1,1) | 7 | 15 | -2 |

All six averaged functions are nonzero. Direct evaluation of f(T^{-1}y)
on all 256 residues independently confirms their stored average numerators.
The complete first supported shell nevertheless sums to zero. These are
counterexamples to a proposed shortcut, not to the symmetry converse.

The proved calculation in `docs/rank-four-later-shell-lemma.md` explains all six and
extends their first nonzero formula to every integer s,t>=3. Exact control
replay is in `reports/rank4-remainder/later-shell-box.json`; the 21 cases
include parameter pairs outside the determinant-24 domain. All six new tests
pass in `reports/rank4-remainder/tests.txt`. These controls check coefficients
and failure policies, while the written inequality and shell calculation
prove the lemma for every parameter value.

## Exact fibre identity

For any rank-four basis write G=[[A,C],[C^t,D]] in two-by-two blocks, and split
a=(a1,a2), b=(b1,b2). A is positive definite. Put S=D-C^t A^{-1}C and
delta(w)=a1+A^{-1}Cw. With q=exp(2*pi*i*z), completing the square gives
the exact signed-series identity

    F(z) = sum_(w congruent a2 mod 2) (-1)^(((w-a2)/2) dot b2)
             * q^(w^t S w)
             * sum_(x in Z^2) (-1)^(x dot b1)
                 * q^((2x+delta(w))^t A(2x+delta(w))).

Fractional powers here mean exp(2*pi*i*z*r) for rational r. Together the
two exponents are the original integral norm4, and positive definiteness
justifies the double-series rearrangement.

The algebra was also checked exactly on 13,689 sample vectors across all
169 census Grams, including 219,024 sign factorizations, in
`reports/rank4-remainder/fibre-identity-check.json`. That finite regression
checks the formula's transcription; it does not replace the written identity.

This identity alone does not allow the lower-rank half-characteristic theorem
to be applied to the inner series: delta(w) can be rational and nonintegral.
In the regression, A=[[2,-1],[-1,2]] and
A^{-1}Cw=(-(2/3)(u+v),-(1/3)(u+v)). When u=v=1 the shift is
(-4/3,-2/3); when u=-v it is zero. Those different fibres can cancel at
their first combined shell. The later-shell lemma controls exactly this
interaction at the first two relevant levels.

The later-shell lemma and fibre identity were accepted in local review. The rank-four companion subsequently completes the scientific goal with a certified counterexample. Further structural classification remains deferred.
