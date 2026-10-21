# Stage 2: strategy review and proposed reductions

2026-10-21, Codex. This responds to Claude's proposed combined conjecture
after the completion of the [minus-type report](rank-four-minus-type.md) and the [census report](rank-four-minimum-determinant.md).
It leaves those reports unchanged.

**Status.** The combined conjecture and the minus-type conjecture remain
open. Sections 4-6 give new proof arguments for narrower reductions, written
out for review; they have not received independent proof review. The exact
checks in section 8 support the constructions but do not replace the proofs.

**Review, 2026-10-21.** Sections 2-6 were checked line by line, with an
independent computation, in [rank-four-stage-two-review.md](../reviews/rank-four-stage-two-review.md),
which also checks Lemma A and Theorem 5 of the minus-type report and runs the
plus-type test of section 1 (determinants 281 to 497). No errors were found.
This is an AI-assisted local review, not peer review.

## 1. Recommended target and order of work

Keep Claude's combined conjecture:

> Every zero in an even, orthogonally indecomposable, positive definite
> rank-four lattice of odd determinant has a symmetry witness.

By Theorem 5 this is equivalent to saying that every zero is one of the six
ordered pairs in the fixed plane of some automorphism with square -I. It
would imply both minus-type parity and the plus-type companion statement.
However, proving the minus-type part first is worthwhile even if the
plus-type part remains open.

The suggested computational test is appropriate, with two corrections:

* The new gap census already settles every even odd-determinant class in
  202..279, including plus type. The untested plus-type portion of the
  proposed bound is therefore **det = 1 mod 8 in 280..501**, namely 281..497
  (28 determinant values). Reuse the existing certificates below 280.
* Compare the **entire set of zero characteristics** with the union of the
  six-pair planes of **all** square roots of -I. Merely asking whether a
  lattice has some square root of -I could miss an additional unexplained
  zero on a lattice that also has Gaussian symmetry zeros. The existing
  characteristic-by-characteristic witness decision already avoids this
  mistake. Keep unresolved cases distinct from certified zeros.

No new census was run for this strategy review. Once this bounded test and
its independent replay are done, move to proofs rather than raising the
bound again. In particular, the next proof milestones should be:

1. Establish plane invariance correctly (section 2).
2. Review the mixed-kernel proof and its level-2 consequence (sections 4-5).
3. Review the exclusion of level 3 in minus type (section 6).
4. Attack the resulting restricted isospectrality problem (section 7).

## 2. Plane invariance is valid, but ordinary S is not its proof

Write D = det G, V = L/2L, Q(x) = x^t G x / 2 mod 2, and
B(x,y) = x^t G y mod 2. Put c = G^{-1}b mod 2. For independent a,c with
B(a,c)=0, the zero/nonzero status depends only on P = span(a,c).

The translation tau -> tau+1 sends (a,c) to (a,c+a), up to a nonzero scalar
factor in the characteristic theta function. The ordinary inversion
tau -> -1/tau instead changes G to G^{-1}; it does not by itself give a
swap inside the same lattice. There is a direct repair.

Since G is even of rank four with odd determinant, D = 1 mod 4, and
H = D G^{-1} is integral. Set delta = (3D+1)/2, an even integer. The integral
symplectic matrix

    Gamma = [[2I, 3G],
             [ H, delta I]]

preserves the scalar family of period matrices:

    tau G -> ((2 tau + 3)/(D tau + delta)) G.

Indeed 2 delta - 3D = 1 and GH = HG = D I. The characteristic
transformation gives, modulo 2,

    a' = delta a - H b = c,
    b' = -3G a + 2b = G a.

The diagonal correction terms vanish modulo 2: they are the diagonals of
delta H and 6G. Thus (a,c) -> (c,a), with a nonzero automorphy factor.
Together with translation, this generates all six ordered bases of P.

This is a deduction from the general theta transformation formula,
[DLMF 21.5.9](https://dlmf.nist.gov/21.5.E9), rather than an application of
ordinary S on the same lattice. Integer changes of characteristic only
change a nonzero phase and do not affect vanishing.

In a four-dimensional minus-type quadratic space, an isotropic plane for B
cannot be totally singular for Q. Since Q restricted to that plane is a
nonzero linear form, it has exactly one singular nonzero vector.

## 3. A required control for the proposed cusp argument

Consider the decomposable even lattice

    G = [[2,1,0,0], [1,2,0,0], [0,0,2,1], [0,0,1,4]],  det G = 21,
    a = (1,0,1,0),  c = (0,1,0,1),  b = Gc mod 2 = (1,0,1,0).

Here B(a,c)=0 and

    Q(a)=0,  Q(c)=1,  Q(a+c)=1.

The characteristic is even globally, but each rank-two factor has an odd
characteristic. Their product vanishes identically. Equivalently,
diag(-1,-1,1,1) is a stabilizer witness with epsilon=1.

Thus a zero and a plane with exactly one singular vector really can
coexist. In particular, the equalities obtained by splitting

    (a+2x)^2 = a^2 + 4(a,x) mod 8

cannot contradict the single-singular-vector condition alone. This example
must pass every necessary Gauss-sum or cusp test for vanishing.

It does not refute the conjecture: it is decomposable. It shows precisely
what the proposed proof still needs: a step that uses **global orthogonal
indecomposability**. A calculation depending only on the displayed
mod-2 plane cannot provide that step. More refined cusp calculations may
still help, but matching local data is a necessary condition, not a proof
that the relevant global lattices are isometric.

## 4. Proposed lemma: mixed order-four kernel maps admit swaps

**Lemma.** Let M be an even positive definite rank-four lattice with
v_2(det M)=2. Let g in Aut(M) have order four and g^2 != -I. If a class
v in M/2M has an orbit of length four under g, there is an involution
R in Aut(M) with Rv=gv. If its orbit has length two, g itself swaps the
two classes.

This upgrades the finite observation in section 3 of the minus-type report.

**Proof.** Put h=g^2. Its rational eigenspaces have dimensions two and two:
g is a quarter-turn on the minus eigenspace and an involution r on the
plus eigenspace. Let N and P be their intersections with M. They are
primitive, orthogonal, even binary lattices.

The lattice N is a torsion-free rank-one module over Z[i], hence free.
Choose e,f with f=ge. Orthogonality of g gives

    Gram(N) = c I_2,   c positive and even.

The glue group M/(N perpendicular P) is elementary abelian of order 2^k:
twice each vector of M lies in N perpendicular P. Its projections into the
discriminant groups A_N and A_P are injective, and land in their 2-torsion.
Thus 0 <= k <= 2. Put s=v_2(c) >= 1 and t=v_2(det P) >= k. Determinants give

    2s + t - 2k = 2.                                      (1)

The case k=2 is impossible. Equation (1) leaves two possibilities:

* (s,t)=(1,4). The projected glue in A_N is all its 2-primary part,
  (Z/2)^2, with nondegenerate discriminant bilinear pairing. Integrality of
  M makes its image in A_P another nondegenerate copy of (Z/2)^2. It splits
  off orthogonally from the 2-primary part of A_P, leaving a group of order
  four. This requires at least three generators, whereas the discriminant
  group of a binary lattice has at most two.
* (s,t)=(2,2). The bilinear pairing on A_N[2] is identically zero, since
  c/4 is integral. Its image would be the whole 2-primary part of A_P,
  whose pairing is nondegenerate. This contradicts the pairing relation
  imposed by integral glue.

Therefore k <= 1. If k=1, equation (1) forces s=1,t=2. Stability under g
forces the N component of the nonzero glue element to be the unique
nonzero quarter-turn-fixed class in N[2], namely (e+f)/2 modulo N. Its
P component is fixed by r.

Let kappa(e)=e, kappa(f)=-f. Both

    R_plus = kappa perpendicular r,
    R_minus = -kappa perpendicular r

are involutions preserving M: they preserve N perpendicular P and fix the
glue class. This also holds when k=0.

If x belongs to N perpendicular P, g^2 x-x=-2x_N lies in 2M, so its class
does not have a four-cycle. Otherwise k=1 and

    x_N = alpha e + beta f,   alpha,beta in Z + 1/2.

The P components of R_plus x, R_minus x and gx agree. On N the differences
are

    R_plus x_N - gx_N = (alpha+beta)(e-f),
    R_minus x_N - gx_N = (beta-alpha)(e+f).

The integers alpha+beta and beta-alpha have opposite parity. One difference
therefore lies in 2N, proving that one of the two involutions sends the
class of x to gx. Since it is an involution, it swaps the pair. QED.

## 5. Proposed consequence: level 2 is impossible at odd nonsquare determinant

Let L be even of rank four with odd nonsquare determinant d. For an even
characteristic set

    M = {x in L : (b,x) is even},
    xi = a/2,   t in L minus M,
    O1 = M + Z xi,   O2 = M + Z(xi+t).

Here (b,x) means the coordinate pairing b^t x. All three of L,O1,O2 are
index-two extensions of M, and det M=4d. Write v1,v2 for their doubled
coset representatives in M/2M, and w=v1+v2 for the representative of L.

**Claim.** If any element of Aut(M) maps O1 to O2, some element swaps
O1 and O2. Thus there are no strict level-2 configurations in this setting,
whether or not L is indecomposable.

The 2-primary discriminant group of M has order four. It has at most two
isotropic subgroups of order two, hence M has at most two even index-two
extensions. To see the bound, a cyclic group of order four has only one
subgroup of order two. In the elementary abelian case, if all three
nonzero elements were isotropic, the discriminant pairing would vanish
identically, contradicting nondegeneracy.

Suppose gv1=v2. If g fixes L, it fixes w, so gv2=v1 and we are done.
Otherwise g exchanges the two even extensions of M. Write their classes
as w,u. Then gw=u, gu=w, and the orbit of v1 is

    v1, v1+w, v1+w+u, v1+u.

Choose an odd power of g that has pure 2-power order and whose exponent is
1 mod 4. It still sends v1 to v2. Its order is four or eight: a finite-order
rational representation of dimension four cannot contain a primitive
16th root of unity.

An order-eight automorphism has characteristic polynomial Phi_8, so its
square has square -I. An order-four automorphism whose square is -I is
already a complex structure. Both cases are impossible here: an integral
rank-four lattice with such a structure is a free rank-two Z[i]-module;
its Gram determinant is the square of the determinant of an integral
Hermitian matrix. But det M=4d is not a square.

The remaining automorphism is mixed order four. The lemma in section 4
supplies an involutive swap. This proves the claim.

A swap is exactly a level-1 witness by the project's kernel criterion.
Combining this claim with Theorem 5, an indecomposable even rank-four
lattice of odd nonsquare determinant has **no level-1 or level-2 zeros**.
In particular this applies to minus type. It also applies to plus type
when the determinant is nonsquare; the square-determinant plus case is
not covered by this argument.

## 6. Proposed consequence: minus type has no level-3 zeros either

Assume now L is even, indecomposable, rank four, and d=5 mod 8. The case
a=c is nonvanishing by Theorem 1 of the parity report, so take independent
a,c. At least one of Q(a),Q(c) equals one.

If Q(a)=1, every vector in O1 minus M or O2 minus M has half-integral norm:

    (a/2+x)^2 = a^2/4 + (a,x) + x^2,   x in L.

Conversely the vectors of M have integral norms. Hence M is intrinsically
recognizable inside either Oi as the set of vectors of integral norm.
Every isometry O1 -> O2 must preserve M. Section 5 then turns it into a
swap, contradicting Theorem 5 for the indecomposable minus-type lattice L.
Thus O1 and O2 are not isometric in this case.

If Q(a)=0, then Q(c)=1. Dualize the entire configuration, not just its theta
identity. Put S=L+Z(a/2). The common sublattice of O1^*,O2^* is S^*, their
joint superlattice is M^*, and the third intermediate lattice is L^*.
In the dual basis of L^*, the characteristic is (b,a), up to signs modulo 2.
Scale its bilinear form by d. The resulting third lattice has Gram matrix

    H = d G^{-1},   det H=d^3=5 mod 8.

This is integral and even: each diagonal cofactor of G is the determinant
of a 3-by-3 even symmetric matrix, hence is even. Duality and scalar
rescaling preserve orthogonal indecomposability. If b=Gc+2z, then

    b^t H b / 2 = d c^t G c / 2 + 2d c^t z + 2 z^t H z,

so Q_H(b)=Q_G(c)=1. We are in the preceding case for this dual, scaled
configuration. An isometry O1 -> O2 would give an isometry O1^* -> O2^*,
which that case excludes. Therefore O1 and O2 are not isometric here either.

**Result of the proposed reduction.** For every relevant characteristic of
an even indecomposable minus-type lattice, the two associated overlattices
are non-isometric. Consequently **any hypothetical zero must have level 4**.
This is a restriction on possible counterexamples, not a proof that there
are none.

## 7. The remaining hard problem, stated precisely

With the reductions above, the minus-type conjecture becomes:

> In a Klein configuration M subset L,O1,O2 subset S arising from an even,
> indecomposable rank-four lattice L with det L=5 mod 8, the non-isometric
> pair O1,O2 cannot have equal theta series.

All the configuration data matter: [L:M]=[O1:M]=[O2:M]=2, the three
extensions are distinct, and the third lattice L has the specified parity,
determinant and indecomposability. General rank-four theta rigidity is
false, so a proof must use these additional restrictions.

Plane invariance allows choosing Q(a)=1 for the analytic attack. Then the
integral-norm part of each Oi is exactly M, and all the extra vectors have
half-integral norms. The target is equality of the two remaining coset
theta series. That supplies a much narrower starting point than all
rank-four isospectral pairs.

Two useful proof directions, with explicit stopping conditions:

* **Relative rigidity.** Try to reconstruct enough of the embedding of M
  from the two equal coset spectra to force a kernel-preserving isometry,
  or force L to split. The already proved rank-three rigidity suggests
  looking for a canonical lower-rank slice or orthogonal projection. Such
  a slice must be recovered from the theta data; assuming one exists would
  beg the question.
* **Modular decomposition.** Derive the two index-four coset equalities and
  calculate their local/cusp constraints, keeping the det-21 control in
  every calculation. Identify explicitly what survives after those
  constraints: a nontrivial global coefficient or cusp-form component
  that must be shown nonzero. If all proposed obstructions are also zero
  on the control and never use indecomposability, they have not reached
  the decisive step.

A useful partial result would already be a theorem for an identifiable
subfamily, such as lattices for which a shortest shell determines a
canonical lower-rank slice. Any such subfamily needs a precise hypothesis
and proof; it should not be described as covering the general case.

For plus type, keep the combined conjecture as the organizing question,
but do not assume these minus-type reductions settle it. Totally singular
planes and odd square determinants are genuine additional cases.

## 8. Supplementary exact checks

Run:

    python code/rank4_stage2_checks.py --output scratch/rank4-stage2/strategy-checks.json

Preserved output: `reports/rank4-stage2/strategy-checks.json`, including hashes of the
two source censuses. The checks:

* Verify the integral symplectic matrix and its swap on all 105 relevant
  characteristics of all 300 stored even odd-determinant classes through
  201; also check the dual quadratic-value identity (31,500 pairs).
* Construct the two reflections in section 4 from saturated eigenlattices,
  rather than merely searching an automorphism list for a lucky swap. Test
  every mixed order-four element and every four-cycle in the 32 stored
  kernels with v_2(det M)=2 through determinant 64. This includes the 540
  elements and 3,936 four-cycle steps in the original observation.
* Add the computed det-84 sign kernel of the det-21 control, and check the
  mapping-to-swapping consequence on the nonsquare-determinant kernels in
  this finite universe. This supplies a nonvacuous decomposable control.
* Verify the det-21 characteristic and its epsilon-one witness exactly.

These are checks of algebra, constructions and small examples. They do not
independently verify the new general proofs or decide the conjecture.
