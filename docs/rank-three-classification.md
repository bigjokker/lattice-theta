# Uniform rank-three parity classification

New structural phase, 2026-09-27. The proof below is self-contained for
positive definite integral lattices. Grok's adversarial review accepts the
theorem; its requested proof insertions are incorporated below. Independent
certificate replay and source-check limits are recorded separately; no
novelty claim is made. The old manuscript checkpoint remains unchanged.

**Theorem.** A positive definite integral rank-three lattice has nonzero
theta for every even half-integral characteristic if and only if it is
orthogonally indecomposable. Every identically zero half-integral
characteristic in rank three has an epsilon-one stabilizer element.
Both statements hold at every determinant. Nonzero means not identically
zero as a function of tau, rather than nonzero at each individual tau.

## Exact orthogonal-decomposition test

A rank-three lattice is decomposable if and only if it has a rank-one
orthogonal direct summand Zv. Put d=v^t G v. Such a summand exists exactly
when v is primitive and every entry of Gv is divisible by d. Indeed the
integral row lambda=(Gv)^t/d satisfies lambda v=1, and

    x=(lambda x)v + (x-(lambda x)v)

is an integral orthogonal direct-sum decomposition for every x. Conversely
an orthogonal summand's primitive generator has these divisibilities.

Every possible summand generator has d<=max_i G_ii in any supplied basis.
Write each basis vector as m_i v+w_i. At least one m_i is nonzero, since
the basis spans L; then G_ii=m_i^2 d+||w_i||^2>=d. Thus enumerating the
complete sphere through max G_ii and testing primitive columns proves
both positive and negative decomposition results. An interrupted sphere
cannot prove indecomposability.

To exhibit a saturated complement, extend v to a unimodular basis by
Bezout row operations, then replace each other column u by
u-(lambda u)v. The resulting unimodular basis gives diag(d,H2).
Binary Gauss reduction of H2 distinguishes a sum of three lines from
a line plus an indecomposable binary lattice, by RANK-TWO.md.

The uniform decomposable cases already follow from the reviewed binary
theorem: a sum of three lines has 3^3=27 even nonzeros and nine even
zeros, giving 27+9=36. A line plus an indecomposable binary lattice
has 3*10=30 even nonzeros and 1*6=6 even zeros, giving 30+6=36;
the zeros have odd characteristics on both factors. Products of nonzero
factor series are nonzero because their least nonzero coefficients
multiply to a nonzero coefficient at the sum of their indices. Every zero has
a witness: negate an odd factor; for an even zero in a decomposable
binary factor use its known coordinate flip. Extend by identity on
the other factors. This proves the symmetry converse for every
decomposable rank-three lattice, without a determinant bound.

## 1. Obtuse superbases by a terminating integer reduction

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

## 2. The positive-edge graph

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

## 3. Characteristics and a complete lower bound

A primal class is a parity vector a on the four vertices modulo adding
the all-one vector. After choosing root r with ar=0, fix xr=0 and require
xi congruent ai mod 2. A dual functional is an integer vector beta with
sum beta_i=0, where beta_i=<vi,b>. Its sign is

    (-1)^(sum_i beta_i (xi-ai)/2).

If the vector sum ai vi differs from the original primal representative
by 2ell, the displayed sign differs from the original engine's sign by
the constant (-1)^<ell,b>. Changing the chosen lift of a therefore only
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

## 4. All indecomposable cases

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
hypothesis. The certificate's selected_vectors field records the vectors
determining the sign; independent replay sums the entire shell.

For example, G=[[2,0,-1],[0,2,-1],[-1,-1,6]] has determinant 20 and
is census class 97 (zero-based indexing). For a=(1,1,0), its norm4=4
shell has coefficient zero when b=(1,1,0) or (1,1,1). At norm4=20,
the full shell has ten vectors; its signed coefficients are respectively
+2 and -2. Eight of those vectors cancel on the k=0 slice, leaving the
two contributing vectors. Here h=4 and m=4.

## Consequences

The four cases prove that every even characteristic of an indecomposable
rank-three integral lattice has a complete nonzero coefficient. Orthogonal
sums fail P by the established product obstruction. Therefore

    P(L) holds in rank three iff L is orthogonally indecomposable.

Combining this with the decomposable-factor argument above proves the
symmetry converse for all even characteristics in rank three, at every
determinant. Every odd characteristic has the witness -I, with
epsilon(-I)=-a^t b=p mod 2, by docs/foundations-and-exact-enumeration.md section 2. Thus the converse
holds for ALL half-integral characteristics in rank three. The graph
reasoning also implies that a positive-edge graph
without articulation cannot represent a decomposable lattice; otherwise
its proved P would contradict the product obstruction. This last
implication is a consequence, not an assumed graph classification.

There is no extrapolation from the 120-class census in this argument.
The shell index m or m+4h is computed in an obtuse superbase and is the
original norm4 index. Degenerate conorms, equal norms, arbitrary positive
integer scalings and boundary symmetries are included explicitly.

## Source and review scope

The elementary reduction and graph proof above are deductions written in
this phase. They are consistent with the known low-dimensional obtuse-
superbase framework; a preliminary primary-author source is
[Conway, The Sensual (Quadratic) Form, third lecture](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/conwaysens.pdf).
No existing theta-classification theorem was identified in the preliminary
search, and that search does not establish novelty. Conway's printed
pp.69--76 were inspected for the obtuse-superbase framework; the termination
proof used here is the explicit integer energy descent above.
[Grok's completed review](../reviews/rank-three-proof-review.md) accepts the
uniform argument and checks Conway's printed pp.69--76 and Kurlin's
Theorem 2.8 / Appendix A. Those sources support the superbase framework,
not this signed-characteristic classification. The report found no
classification statement in the sources it opened; this limited search
does not establish novelty. Selling's original paper and Delone's book
were not read in that review.

## Verification record and current status

The uniform argument is now reviewed and its requested insertions are
incorporated. The finite structural claims below are independently
certified. Computation checks the supplied instances; it does not replace
the determinant-independent proof in sections 1--4. The saved initial
box report's symbolic_theorem_externally_reviewed=false records its
original run before the review; that historical report remains unchanged.
The current acceptance and artifact bindings are recorded separately in
[the review resolution](../reports/census-rank3/structure-review-resolution.json).

- [Structure pack](../data/census-rank3/structure-det24.json.gz): all 120
  existing isometry classes; 51 sums of three lines, 44 line-plus-
  indecomposable-binary sums, and 25 indecomposable lattices. Every split
  has an integral projection row and saturated unimodular block basis.
  All 7,680 pairs are explicitly recorded, including the 28 odd pairs
  per class. There are 723 even zeros plus 3,360 odd zeros and 3,597
  even nonzeros. Every zero has an explicit checked witness.
- [Independent box report](../reports/census-rank3/structure-det24-box.json):
  all decomposition spheres, bases, factor transports, Selling traces,
  characteristic coverage and signed coefficients replay successfully;
  no limits. The 25 positive graphs are 12 cycles, 11 diamonds and two
  K4s. All 900 positive even pairs have exact symbolic shell claims.
  All 22 cancelled-first-shell cases use the diamond-later coefficient.
  Source-domain coverage relies on the previous complete census replay;
  this report checks that the same complete class list is covered and
  binds its exact source-pack SHA-256, rather than regenerate lattices.
- [Tests](../tests/test_rank3_structure.py): all seven pass. They include equal and
  unequal conorms, positive scaling, noncanonical bases, primitive Bezout
  completion, zero/interrupted searches, and altered certificates. The
  existing determinant-29 first-shell trap is an isolated regression:
  its complete coefficient at norm4=11 is -2 in the diamond-later case.
  Other isolated scaled graph controls test the uniform argument. The new
  explicit collision regression checks both the h=4 ten-vector shell
  and the h=8 six-vector shell, each with only two selected vectors and
  complete coefficient +2, as well as both signs in census class 97.
  No larger determinant census was conducted.

Structure pack SHA-256:

    ae8288b4f1ae0bbddf3f6161200ea9e7e7c52eab28df772146b84aafb1a762f6

Source census SHA-256:

    1083147c76d6b66aeed3491d35bcb1f83ebed65c9cc0c630e95eaac90b32816a

Print-only replay and tests:

```powershell
python code/verify_rank3_structure.py data/census-rank3/structure-det24.json
python scripts/run_tests.py --all
python scripts/verify_integrity.py
```

Rebuild to a fresh destination:

```powershell
python code/build_rank3_structure.py --output reports/rebuilt-rank3/structure-det24.json
```

The determinant-independent proof and structural certificates are incorporated in the lower-rank companion manuscript.
