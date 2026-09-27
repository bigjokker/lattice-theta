# The symmetry converse fails in rank four

Completed counterexample manuscript, 2026-09-27. The exact proof and certificate were independently accepted
in [accepted local review](../reviews/rank-four-counterexample-review.md). This is a local review manuscript;
no publication, priority or minimum-determinant claim is made.

## Abstract

The characteristic stabilizer of an integral lattice carries a sign character
whose nontriviality forces identical vanishing of its theta series. The
converse holds in ranks one through three. We give a primitive positive
definite integral rank-four lattice of determinant 23,716 with an identically
zero even half-characteristic, although its full integral automorphism group
is {I,-I} and the sign is trivial on both elements. Vanishing follows from
an exact norm-preserving bijection between two entire cosets of an ambient
lattice. An exhaustive 315-point calculation proves the full original group.
Thus rank four is the first possible rank of failure of the symmetry converse.

## 1. Definitions and normalization

Let G be a positive definite integral Gram matrix and a,b binary columns.
The primal shift has coordinates xi=a/2, and the dual shift has coordinates
delta=G^{-1}b/2. The original theta function is

    Theta_G[a,b](tau) = sum_(x in Z^n)
        exp(pi*i*tau*(x+a/2)^t G (x+a/2) + pi*i*(x+a/2)^t b).

An integral automorphism T satisfies T^t G T=G. The characteristic stabilizer
consists of such T with Ta=a modulo two and T^{-t}b=b modulo two. Its sign is
(-1)^epsilon(T), where

    epsilon(T)=a dot (T^{-t}b-b)/2 modulo 2.

An epsilon-one element forces identical vanishing. The symmetry converse
asks whether every identically zero half-characteristic has such an element.
Odd a dot b has -I as a witness. Our example has a dot b=0.

For y=2x+a define the signed coefficient series

    F(z)=sum_(y congruent a modulo 2) (-1)^(((y-a)/2) dot b)
            * exp(2*pi*i*z*y^t G y).

Then Theta_G[a,b](tau)=i^(a dot b)*F(tau/8). A norm4 index is y^t G y,
four times the squared norm in the original shifted-lattice sum. Positive
definiteness gives absolute convergence and finite coefficients at each index.

## 2. The counterexample theorem

**Theorem.** For

    G = [[12,12,-4,-1],
         [12,16,-6, 2],
         [-4,-6,28,-14],
         [-1, 2,-14,28]],
    a=(0,1,0,0),  b=(1,0,0,0),

G is primitive, integral and positive definite with determinant 23,716;
Theta_G[a,b] is identically zero; and Aut(G)={I,-I}, with epsilon zero
on both elements. In particular, the unrestricted rank-four symmetry
converse is false.

## 3. Positive definiteness and normalization

The exact LDL pivots of G are 12,4,77/3,77/4, all positive, and their product
is 4*77^2=23,716. The entries have gcd one. Also a dot b=0, so this is even
and its removed phase i^(a dot b) is one.

Write y=2x+a and Q(y)=y^t G y. The coefficient series is

    F(z)=sum_(x in Z^4) (-1)^x0 * exp(2*pi*i*z*Q(2x+a)).

The project's original theta function is F(tau/8). Positive definiteness
ensures absolute convergence throughout the upper half-plane.

## 4. Equality of two entire coset series

Define the following exact integer matrices:

    M = [[ 8,-4, 1, 2],
         [-4, 8,-3, 1],
         [ 1,-3,14,-7],
         [ 2, 1,-7,14]],

    V = [[-1,0,0,0],
         [ 1,2,0,0],
         [ 0,0,2,0],
         [ 0,0,0,2]],

    R = [[0,-1,0,0],
         [1,-1,0,0],
         [0,0,0,-1],
         [0,0,1,-1]].

Direct multiplication gives

    V^t M V = 2G,  det(V)=-8,
    R^t M R = M,  det(R)=1,  R^3=I.

Let L0={x in Z^4: x0 is even}, with basis K=diag(2,1,1,1). Under y -> Vy,
the positive and negative sign sets of F have images

    C_plus  = V(a+2L0)       = 2e2+4Z^4,
    C_minus = V(a+2e1+2L0)   = 2e1+4Z^4.

Here e1 is the first standard coordinate vector and e2 the second.
Indeed 2VK=4W, where

    W = [[-1,0,0,0],
         [ 1,1,0,0],
         [ 0,0,1,0],
         [ 0,0,0,1]],  det(W)=-1.

Also Va=2e2 and V(a+2e1)=-2e1+4e2, congruent to 2e1 modulo 4Z^4.
These identities establish the complete cosets, including every vector;
V is injective and W is unimodular.

The ambient norm is w^t M w/2, since Q(y)=(Vy)^t M(Vy)/2. As R is integral
unimodular and Re1=e2, it maps C_minus bijectively onto C_plus and preserves
this norm. Therefore the unsigned counts in the two cosets agree at every
norm. Their entire theta series agree, so F is identically zero. This proves
identity at all coefficients, without a truncated cancellation argument or
a modular cutoff.

The rotation does not swap the two cosets together: it next sends C_plus
to 2e1+2e2+4Z^4, a third coset. The equality is obtained in the ambient
lattice, while the symmetry converse asks for an isometry of the original
integral lattice stabilizing its characteristic.

For a direct basis check, the ambient rotation has original coordinates

    V^{-1} R V = [[ 1,   2,0, 0],
                 [-3/2,-2,0, 0],
                 [ 0,   0,0,-1],
                 [ 0,   0,1,-1]].

Its nonintegral entry prevents it from being an automorphism of the original
integral lattice. The full group argument below excludes every other possible
characteristic symmetry witness as well.

## 5. The full original group is only {I,-I}

An automorphism T is determined by its four columns v_i=T e_i. Their norms
are respectively 12,16,28,28 and their mutual products must be G_ij.
For Q(v)<=28, the bound v_i^2<=28*(G^{-1})_ii uses

    diagonal(G^{-1})=(4/11,3/11,4/77,4/77).

It confines every relevant integer vector to

    |v1|<=3, |v2|<=2, |v3|<=1, |v4|<=1.

There are just 315 integer points in this box. Substitution in Q gives these
complete lists, where each displayed vector includes its negative:

| Norm | Representatives up to sign |
|---|---|
| 12 | e1 |
| 16 | e2, 2e1-e2, 2e1-2e2 |
| 28 | 3e1-2e2, e1-2e2, e3-e1+e2, e3, e4, e3+e4 |

Thus v1=+/-e1. Multiplying T by -I lets us assume v1=e1. The required
products with v1 then force

    v2 in {e2,2e1-e2},
    v3 in {e3,e3-e1+e2},
    v4=e4.

Now <v2,e4>=G24=2 excludes 2e1-e2, whose product with e4 is -4.
Also <v3,e4>=G34=-14 excludes e3-e1+e2, whose product is -11.
The sole possibility is T=I. Reinstating the global sign gives exactly
I and -I. Both are integral automorphisms, so the group is complete.

This argument is also replayed exhaustively with two vector searches:
LDL sphere production and independent inverse-Gram boxes. Both give precisely
these two matrices. A decomposable rank-four lattice would have the separate
sign changes of two positive-rank summands and at least four automorphisms;
therefore this example is also orthogonally indecomposable.

## 6. The characteristic stabilizer has no negative sign

Both I and -I stabilize the binary labels. For I, epsilon is zero. For -I,

    epsilon(-I)=a dot ((-b-b)/2)=-a dot b=0 modulo 2.

These are all automorphisms. The identically zero series has no epsilon-one
characteristic-stabilizing automorphism, proving the proposition.

## 7. Relation to the completed lower-rank results

The uniform rank-one through rank-three converse in [lower-rank manuscript](low-rank-and-root-lattice-results.md)
remains valid. Together with the theorem above it implies that four is the
first rank at which the converse can fail. Those earlier proofs were accepted
before this stage; Astra's review checked this counterexample rather than
re-reviewing every lower-rank result.

Every decomposable integral rank-four lattice still satisfies the converse:
its theta series factors over positive-rank orthogonal summands, an identically
zero holomorphic product has a zero factor, and the rank-at-most-three witness
extends by the identity. Consequently a rank-four counterexample must be
indecomposable, as our group calculation also shows directly.

The exact determinant<=24 rank-four census contains 169 integral isometry
classes and resolves all 43,264 binary characteristics. Every zero in that
domain has a symmetry witness. The determinant-23,716 example is outside
that domain. The census does not set the smallest counterexample determinant;
we make no such claim. Its certificate and independent replays are listed in
[reproduction guide](../docs/reproduction.md).

The earlier parity-property and root-lattice classifications are unchanged.
An even zero with trivial stabilizer sign disproves the symmetry converse;
it does not classify all rank-four lattices or their even zeros.

## 8. Exact evidence and closing status

The all-level coset identities and the complete original group calculation
are independently replayed from `data/rank4-counterexample/det23716.json`.
The replay verifies lattice/metric/offset identities before the exhaustive
group search. Five adversarial counterexample tests include successful replay
with no positive-norm coefficient corroboration, demonstrating that truncated
cancellation is not the identity proof.

[accepted local review](../reviews/rank-four-counterexample-review.md) accepts the proof, all normalization
conventions, the separate 315-point enumeration and the full stabilizer signs.
Its optional display of V^{-1}RV is incorporated above. The historical research workspace retains the earlier proof checkpoint.

Reproduction commands, certificate hashes, census context and the review
inventory are in [reproduction guide](../docs/reproduction.md). The accepted scientific
goal is achieved through a rigorously certified counterexample. The closing manuscript and reproducible evidence are complete.
Potential minimum-determinant or higher-rank investigations remain deferred
in [deferred research questions](../docs/deferred-research.md).
