# Resolving the six cancelled-first-shell regressions

This is a proved coefficient lemma supporting the main rank-four question.
It handles the mechanism in Grok's six indecomposable regressions, uniformly
in two diagonal parameters. It does not prove the unrestricted converse or
classify every characteristic in this family. The deferred family-classification
side quest is not being promoted.

## Statement

For integers s,t>=3, put

    G(s,t) = [[ 2,-1,-1,-1],
              [-1, 2, 0, 0],
              [-1, 0, s, 0],
              [-1, 0, 0, t]],
    a = (0,0,1,1),  m = s+t.

The matrix is positive definite and has determinant 3st-2s-2t. For the three
binary columns below, the complete signed norm4 coefficient is zero at every
index N<m+8, and its first nonzero value is:

| b | First nonzero norm4 | Coefficient |
|---|---:|---:|
| (0,1,1,1) | s+t+8 | -2 |
| (1,0,0,0) | s+t+8 | +2 |
| (1,1,1,1) | s+t+8 | -2 |

In each case the minimal shell in the original characteristic support is m
and cancels. For the two census cases the group average also has support
at m, as the separate averaged-support scan verifies. The three pairs
are even. No claim about the automorphism group is needed for the coefficient
proof. Nonvanishing itself rules out an epsilon-one stabilizer witness.

The cases (s,t)=(3,3) and (3,4) are exactly census representatives 131 and
132 and account for all six regressions found by the averaged-support scan.
The norm4 first nonzeros are 14 and 15 respectively.

## Excluding every larger odd coordinate

A vector y congruent to a has the form y=(x0,x1,u,v), with x0,x1 even and
u,v odd. Write A=[[2,-1],[-1,2]] and w=(x0,x1). Completing the square in w
gives

    Q(y) >= s*u^2+t*v^2 - (2/3)*(u+v)^2.

The exact completed square has centre A^{-1}(u+v,0)^t. Since A is positive
definite and (u+v)^2<=2(u^2+v^2), the residual form is positive definite for
s,t>=3. This also proves positive definiteness of G(s,t).

If either |u| or |v| is at least three, u^2+v^2>=10 and

    Q(y)-m >= 3*(u^2+v^2-2) - (4/3)*(u^2+v^2)
            = (5/3)*(u^2+v^2)-6 >= 32/3 > 8.

Thus the entire sphere Q<=m+8 contains only u,v in {+1,-1}. This bound
excludes all large-coordinate vectors, rather than inspecting a subset of
the prospective witness shell.

## Exact two-variable shell calculation

Put x0=2k, x1=2l. If u=v=+1,

    Q=m+8*R,  R=k^2+l^2-kl-k.

If u=+1,v=-1,

    Q=m+8*H,  H=k^2+l^2-kl.

The antipodal cases give the same signs because a dot b is even. Both
H and R are nonnegative on integer pairs: H is positive definite, and
R=H(k-2/3,l-1/3)-1/3, so R>=-1/3 and its integer value is at least zero.
There are no intermediate norm indices between m and m+8.

The complete solutions at these two levels are:

| Form and value | Integer pairs (k,l) |
|---|---|
| H=0 | (0,0) |
| H=1 | (1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1) |
| R=0 | (0,0),(1,0),(1,1) |
| R=1 | (0,1),(0,-1),(2,1) |

Completeness follows from H=(l-k/2)^2+3k^2/4 and
R=(l-k/2)^2+3k^2/4-k. At value 0 or 1 these bounds respectively force
k in {-1,0,1} for H and k in {0,1,2} for R; substitution gives precisely
the listed solutions (k=2 supplies no R=0 solution).

For b=(0,1,1,1), the signs on the same-sign and opposite-sign cases are
(-1)^l and -(-1)^l. The Walsh sums of (-1)^l on R=0,H=0 are 1,1; on
R=1,H=1 they are -3,-2. The two antipodal copies therefore give
c_m=2*(1-1)=0 and c_(m+8)=2*(-3+2)=-2.

For b=(1,0,0,0), both cases have sign (-1)^k. The sums on R=0,H=0 are
-1,1, and on R=1,H=1 are 3,-2. Thus c_m=0 and c_(m+8)=2*(3-2)=+2.

For b=(1,1,1,1), the signs are (-1)^(k+l) and -(-1)^(k+l).
The sums are again 1,1 at level zero and -3,-2 at level one, yielding
c_m=0 and c_(m+8)=-2. Below m there are no supported vectors. This proves
the claimed complete coefficients and first nonzero indices in all cases.

## Role in the main proof

The failed shortcut was that a nonzero group average must have a nonzero
first supported shell. The lemma instead exhibits a cancelled initial slice
and an exactly controlled later shell. Its proof survives the regression
without assuming that equal norms imply isometries or that every invariant
theta function is independent.

The remaining challenge is an analogous mechanism for arbitrary indecomposable
rank-four G and epsilon-trivial pairs, which may have more interacting slices.
This calculation supplies a necessary regression and a working local model;
it does not supply the missing uniform structural reduction.
