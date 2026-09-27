# Uniform rank-two parity classification

Reviewed and incorporated 2026-09-27 from the argument in
[reviews/initial-mathematical-review.md section 2](../reviews/initial-mathematical-review.md#2-rank-two-the-parity-property-is-indecomposability).
The proof below makes the unequal-axis shell estimate and equal-summand reduction
explicit. It uses no determinant bound or finite census assumption.

**Theorem.** A positive definite integral lattice of rank two has nonzero theta
for every even-p half-integral characteristic if and only if it is orthogonally
indecomposable. In a decomposable rank-two lattice exactly one of the ten even
characteristic classes vanishes. Every vanishing characteristic in rank one or
two has an epsilon-one witness on its characteristic stabilizer. Thus the
symmetry converse holds in these ranks at every determinant.

Here vanishing means identical vanishing as a function of tau in the upper half
plane. Coordinates and signed coefficients are those of docs/foundations-and-exact-enumeration.md sections 1–4.

## Reduced bases and intrinsic decomposition

Binary reduction gives a basis with

\[
G=\begin{pmatrix}A&B\\B&C\end{pmatrix},\qquad
1\le A\le C,\quad 0\le2B\le A.
\]

Subtracting multiples of the first vector from the second makes |2B|<=A.
If C<A, interchange vectors and repeat; each interchange strictly decreases
the positive integer A. A sign change makes B nonnegative. All steps are
unimodular, so characteristic coverage and vanishing are preserved.

Write N(m,n)=Am^2+2Bmn+Cn^2. Its minimum on nonzero integer vectors is A.
For mn>=0 or a zero coordinate this follows directly. For mn<0 set
u=|m|, v=|n|. Reduction gives

\[
N(m,n)\ge A u(u-v)+Cv^2.
\]

If u>=v this is at least Cv^2>=A. If u<v, use C>=A to obtain
N>=A(u^2-uv+v^2); since v>=2, this is at least 3A, by
u^2-uv+v^2=(u-v/2)^2+3v^2/4. The vector (1,0) attains A.

We claim that B=0 is equivalent to orthogonal decomposability. One direction
is immediate. For the converse let L=Ze perpendicular Zf with norms d<=h.
The shortest vectors are +/-e, and also +/-f if d=h. Since the first vector of
any reduced basis is shortest, it is one of these. Relabel the two orthogonal
axes if necessary and change signs so that the first vector is e. The second
has the form xe+sf with s=+/-1, by the basis determinant. Its pairing with e
has absolute value |x|d. The reduced bound |2B|<=A=d forces x=0. This argument
also covers d=h; it does not assume the second reduced vector is shortest.
Consequently every reduced basis of a decomposable lattice has B=0, and every
indecomposable lattice has B>0.

## The ten even pairs

Take a,b in {0,1}^2 in primal and dual coordinates, with p=a dot b even.
Four pairs have a=0: their unique norm-zero vector has signed coefficient 1.
Three further pairs have b=0: all signs are positive, so their series are
nonzero. This leaves exactly

\[
(a,b)=((1,0),(0,1)),\quad((0,1),(1,0)),\quad((1,1),(1,1)).
\]

For y=2x+a, norm4=N(y), and the sign after removing i^p is (-1)^(x dot b).
The following table gives the complete minimal shells for those pairs.

| a | b | Minimal norm4 | Complete minimal vectors y | Signed coefficient |
|---|---|---|---|---:|
| (1,0) | (0,1) | A | +/-(1,0) | 2 |
| (0,1) | (1,0) | C | +/-(0,1) | 2 |
| (1,1), B>0 | (1,1) | A+C-2B | (1,-1), (-1,1) | -2 |
| (1,1), B=0 | (1,1) | A+C | (+/-1,+/-1), all four choices | 0 |

Here is the completeness proof; in particular the two axis estimates do not
assume that exchanging axes preserves A<=C.

**Odd/even coset.** If n=0, the least norm is A at m=+/-1. If |n|>=2,
the same-sign case has N>=4C>A. In the opposite-sign case use the bound above.
For u>=v it gives N>=Cv^2>=4C; for u<v it gives
N>=A((u-v/2)^2+3v^2/4)>=3A. Thus no other vector has norm <=A.
Both displayed vectors have x_2=0 and therefore sign +1.

**Even/odd coset.** If m=0, the least norm is C at n=+/-1. Otherwise u>=2
and v>=1. For mn>=0, N>=4A+C>C. For mn<0,

\[
N-C\ge A u(u-v)+C(v^2-1).
\]

If u>=v, their different parities imply u>v, so the right side is positive.
If u<v, then v>=3 and C>=A gives
N-C>=A(u^2-uv+v^2-1)>=A(3v^2/4-1)>0.
Thus these are the only minimal vectors. Both have x_1=0 and sign +1.

**Odd/odd coset.** Put M=A+C-2B. Same-sign vectors have N>=A+C,
strictly above M when B>0. For opposite signs u,v are positive odd integers and

\[
N-M=A(u^2-1)-2B(uv-1)+C(v^2-1)
\ge A u(u-v)+C(v^2-1).
\]

If u>=v, this bound is nonnegative, and is zero only when u=v=1.
If u<v, v>=3 and the same estimate as above gives a strictly positive bound.
The only vectors at M for B>0 are therefore (1,-1) and (-1,1). They give
x=(0,-1) and (-1,0), so both signs are -1. If B=0, independent squares give
exactly four vectors at A+C, with two positive and two negative signs.

## Conclusions and symmetry witnesses

When B>0 all ten even pairs have a complete nonzero coefficient, proving the
positive direction. When B=0 the only remaining pair is a=b=(1,1).
Its identical vanishing follows from the two odd rank-one factors, not from
first-shell cancellation. Equivalently T=diag(-1,1) preserves G, fixes both
classes modulo 2, and has dual shift (-1,0); hence epsilon(T)=-1=1 modulo 2.
The other nine even pairs are nonzero by the arguments above. This proves both
the negative direction and the exact count of one even zero.

For rank one, the even pairs are (a,b)=(0,0),(0,1),(1,0): the first two have
constant coefficient 1, and the last has complete minimal coefficient 2 at
norm4=d for G=[d]. All odd pairs in either rank vanish by -I. Every zero in
ranks one/two therefore has a symmetry witness, while an epsilon-one witness
always forces a zero by docs/foundations-and-exact-enumeration.md section 2. This proves the symmetry equivalence
in these ranks without requiring a full-group classification.

The proof permits A=C, 2B=A, odd entries and nonprimitive forms. Its special
cases include A2 and its positive integer scalings. The determinant<=24 census
is now an independent finite regression of this uniform theorem, rather than
a hypothesis of the result. Historical census artifacts remain unchanged.
