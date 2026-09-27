# Doubled cosets and rational rescaling

Reviewed supplement, 2026-09-27. Sections 1-3 incorporate the comparison in
[reviews/theta-transformation-review.md](../reviews/theta-transformation-review.md); section 4 incorporates
[reviews/coset-rescaling-review.md](../reviews/coset-rescaling-review.md), with the explicit
valuation and cusp normalization below. No minimum modular level or
coefficient cutoff is asserted, and no census computation is changed.

## 1. Series and source conventions

Let L be positive definite and integral of rank n. Write a=2xi in L,
b=2delta in L*, assume p=<a,b> is even and b is not in 2L*, and put
L0={x in L:<x,b> even}. Choose t in L with <t,b> odd. Then a and 2t
belong to L0, and 2t does not belong to 2L0. For the distinct cosets

\[
C_0=a+2L_0,\qquad C_1=a+2t+2L_0,\qquad
\Phi_C(z)=\sum_{y\in C}e^{2\pi i z\|y\|^2},
\]

the identity from docs/foundations-and-exact-enumeration.md section 8 is

\[
\Theta_L[\xi,\delta](\tau)=i^p F(\tau/8),\qquad
F(z)=\Phi_{C_0}(z)-\Phi_{C_1}(z).
\tag{1}
\]

Indeed y=2x+a gives ||x+xi||^2=||y||^2/4. The unsigned Fourier variable
is q_z=exp(2 pi i z), so the denominator in the substitution is eight.

[Kane-Kim, section 2.1 and Proposition 2.3](https://arxiv.org/html/2211.03987v2)
use Q(x)=B(x,x), discriminant det G, and level the least positive integer
N for which N G^-1 is integral. Here B is our inner product. Their coset
cM+nu has nu in M and exact conductor c relative to cM. Their modularity
statement applies in every rank; the later spinor-genus theorem is ternary.

## 2. Bases and level scaling

If nu is outside 2L0, then 2L0+nu has conductor two with base L0.
If nu belongs to 2L0, the coset is the ordinary lattice 2L0, represented
with base 2L0 and conductor one. The distinct cosets cannot both be trivial.
For L=Z, a=0, b=1, the sets are 4Z and 2+4Z.

Let G be a Gram matrix of L0, d=det G and N=N_L0. Doubling the basis gives
Gram matrix 4G and discriminant 4^n d. Its level is exactly 4N:
put H=N G^-1, an integer matrix. The gcd g of its entries is one. Suppose
a prime ell divides g. It cannot divide N, since then N/ell would clear
the denominators, contradicting minimality. But ell^n divides det H and
d det H=N^n, so

\[
0=n\,v_\ell(N)=v_\ell(d)+v_\ell(\det H)\ge n,
\]

a contradiction. Therefore g=1.
Now s(4G)^-1=sH/(4N) is integral exactly when 4N divides s, since the
entries of H have gcd one. This proves the level assertion.

## 3. Common space in z

Put k=n/2 and define the source character

\[
\chi=\begin{cases}
\chi_{4d},&n\text{ odd},\\
\chi_{(-1)^{n/2}4d},&n\text{ even},
\end{cases}\qquad \chi_D(s)=\left(\frac{D}{s}\right).
\]

Conductor two on L0 gives group Gamma0(16N) intersect Gamma1(2).
Conductor one on 2L0 gives Gamma0(16N), with discriminant multiplied
by 4^n. On this group s is odd, and the square factor does not change
the character. Gamma1(2) is redundant here: lower-left divisibility by
16N and determinant one force both diagonal entries to be odd.
Consequently both unsigned summands and F lie in

\[
M_k(\Gamma_0(16N),\chi).
\tag{2}
\]

This adopts the corrected comparison in the submitted review. The difference
need not be cuspidal; cancellation at infinity would not establish that
condition at every cusp. There is no conclusion about symmetry from (2).

## 4. A sufficient congruence group after rescaling

**Rescaling proposition.** Set g(tau)=F(tau/8) and

\[
\mathcal H=\Gamma_0(16N)\cap\Gamma^0(8),\qquad
\Gamma^0(8)=\left\{
\begin{pmatrix}r&s\\u&v\end{pmatrix}\in\mathrm{SL}_2(\mathbb Z):8\mid s
\right\}.
\]

This congruence subgroup contains Gamma(16N). Its cusp width at infinity
is eight; it does not contain translation by one. We use the usual
modular-form definition allowing arbitrary cusp widths. Kane-Kim's section
2.2 explicitly requires translation by one in its group; that particular
convention cannot be copied to H without this qualification.

For gamma=[[r,s],[u,v]] in H, set

\[
\gamma'=\begin{pmatrix}r&s/8\\8u&v\end{pmatrix}.
\]

This is integral with determinant one and belongs to Gamma0(16N).
With z=tau/8 we have gamma'(z)=(gamma(tau))/8 and
(8u)z+v=u tau+v. The transformation of F at gamma' therefore gives the
transformation of g at gamma. For even n it is

\[
g(\gamma\tau)=\chi(v)(u\tau+v)^k g(\tau).
\tag{3}
\]

For odd n, use the half-integral slash convention of Kane-Kim section 2.2:

\[
(f|_k\gamma)(z)=
\left(\frac{u}{v}\right)\varepsilon_v^{2k}(uz+v)^{-k}f(\gamma z),
\quad
\varepsilon_v=\begin{cases}1&v\equiv1\ (4),\\i&v\equiv3\ (4).
\end{cases}
\]

Both matrices belong to Gamma0(4). Multiplicativity gives
(8u/v)=(8/v)(u/v), and the symbol's square on coprime entries is one.
The epsilon and analytic power factors agree because the denominator
u tau+v is identical. Therefore

\[
g|_k\gamma=\chi(v)\chi_8(v)g,\qquad
\chi_8(v)=\left(\frac{8}{v}\right),\qquad n\text{ odd}.
\tag{4}
\]

Use the principal branch with argument in (-pi,pi]. The same argument and
power occur on both sides of the conjugation, also for negative v. At -I,
the half-integral slash factor is i^n exp(-pi i n/2)=1; both chi(-1) and
chi_8(-1) are one in odd rank. In even rank chi(-1)=(-1)^k cancels
the factor (-1)^k. Thus the central element imposes no further restriction.
The character in (3) or (4) is well defined on this congruence subgroup.

The extra character in (4) accounts for the change in the theta multiplier.
Equations (3)-(4) hold for the signed theta after multiplication by i^p.
The chosen group is sufficient, with no assertion that it is optimal.

Holomorphy on the upper half-plane is immediate. For cusp holomorphy,
the positive rational dilation R=[[1,0],[0,8]] sends rational cusps to
rational cusps and intertwines gamma with gamma'=R gamma R^-1.
Given a cusp frame sigma in SL2(Z), choose a frame sigma' for R sigma(infinity)
and put

\[
U=\sigma'^{-1}R\sigma=
\begin{pmatrix}a_c&b_c\\0&d_c\end{pmatrix},\qquad a_c d_c=8.
\]

Replace sigma' by -sigma' if necessary to make d_c>0, hence a_c>0.
Write j(sigma,tau)=c_sigma tau+d_sigma and
phi(w)=j(sigma',w)^(-k) F(sigma' w), using the chosen analytic power.
The matrix identity gives j(sigma',U tau)=(8/d_c) j(sigma,tau), so the
positive scalar ensures the principal powers agree exactly and

\[
j(\sigma,\tau)^{-k}g(\sigma\tau)
=\left(\frac{d_c}{8}\right)^{-k}
\phi\!\left(\frac{a_c\tau+b_c}{d_c}\right).
\tag{5}
\]

The prefactor is a positive real constant. Nonnegative exponents in the
holomorphic cusp expansion of F remain nonnegative because a_c/d_c>0;
the translation changes only their scalar coefficients. Clearing rational
denominators gives a local parameter of suitable width. Thus there are
no poles at the cusps. This also preserves the polynomial-growth condition
used by the source, but does not imply vanishing at the cusps.

At infinity, if F(z)=sum_{m>=0} A_m exp(2 pi i m z), then

\[
g(\tau)=\sum_{m\ge0} A_m e^{2\pi i m\tau/8}.
\]

This is an expansion in the width-eight parameter, consistent with the
original norm4 coefficients. It need not be an integral-power expansion
in exp(2 pi i tau). No numerical cutoff has been derived here.

## Review and computational scope

The GL(n,Z) transformation audit is complete without correction. The common
z-space and tau-space deductions are now incorporated after the two focused
reviews. The old reports and
census artifacts are preserved. This note requires no new lattice search
and does not change any zero/nonzero verdict.

An exact denominator check confirmed level(4G)=4 level(G) on all 94 stored
census forms. Truncated rank-one and rank-two series, with L=Z and Z^2,
a=0, b=(1,0,...), tau=0.1+0.7i and gamma=[[171,8],[64,3]], agreed with
(3)-(4) to relative error below 6e-13. This matrix has chi_8(3)=-1, so
the odd-rank check exercises the extra sign. Truncation was |x|<=400 in
each rank-one factor. [The check record](../reports/coset-rescaling-checks.json)
is corroboration only, not an exact modularity or cusp certificate.
