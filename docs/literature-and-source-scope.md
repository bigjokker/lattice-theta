# Checked sources and remaining literature gaps

Updated 2026-09-27. [reviews/initial-mathematical-review.md](../reviews/initial-mathematical-review.md) and
[reviews/modularity-follow-up.md](../reviews/modularity-follow-up.md) and
[reviews/theta-transformation-review.md](../reviews/theta-transformation-review.md) and
[reviews/coset-rescaling-review.md](../reviews/coset-rescaling-review.md) are the submitted reviews.
This note records the sources independently checked during incorporation,
and distinguishes source statements from deductions made in this project.
No unsuccessful search establishes novelty.

## Chihara–Stanton: checked polynomial citation

Laura M. Chihara and Dennis Stanton, *Zeros of generalized Krawtchouk polynomials*,
Journal of Approximation Theory 60 (1990), 43–57,
[DOI](https://doi.org/10.1016/0021-9045(90)90072-X),
[author PDF](https://www-users.cse.umn.edu/~stant001/PAPERS/krlaur.pdf).
The downloaded [local copy](https://www-users.cse.umn.edu/~stant001/PAPERS/krlaur.pdf) was checked by
text extraction and visual inspection of author-PDF pages 2 and 5.

Equation (2.2), with q=2 and length 2k, identifies their k_k(s,2,2k) with
the coefficient [t^k](1-t)^s(1+t)^(2k-s) in docs/foundations-and-exact-enumeration.md section 6.
Equation (2.3) gives leading coefficient (-2)^k/k!. The remark immediately
after Theorem 3.4 gives the k odd roots 1,3,...,2k-1. The product formula follows.
The source citation is therefore verified; the An shell evaluation is a separate
lattice argument. Only author-PDF page numbers are used here.

## Parity and the standard Riemann theta specialization

Osorio–Vázquez-Mozo, *Quantum corrections in two-dimensional non-supersymmetric
heterotic strings*, [arXiv:hep-th/9511149v2](https://arxiv.org/abs/hep-th/9511149),
Appendix A, PDF pages 33–34, was independently checked. The running assumptions
there are an even self-dual lattice and 2a,2b in that lattice. The paragraph
after (A.4) states oddness in the elliptic variable and vanishing at its origin
for odd characteristic pairing. This is relevant background under those
hypotheses; the project's substitution proof does not need them.

This does not mean that general odd-characteristic vanishing is a new theorem.
In coordinates the project's series is exactly the Riemann theta constant
with period matrix Omega=tau G and characteristic (a/2,b/2):

\[
\sum_{x\in\mathbb Z^n}\exp\bigl(\pi i(x+a/2)^t\Omega(x+a/2)
+2\pi i(x+a/2)^t(b/2)\bigr).
\]

Positive definiteness puts Omega in Siegel space. This identification is a
direct algebraic substitution, valid for every lattice in our setup.
NIST's [DLMF (21.2.5)](https://dlmf.nist.gov/21.2.E5) gives this definition,
and [DLMF (21.3.6)](https://dlmf.nist.gov/21.3.E6) gives general
half-characteristic parity. These reference formulas have now been read;
the general parity citation is no longer an unresolved task.

The characteristic shift law [DLMF (21.3.4)](https://dlmf.nist.gov/21.3.E4),
basis generator [DLMF (21.5.5)](https://dlmf.nist.gov/21.5.E5),
translation formula [DLMF (21.2.6)](https://dlmf.nist.gov/21.2.E6), and
characteristic transformation [DLMF (21.5.9)](https://dlmf.nist.gov/21.5.E9)
have also been checked. [docs/theta-transformations.md](theta-transformations.md)
derives the exact GL(n,Z) identity directly from the defining series, including
determinant minus one. Its combined prefactor is one for unreduced
characteristics; reducing the transformed characteristic on a lattice
stabilizer introduces precisely the project's epsilon sign. This is an
explicit specialization of standard theta transformations, not a novelty
claim. DLMF cites Igusa and Mumford for the general law; those books remain
unopened and are not represented as checked pages. reviews/theta-transformation-review.md
independently checked this comparison and found no correction necessary.

## Kane–Kim: background, with a qualified applicability statement

Kane and Kim, *Theta series of ternary quadratic lattice cosets*,
[arXiv:2211.03987v2](https://arxiv.org/html/2211.03987v2), sections 1–2 and
Theorem 1.2, was independently checked. It treats unsigned representation
series and their genus/spinor-genus decomposition. It does not directly prove
our symmetry converse or supply a coefficient cutoff for the signed series.

Calling that theory wholly inapplicable to our coset difference would be too
strong. For even p and b not in 2L*, put L0={x in L:<x,b> even}. Then a=2xi
belongs to L0, and 2t belongs to L0. Doubling the two shifted cosets in
docs/foundations-and-exact-enumeration.md section 8 gives a+2L0 and a+2t+2L0, both integral lattice cosets.
Their unsigned theta series at z=tau/8 reproduce the original coset difference
when the unsigned convention is exp(2 pi i z ||y||^2). Thus unsigned coset
results apply to the summands in their variable z under the corrected
conductor conventions below. Modularity alone does not imply vanishing or
the symmetry converse.

There is a base-lattice correction to reviews/modularity-follow-up.md section 3. A shift
nu outside 2L0 gives a conductor-two representation 2M+nu with M=L0.
If nu belongs to 2L0, the coset is the ordinary lattice 2L0; changing the
conductor to one while retaining M=L0 would instead give L0 and change
the series. Use M=2L0 and shift zero for this summand. For example, L=Z,
a=0 and b=1 give L0=2Z and the doubled cosets 4Z and 2+4Z; the first
cannot be replaced by 2Z.

reviews/theta-transformation-review.md now supplies the corrected common space. The level
and discriminant of 2L0 are 4N_L0 and 4^n det(L0). Proposition 2.3 puts
both summands in weight n/2 on Gamma0(16N_L0), with the same character on
the group. Its intersection with Gamma1(2) is redundant. The source's
definition and slash operator in section 2.2 and its Proposition 2.3 were
independently re-read during incorporation. [docs/coset-modularity.md](coset-modularity.md)
sections 1-3 give the convention translation and exact level proof.

Section 4 gives the reviewed derivation for tau/8: on
Gamma0(16N_L0) intersect Gamma^0(8), integral weight retains
the source character and odd rank acquires an extra chi_8 factor in the
half-integral slash convention. The group has cusp width eight at infinity
and does not satisfy the source's convention requiring translation by one.
reviews/coset-rescaling-review.md confirms the multiplier, including negative entries
and -I. The note now writes the valuation proof and the exact positive cusp
factor (d_c/8)^(-n/2). The main manuscript incorporates the resulting
holomorphy statement. No coefficient cutoff or cusp-form claim is made.

## Still unverified as prior results

Prior formulations of the general stabilizer sign character, the root-lattice
parity classification, the An shell interpretation and the uniform binary
classification remain unchecked. The books and additional leads listed in
reviews/initial-mathematical-review.md have not all been read. Ordinary isospectral lattices are not
counterexamples for the specified two cosets of one lattice. Neither general
novelty nor absence of relevant literature is asserted.
