# Standard theta transformations and the lattice sign

Checked 2026-09-27 following GROK-FOLLOWUP.md. This note resolves the general
parity citation and the GL(n,Z) prefactor question with accessible reference
formulas and a direct series proof. It does not claim priority for the lattice
criterion or access to the unopened Igusa/Mumford books.
The subsequent audit in [reviews/theta-transformation-review.md](../reviews/theta-transformation-review.md)
confirmed the identities without correction; its truncated numerical check
is corroboration, not the proof.

## Reference normalization and parity

NIST's [DLMF (21.2.5)](https://dlmf.nist.gov/21.2.E5) defines

\[
\theta\!\begin{bmatrix}\alpha\\\beta\end{bmatrix}(z\mid\Omega)
=\sum_{x\in\mathbb Z^n}
\exp\bigl(\pi i(x+\alpha)^t\Omega(x+\alpha)
+2\pi i(x+\alpha)^t(z+\beta)\bigr).
\]

With z=0, Omega=tau G, alpha=a/2 and beta=b/2 this is exactly the project's
series; b is in dual coordinates. [DLMF (21.3.6)](https://dlmf.nist.gov/21.3.E6)
states parity for half characteristics with sign (-1)^(4 alpha dot beta).
Thus odd p=a dot b forces the nullwert to vanish for every Omega. Restriction
to tau G proves the project's odd-parity implication. There is no narrower
even/unimodular restriction in this reference formula.

## Exact identity under an integral basis change

For U in GL(n,Z), the bijection y=Ux in the defining sum proves

\[
\theta\!\begin{bmatrix}U^{-1}\alpha\\U^t\beta\end{bmatrix}
(U^t z\mid U^t\Omega U)
=\theta\!\begin{bmatrix}\alpha\\\beta\end{bmatrix}(z\mid\Omega).
\tag{1}
\]

The quadratic term becomes (y+alpha)^t Omega(y+alpha), and the linear term
becomes (y+alpha)^t(z+beta). No factor is introduced by a bijection of Z^n.
This proves (1) equally for det U=1 and det U=-1 and for unreduced real
characteristics, not just their binary representatives.

The symplectic block is diag(U^t,U^{-1}). It acts on Omega by U^t Omega U,
on z by U^t z, and on the characteristic by (U^{-1}alpha,U^t beta).
Its off-diagonal blocks vanish, so the diagonal characteristic corrections
in [DLMF (21.5.9)](https://dlmf.nist.gov/21.5.E9) vanish. The basis generator
[DLMF (21.5.5)](https://dlmf.nist.gov/21.5.E5), together with the characteristic
translation formula [DLMF (21.2.6)](https://dlmf.nist.gov/21.2.E6), agrees
with (1): its exponential is unchanged by the basis substitution.
Consequently the combined automorphy prefactor
for this block and these unreduced conventions is 1. This does not assign
separate branch values to the square root and root-of-unity factors of a
general symplectic transformation.

## Reduction to the lattice stabilizer sign

Let T be a lattice isometry, take U=T^{-1} in (1), and set z=0. Since
T^{-t}G T^{-1}=G, the period matrix remains tau G and

\[
\theta\!\begin{bmatrix}T\alpha\\T^{-t}\beta\end{bmatrix}(0\mid\tau G)
=\theta\!\begin{bmatrix}\alpha\\\beta\end{bmatrix}(0\mid\tau G).
\]

If T fixes both half-characteristic classes, write
lambda=(T-I)a/2 and mu=(T^{-t}-I)b/2 in Z^n.
The characteristic shift law [DLMF (21.3.4)](https://dlmf.nist.gov/21.3.E4) gives

\[
\theta\!\begin{bmatrix}\alpha+\lambda\\\beta+\mu\end{bmatrix}
=e^{2\pi i\alpha^t\mu}
\theta\!\begin{bmatrix}\alpha\\\beta\end{bmatrix}
=(-1)^{a^t\mu}\theta\!\begin{bmatrix}\alpha\\\beta\end{bmatrix}.
\]

The exponent is precisely epsilon(T) in PROOFS.md. A nontrivial sign forces
vanishing. Using U=T instead gives the inverse convention with T^t on b,
namely epsilon(T^{-1}), equal to epsilon(T) on the stabilizer.
The homomorphism and representative-independence proofs remain those of
docs/foundations-and-exact-enumeration.md section 2. The connection to standard transformations is now explicit;
whether the particular classification statements have prior formulations is
a separate literature question.
