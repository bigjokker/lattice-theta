# Transformation audit and corrected coset comparison

Prepared 2026-09-27 from the focused brief. The two earlier reports, the proofs,
the code, and the census packs were not changed. No lattice enumeration was
rerun. One truncated series comparison for a determinant-minus-one basis
change was computed directly and is marked below.

**Checked** means a formula was read on the page cited. **Deduction** means
an algebraic step from those formulas. **Unresolved** means a step that was
not established.

## 1. Transformation identity

No correction to `docs/theta-transformations.md` or to `MANUSCRIPT.md` §2 is
needed. The displayed identities agree with `docs/foundations-and-exact-enumeration.md` §§1–3 and with the
DLMF formulas that were read for this check.

### 1.1. What was read

Checked, NIST DLMF 21.2–21.5, version 1.2.8:

- (21.2.5) defines the characteristic series by
  `exp(2π i (1/2 (n+α)·Ω·(n+α) + (n+α)·(z+β)))`, which is the note's
  `π i` and `2π i` form.
- (21.2.6) writes that series as an exponential times the characteristic-free
  theta at `z + Ωα + β`.
- (21.3.4), in the corrected form after the 2012 erratum, is
  `θ[α+m1; β+m2](z|Ω) = exp(2π i α·m2) θ[α; β](z|Ω)` for integer vectors
  `m1, m2`. The factor depends on the lower shift.
- (21.3.6) is `θ[α; β](-z|Ω) = (-1)^{4α·β} θ[α; β](z|Ω)` for half-period
  characteristics.
- (21.5.5) is the basis generator
  `Γ = [A, 0; 0, (A^{-1})^T]`, with `A` any invertible integer matrix, and
  `θ(Az | AΩA^T) = θ(z|Ω)`. There is no extra factor.
- (21.5.9) transforms characteristics by
  `(Dα - Cβ + (1/2) diag(CD^T), -Bα + Aβ + (1/2) diag(AB^T))`
  and multiplies by an unspecified factor `κ(α,β,Γ)` times a determinant
  branch. (21.5.8) says that this determinant branch is a square root.

### 1.2. The unreduced identity

Deduction. For `U ∈ GL(n,ℤ)` and `y = Ux`, the summand of
`θ[U^{-1}α; U^t β](U^t z | U^t Ω U)` is the summand of
`θ[α; β](z|Ω)` at `y`. The map is a bijection of `ℤ^n` for both values of
`det U`. Identity (1) in the note follows, with prefactor 1, for unreduced
real characteristics.

The same identity is the composition of (21.5.5) with (21.2.6). Set the
DLMF matrix `A` equal to `U^t`. Then `(A^{-1})^T = U^{-1}`, so the block is
`diag(U^t, U^{-1})`. It sends `Ω` to `U^t Ω U` and `z` to `U^t z`. The
exponential prefactor in (21.2.6) is unchanged by this substitution, and
(21.5.5) contributes 1, including when `det U = -1`. The product of those
two formulas is identity (1).

Checked calculation. For `U = [[0,1],[1,0]]`, determinant `-1`, and a
positive-definite imaginary period matrix, a truncation of both series to
coordinates `-8,...,8` agreed to `2×10^{-17}`.

In (21.5.9) the same block has `B = C = 0`, so both diagonal corrections
vanish and the new characteristic is `(U^{-1}α, U^t β)`. The printed factor
`κ det(CΩ+D)` is not split into a square-root branch and a root of unity.
Identity (1) shows that the product of those pieces is 1 for this block.
That is what the note claims, and it is the right claim: (21.5.9) does not
by itself name the branch.

### 1.3. The stabilizer sign

Deduction. Let `T^t G T = G` and set `U = T^{-1}`. Then
`T^{-t} G T^{-1} = G`, so `Ω = τG` is fixed. Identity (1) at `z = 0` becomes

```
θ[Tα; T^{-t}β](0 | τG) = θ[α; β](0 | τG).
```

If `T` fixes both half-characteristic classes, `λ = (T-I)a/2` and
`μ = (T^{-t}-I)b/2` are integral. Formula (21.3.4) multiplies by
`exp(2π i α^t μ)`. Since `α = a/2`, this is `(-1)^{a^t μ}`, and
`a^t μ` is the integer whose class modulo 2 is `ε(T)` in `docs/foundations-and-exact-enumeration.md` §3.
A nontrivial sign forces vanishing.

Using `U = T` instead sends `β` by `T^t` and produces
`a^t (T^t - I)b/2`, which is `ε(T^{-1})`. On the stabilizer this equals
`ε(T)`, by the homomorphism in `docs/foundations-and-exact-enumeration.md` §2. The manuscript's sentence
that the reduced characteristic `(Tα, T^{-t}β)` contributes
`(-1)^{ε(T)}` is this calculation.

Odd vanishing is the same specialization. For `α = a/2` and `β = b/2`,
`4α·β = a·b = p`. Formula (21.3.6) makes the nullwert zero for every `Ω`
when `p` is odd. Restriction to `Ω = τG` is the project's odd-parity
implication. Integer shifts do not change that parity, by (21.3.4).

### 1.4. What this does not classify

The transformation is a standard specialization. It does not classify
lattices with the parity property, and it does not prove the symmetry
converse. `docs/theta-transformations.md` keeps that distinction, and
`MANUSCRIPT.md` §2 does not claim those classifications as consequences of
the DLMF formulas.

| Location | Statement | Severity | Reasoning | Correction |
|---|---|---|---|---|
| — | No incorrect identity was found in the transformation note or in manuscript §2. | — | The unreduced identity, both determinant signs, the block `diag(U^t, U^{-1})`, the vanishing diagonal corrections, the prefactor 1, and both the `U = T^{-1}` and `U = T` sign conventions were checked against (21.2.5), (21.2.6), (21.3.4), (21.3.6), (21.5.5) and (21.5.9). | None. |

## 2. Corrected coset comparison

The doubled-coset identity is kept. The conductor packaging in the earlier
follow-up is not. If a shift already lies in `2L_0`, the coset is the
lattice `2L_0`. Calling it a conductor-one coset of the base `L_0` replaces
that series by the theta series of `L_0`.

### 2.1. Checked conventions from Kane–Kim

Read again in arXiv:2211.03987v2, §2.1 and Proposition 2.3.

Their bilinear form satisfies `Q(x) = B(x,x)`. The Gram matrix
`A = (B(e_i, e_j))` is the project's Gram matrix when `B` is the inner
product. The discriminant `d_L` is `det A`. The level `N_L` is the smallest
positive integer `N` such that `N A^{-1}` is integral. A lattice is integral
when `B(L,L) ⊆ ℤ`, which forces `Q(aL+ν) ⊆ ℤ`.

The conductor of `L+v_0` is the smallest positive integer `a` such that
`a v_0 ∈ L`. Their coset notation `aL+ν` assumes `ν ∈ L` and that this
conductor, relative to the lattice `aL`, is exactly `a`. Conductor 1 means
the shift lies in the base lattice, so the coset is that lattice.

Proposition 2.3: for a coset of rank `k`, conductor `a`, level `N_L` and
discriminant `d_L`,

```
Θ_{aL+ν}(z) = Σ_{x ∈ aL+ν} exp(2π i z Q(x))
```

lies in weight `k/2` on `Γ_0(4 N_L a^2) ∩ Γ_1(a)`, with character
`χ_{4 d_L}` when `k` is odd and `χ_{(-1)^{k/2} 4 d_L}` when `k` is even.
The character and the level are those of the base lattice `L`, not of a
sublattice chosen later. The proof evaluates the character at the odd lower
row of a matrix in `Γ_0(4N_L a^2)`.

Theorem 1.2, the genus and spinor-genus splitting, remains a rank-3
statement and is not used below.

### 2.2. The series identity

Deduction. Let `b = 2δ ∉ 2L^*`, `L_0 = {x ∈ L : ⟨x,b⟩ even}`, and
`a = 2ξ ∈ L_0`. Choose `t ∈ L` with `⟨t,b⟩` odd. Then `2t ∈ L_0` but
`2t ∉ 2L_0`, so `a+2L_0` and `a+2t+2L_0` are distinct cosets of `2L_0`.

For `y = 2x+a`,

```
||x+ξ||^2 = ||y||^2 / 4,
π i τ ||x+ξ||^2 = 2π i (τ/8) ||y||^2.
```

Therefore, with `Φ_C(z) = Σ_{y ∈ C} exp(2π i z ||y||^2)`,

```
Θ_L[ξ,δ](τ) = i^p ( Φ_{a+2L_0}(τ/8) - Φ_{a+2t+2L_0}(τ/8) ).
```

Each `Φ_C(z)` is a Kane–Kim series for `Q = ||·||^2`. This fixes `z = τ/8`
as an equality of series. It does not put `Φ(τ/8)` into a modular form in
the variable `τ`.

### 2.3. Which base to use

Deduction. Let `n = rank L`.

- If the shift `ν` lies in `2L_0`, then `2L_0+ν = 2L_0`. Package it as
  conductor 1 on the base `M = 2L_0`, with shift 0. The series is the
  ordinary theta series of the lattice `2L_0`.
- If the shift does not lie in `2L_0`, package it as conductor 2 on the
  base `M = L_0`: the set is `2L_0+ν` with `ν ∈ L_0`.

In the example `L = ℤ`, `a = 0`, `b = 1`, one has `L_0 = 2ℤ` and the doubled
cosets are `4ℤ` and `2+4ℤ`. The first is conductor 1 on `M = 4ℤ`. Conductor
1 on `M = 2ℤ` is the different series `Σ_{m ∈ 2ℤ} q^{m^2}`. The second is
conductor 2 on `M = 2ℤ`, since the smallest positive `c` with `c·2 ∈ 4ℤ` is
`c = 2`.

Both cosets nontrivial means both shifts lie outside `2L_0`, so both use
`M = L_0` and conductor 2. One trivial means one uses `M = 2L_0` and
conductor 1, and the other uses `M = L_0` and conductor 2.

### 2.4. Level and character from `L_0` to `2L_0`

Deduction. Let `G` be a Gram matrix of `L_0`, with discriminant `d` and
level `N`. On the doubled basis the Gram matrix of `2L_0` is `4G`, so its
determinant is `4^n d` and its inverse is `G^{-1}/4`.

Let `H = N G^{-1}`, an integer matrix, and let `g` be the gcd of its
entries. No prime dividing `g` can divide `N`: otherwise `N/p` would clear
the denominators of `G^{-1}`. But `det H = N^n / d` is divisible by `g^n`,
so `g^n` divides `N^n`. The only possibility is `g = 1`. Therefore the level
of `2L_0` is exactly `4N`, and its discriminant is `4^n d`.

Proposition 2.3 then gives:

- conductor 2 on `L_0`: weight `n/2`, group
  `Γ_0(4 N · 4) ∩ Γ_1(2) = Γ_0(16N) ∩ Γ_1(2)`, character `χ` of `L_0`;
- conductor 1 on `2L_0`: weight `n/2`, group
  `Γ_0(4·(4N)) = Γ_0(16N)`, character `χ` of `2L_0`.

Here `χ` means `χ_{4d}` for odd rank and `χ_{(-1)^{n/2} 4d}` for even rank,
with `d` the discriminant of the base being used. The two labels differ by
the square factor `4^n`. On odd integers, which are the only arguments at
which these characters are evaluated for matrices in `Γ_0(16N)`, the
Kronecker symbol of that square factor is 1. The two characters agree on
the group.

Weight `n/2` is half-integral for odd `n` and integral for even `n`.
Proposition 2.3 states both cases. For odd `n` the level `16N` is divisible
by 4, as the half-integral slash operator requires.

### 2.5. The common space in the variable `z`

Deduction. In all three configurations the two unsigned series lie in

```
M_{n/2}( Γ_0(16 N_{L_0}) ∩ Γ_1(2), χ_{L_0} ),
```

with `χ_{L_0}` the character of Proposition 2.3 for `L_0`. Their difference
lies in that space. This is a statement about `Φ(z)`, not about `Φ(τ/8)`.

Unresolved. The slash operator for the substitution `z = τ/8` was not
computed. A form of level `16 N_{L_0}` in Kane–Kim's variable need not have
`f(τ/8)` of any particular level in `τ` until that cocycle is calculated.
No coefficient cutoff is asserted. Cusps, the Fourier normalization after
the substitution, and the rank-3 spinor decomposition are not inputs to this
comparison.

## 3. Remaining gaps

- The modular form in the variable `τ`, as opposed to the series identity at
  `z = τ/8`, is unresolved.
- Nebe, Rains, and Sloane, *Self-Dual Codes and Invariant Theory*, Springer,
  2006, Chapter 9, pages 249–284, remains unread. The 2008 lecture PDF in
  `sources/clifftype.pdf` does not replace it.
- No checked antecedent was found in this pass for the root-lattice list,
  the `A_n` shell evaluation, or the rank-two indecomposability criterion.
  That is not a novelty claim.
- Igusa's book and Mumford's book were not opened. The DLMF formulas used
  here cite those books and were themselves read; the citations in the
  transformation note do not present the books as checked pages.
