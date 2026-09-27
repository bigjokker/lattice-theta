# Rescaling and cusp review

Prepared 2026-09-27 from `docs/coset-modularity.md`. Proofs, the manuscript, code, and census files were not changed. The completed `GL(n,ℤ)` audit is not repeated.

**Checked** means a source formula or a finite arithmetic check. **Deduction** means an algebraic step from those formulas. **Unresolved** means a step that was not taken.

The truncated series in `reports/coset-rescaling-checks.json` were read and not recomputed. Their recorded relative errors are about `4.6×10^{-13}` and `5.5×10^{-13}`. The matrix arithmetic and the Kronecker value `(8/3)=-1` were recomputed directly.

## 1. Issue table

| Location | Statement | Severity | Reasoning | Correction |
|---|---|---|---|---|
| §2, the sentence “so `g^n` divides `N^n`” | From `g^n` dividing `det H = N^n/d`, the note concludes that `g^n` divides `N^n`. | missing justification | `g^n` divides `N^n/d` does not by itself put `g^n` into `N^n`. The conclusion `g=1` is nevertheless true. If a prime `p` divides every entry of `H=N G^{-1}`, then `p^n` divides `det H`. The same prime cannot divide `N`, or `N/p` would clear the denominators of `G^{-1}`. Then the `p`-adic valuation gives `0 = v_p(N^n) ≥ n + v_p(d) ≥ n`, which is impossible. | Replace that sentence by this valuation. The next sentence, that `s(4G)^{-1}` is integral exactly when `4N` divides `s`, is correct once `g=1`. The level of `4G` is exactly `4N`. |
| §4, the cusp paragraph beginning “Given a cusp frame” | The expansion of `g` is the expansion of `F` up to a nonzero constant weight factor, after a positive rational dilation. | missing justification | The constant is real and nonzero only after the frame is normalized. Let `R=[[1,0],[0,8]]`, so `R(τ)=τ/8`. For `σ'^{-1} R σ = [[a,b],[0,d]]`, one has `ad=8`. Replacing `σ'` by `-σ'` flips both diagonal entries and preserves the cusp. Choose the sign so that `d>0`; then `a>0` as well. The cocycle gives `j(σ,τ)^{-k} g(σ τ) = (d/8)^{-k} φ(Uτ)`, where `φ` is the holomorphic cusp expansion of `F` and `Uτ=(aτ+b)/d`. For `d>0` the factor `(d/8)^{-k}` is the positive-real power. It introduces no pole. Nonnegative exponents stay nonnegative under `τ ↦ (a/d)τ + b/d` with `a/d>0`. | Insert this normalization and the factor `(d/8)^{-k}`. The holomorphy claim, and the refusal to claim vanishing at the cusps, then stand. |

No incorrect transformation law was found in (3) or (4). `Γ_1(2)` is redundant in §3, as claimed.

## 2. The rescaling proposition

Checked source formulas. Kane–Kim, arXiv:2211.03987v2, §2.1 and Proposition 2.3, give the weight, group, and character of each unsigned coset series in the variable `z`. Their §2.2 defines the half-integral slash

```
(f|_k γ)(z) = (c/d) ε_d^{2k} (cz+d)^{-k} f(γz),
```

with `ε_d = 1` for `d ≡ 1 (mod 4)` and `ε_d = i` for `d ≡ 3 (mod 4)`, and it requires the group to contain translation by one. The Kronecker symbol supplies `χ_D`, including at negative arguments.

Deduction. Let `N = N_{L_0}`, `d = det G_{L_0}`, `k = n/2`, and let `χ` be `χ_{4d}` for odd `n` and `χ_{(-1)^{n/2} 4d}` for even `n`. Section 2 shows that the level of `2L_0` is `4N` and its discriminant is `4^n d`. On odd integers the square factor `4^n` does not change `χ`. Every matrix in `Γ_0(16N)` has odd diagonal entries, so it already lies in `Γ_1(2)`. Thus

```
F ∈ M_k(Γ_0(16N), χ)
```

in Kane–Kim’s normalization `q = exp(2π i z)`. This is the incorporated statement (2). The two cosets are distinct, so at most one is the lattice `2L_0` itself.

Let `g(τ) = F(τ/8)` and

```
H = Γ_0(16N) ∩ Γ^0(8),
```

where `Γ^0(8)` consists of the matrices in `SL_2(ℤ)` whose upper-right entry is divisible by 8. Then `Γ(16N) ⊂ H`. The translations in `H` are exactly the translations by multiples of 8, so the cusp width at infinity is 8. In particular `[[1,1],[0,1]]` is not in `H`, and Kane–Kim’s §2.2 does not apply to `H` until the definition is widened to arbitrary cusp widths.

For `γ = [[r,s],[u,v]] ∈ H`, the conjugate

```
γ' = [[r, s/8], [8u, v]] = R γ R^{-1},  R = [[1,0],[0,8]],
```

is an integral matrix of determinant 1 in `Γ_0(16N)`. Both `γ` and `γ'` lie in `Γ_0(4)`. Also `γ'(τ/8) = γ(τ)/8` and `8u·(τ/8)+v = uτ+v`. The lower-right entry `v` is odd.

For even `n`, integer powers are single-valued on `ℂ*`, and

```
g(γτ) = χ(v) (uτ+v)^k g(τ).
```

Equivalently, with the integral slash `(cτ+d)^{-k}`, one has `g|_k γ = χ(v) g`. Negative `v` are included: `(-I)τ = τ`, `χ(-1) = (-1)^{n/2}`, and `χ(-1)(-1)^k = 1`, so `-I` acts consistently.

For odd `n`, the half-integral slash uses one common value of `ε_v` and one common branch of `(uτ+v)^{-k}`, because both matrices have the same lower-right entry and the same linear form `uτ+v`. Multiplicativity for odd `v` gives `(8u/v) = (8/v)(u/v)`, and `(8/v)^2 = 1`, so

```
g|_k γ = χ(v) χ_8(v) g,  χ_8(v) = (8/v).
```

The principal argument `Arg ∈ (-π, π]` gives `Arg(-1) = π`. For that branch, `-I` acts by `+1`, which equals `χ(-1)χ_8(-1)` because both `4d` and `8` are positive. No further restriction on `H` is required. Both `χ` and `χ χ_8` are Dirichlet characters modulo a divisor of `16N`, hence are well-defined on the lower-right entries of `H`.

The Fourier expansion at infinity is

```
g(τ) = Σ_{m≥0} A_m exp(2π i m τ / 8).
```

It is an expansion in the width-eight parameter. It need not be an integral series in `exp(2π i τ)`.

At every cusp, after the frame normalization above, `g` has nonnegative Fourier exponents in a local parameter and at most polynomial growth. It is holomorphic at the cusps. The constant term need not vanish, so this is not cuspidality. No smallest group and no coefficient cutoff are claimed.

Checked arithmetic. For `γ = [[171,8],[64,3]]`, the determinant is 1, the upper-right entry is divisible by 8, and `(8/3) = -1`. For the standard lattice `ℤ` with `b = 1`, one has `L_0 = 2ℤ`, Gram `(4)`, level `4`, and `16N = 64`, which divides the lower-left entry `64`. The same level pattern holds for the standard `ℤ^2` example with `b = (1,0)`. These are the matrices recorded in the check file. The series truncation itself was not rerun.

## 3. What remains unresolved

- Nothing in this rescaling produces a Sturm bound or a minimum level.
- The genus and spinor-genus splitting in Kane–Kim remains a rank-3 statement and was not used.
- Whether every form obtained this way is already modular for a subgroup of smaller level was not investigated.
