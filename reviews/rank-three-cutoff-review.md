# Coefficient cutoff for the rank-three signed coset difference

Prepared 2026-09-27 from `Grok-rank3.md`. The manuscript, earlier reviews, census packs, and saved reports were not changed. The rank-three determinant-`≤24` census is already closed by complete shells; this note is a general identity criterion for the same series. It does not enlarge that census and it does not review the manuscript again.

**Checked** means a source statement read in this pass, or an exact integer or shell computation run in this pass. **Deduction** means a step proved from those statements. The `GL(n,ℤ)` identities in `docs/theta-transformations.md` are not used.

## 1. Theorem

Let `G` be a positive definite integral `3×3` Gram matrix and let `b` be a nonzero column in `{0,1}^3`. Put `L_0 = {x ∈ ℤ^3 : x^t b even}`, choose a `ℤ`-basis matrix `K` of `L_0`, and write `G_0 = K^t G K`, `d_0 = det G_0`, and `N_0` for the least positive integer such that `N_0 G_0^{-1}` is integral. Let

```
M = 16 N_0,
μ(M) = [SL_2(ℤ) : Γ_0(M)] = ∏_{p^e || M} p^{e-1}(p+1).
```

Then `8` divides `μ(M)`, and the positive integer

```
B(N_0) = μ(16 N_0)/8
       = 2^{v_2(N_0)} · 3 · ∏_{p^a || N_0, p odd} p^{a-1}(p+1)
```

depends only on `N_0`. Let `c_N` be the norm4 coefficient in the brief, the coefficient of `q_z^N` in `F`, with `q_z = exp(2π i z)`.

**Theorem.** If `c_N = 0` for every integer `N` with `0 ≤ N ≤ B(N_0)`, then `F` is identically zero, and therefore `Θ(τ) = i^p F(τ/8)` is identically zero.

The test includes the constant term. One complete coefficient `c_N ≠ 0` already proves that `F` is not identically zero; the bound is used only in the vanishing direction. The group `Γ_0(16 N_0)` is a sufficient level, carried forward from `docs/coset-modularity.md`. This theorem does not assert that the level is minimal, that `F` is a cusp form, or that `B(N_0)` is the largest integer with the same implication.

## 2. The modular input, including the trivial coset

Checked source. Kane–Kim, arXiv:2211.03987v2, §2.1, define the discriminant as `det G` and the level as the least positive integer `N` for which `N G^{-1}` is integral. Their §2.2 slash operator on `Γ_0(4)`, for weight `κ ∈ ℤ + 1/2`, is

```
(f |κ γ)(z) = (c/d) ε_d^{2κ} (cz+d)^{-κ} f(γz),
ε_d = 1 if d ≡ 1 (mod 4),   ε_d = i if d ≡ 3 (mod 4).
```

A holomorphic modular form on a congruence subgroup `Γ ⊆ Γ_0(4)` which contains the translation `T = [[1,1],[0,1]]`, with character `χ`, satisfies `f |κ γ = χ(d) f`, is holomorphic on the upper half-plane, and grows at most polynomially in `y` as `z = x+iy` approaches `ℚ ∪ {i∞}`. Proposition 2.3 places the theta series of an integral coset of conductor `a` and rank `k` in weight `k/2` on `Γ_0(4 N_L a^2) ∩ Γ_1(a)`, with character `χ_{4 d_L}` when `k` is odd. The proof’s identity (2.9) is the same slash relation before the conductor forces the shift `ν` to be fixed.

Checked project note. `docs/coset-modularity.md` §§1–3, re-read in this pass, package the two doubled cosets of the signed difference as follows. For `b ≠ 0` the kernel `L_0` has index two, `d_0 = 4 det G`, and the two doubled cosets are distinct. A shift already in `2L_0` is the lattice `2L_0` itself, with base `2L_0` and conductor one. A shift outside `2L_0` has base `L_0` and conductor two. Replacing the trivial coset by the conductor-one coset of `L_0` changes the series. Both summands, and therefore `F`, lie in

```
M_{3/2}(Γ_0(16 N_0), χ_{4 d_0}),   χ_D(v) = (D/v).
```

Here `4 d_0 = 16 det G`. The factor `Γ_1(2)` is redundant: every matrix of `Γ_0(16 N_0)` has odd diagonal entries. The Fourier index of `F` in `q_z` is exactly norm4, and the expansion contains `N = 0`. The cusp width of `Γ_0(16 N_0)` at infinity is `1`, because `T` lies in the group. This is the variable in which the cutoff is stated. The rescaled form `g(τ) = F(τ/8)` is not needed.

The same section’s valuation shows that if `H = N_0 G_0^{-1}`, then `d_0 det H = N_0^3`. Every prime divisor of `d_0` divides `N_0`. For `γ = [[*,*],[c,d]] ∈ Γ_0(M)`, one has `gcd(d, M) = 1`. Thus `gcd(d, 4 d_0) = 1`, the value `χ_{4 d_0}(d)` is `±1`, and its fourth power is `1`. The symbol `(4 d_0 / -1)` equals `1` because `4 d_0 > 0`, so the character is even. Both facts were checked on the sample discriminants below with an exact Kronecker implementation, including `(0/±1) = 1`.

## 3. Fourth power, character, and central sign

Deduction. Let `κ = 3/2` and let `f` stand for `F`. The modularity relation `f |κ γ = χ(d) f` rearranges to

```
f(γz) = χ(d) (c/d)^{-1} ε_d^{-3} (cz+d)^{3/2} f(z),
```

where `(cz+d)^{3/2}` is the branch fixed by Kane–Kim’s slash operator. For `γ ∈ Γ_0(M)` the lower-right entry `d` is odd and `gcd(c,d) = 1`, so `(c/d) = ±1` and `(c/d)^{-1} = (c/d)`. Raise the identity to the fourth power. Then `χ(d)^4 = 1` and `(c/d)^{-4} = 1`. Also `ε_d^4 = 1`, so `ε_d^{-12} = 1`. The branch satisfies

```
[(cz+d)^{3/2}]^4 = (cz+d)^6,
```

the single-valued integer power: if `Log` is the logarithm chosen by the slash operator, then `exp(6 Log(cz+d)) = (cz+d)^6`. Therefore

```
F(γz)^4 = (cz+d)^6 F(z)^4,   γ ∈ Γ_0(M).
```

Equivalently, with the integral-weight slash `(φ |_6 γ)(z) = (cz+d)^{-6} φ(γz)`,

```
F^4 |_6 γ = F^4.
```

The character of `F^4` on `Γ_0(M)` is trivial. The weight `6` is even. At `-I`, the factor `(-1)^6` equals `1`, which agrees with the trivial character at `d = -1` and with `(-I)z = z`.

The square `F^2` is not the form used here. Its weight would be `3`. Brunault’s §0.1 records that an odd-weight form for a congruence subgroup containing `-I` is zero. Since `-I ∈ Γ_0(M)`, there is no nonzero element of `M_3(Γ_0(M))` with trivial character. The even positive discriminant character of `F` does not supply the sign `χ(-1) = -1` that an odd-weight nebentypus would need. The fourth power is the integral-weight form whose character and central sign match.

Holomorphy. Each coset theta series lies in Kane–Kim’s space, so it is holomorphic on the upper half-plane and of polynomial growth at the cusps. The same is true of their difference `F`, by the vector-space property recorded in `docs/coset-modularity.md` (2). Hence `F^4` is holomorphic on the upper half-plane. Polynomial growth passes to the fourth power. At a cusp `σ(∞)`, the function

```
φ(z) = (c_σ z + d_σ)^{-6} F(σz)^4
```

is invariant under a positive translation, the width of that cusp, because `F^4` has trivial character and even weight. Its Fourier series in the corresponding `q^{1/w}` therefore has only finitely many negative powers, and polynomial growth as the imaginary part tends to infinity forces those coefficients to vanish. Thus `F^4` is holomorphic at every cusp. It lies in `M_6(Γ_0(M))` in the sense required by the integral-weight theorem below. No coefficient is asserted to vanish at any cusp: holomorphy allows a nonzero constant term in the local parameter.

## 4. The integral-weight vanishing theorem

Checked source. François Brunault, *Sturm bounds for general congruence subgroups*, 28 May 2021, Theorem 1, read from the author’s PDF. Let `Γ` be a congruence subgroup of `SL_2(ℤ)` and let `m` be the index of `±Γ` in `SL_2(ℤ)`. If `f ∈ M_k(Γ)` and

```
ord_∞(f) > k m / 12,
```

equivalently `a_n(f) = 0` for all integers `n` with `0 ≤ n ≤ ⌊k m / 12⌋`, then `f = 0`. The order is the least exponent in the expansion `∑ a_n q^{n/w}`, where `w` is the width of the cusp infinity. The note cites Jacob Sturm, Lecture Notes in Mathematics 1240 (1987), 275–280, for the classical case. That Springer chapter was not opened in this pass. For `Γ_0(M)` one has `-I ∈ Γ`, so `±Γ = Γ` and Brunault’s `m` equals `μ(M)`. The width at infinity is `1`. For even weight the expansion is an ordinary integer power series in `q`. In this width-one situation the numerical threshold is the classical one; the note’s improvement for width greater than one is not what produces the factor `1/8` below.

Brunault’s §0.1 also records the central-sign obstruction used in §3: odd weight and `-I ∈ Γ` force the zero form.

Deduction. Apply Theorem 1 to `F^4`, in weight `6`. The vanishing threshold is

```
⌊6 μ(M) / 12⌋ = ⌊μ(M)/2⌋ = μ(M)/2,
```

the last equality because `8` divides `μ(M)`, as proved in §5. Suppose `c_N(F) = 0` for `0 ≤ N ≤ B` with `B = μ(M)/8`. In the product of four series, a coefficient of `q^n` in `F^4` with `n ≤ 4B+3` has a factor `c_j` with `j ≤ B`, because four indices each at least `B+1` would sum to at least `4B+4`. Those coefficients vanish. In particular they vanish for `0 ≤ n ≤ μ(M)/2`, since

```
4B+3 = μ(M)/2 + 3 ≥ μ(M)/2.
```

Brunault’s theorem gives `F^4 = 0`.

Characteristic zero. Let `m` be the least index with `c_m ≠ 0`, or `m = +∞` if there is none. The coefficient of `q^{4m}` in `F^4` is `c_m^4`. In `ℂ`, a nonzero complex number has nonzero fourth power, so

```
ord_∞(F^4) = 4 ord_∞(F)
```

whenever `F ≠ 0`. The same conclusion `F = 0` follows because the ring of formal power series `ℂ[[q]]` is an integral domain. This is an identity of complex Fourier series. It is not a congruence statement: modulo a prime, a nonzero coefficient can have vanishing fourth power, and Brunault’s theorem as used here is the complex form of the bound.

The resulting test on `F` is exactly `0 ≤ N ≤ B(N_0)`. If `μ/8` were not an integer, the same argument would use `⌊μ/8⌋`; divisibility makes the floor unnecessary. The inequality is strict on the order of `F^4`: that order is at least `4(B+1) = μ/2 + 4 > μ/2`. Vanishing only through `B-1` would give order at least `μ/2` for `F^4`, which the valence threshold still permits. The inclusive bound `B` is the one proved here.

## 5. Exact index arithmetic

Deduction, integer factorization only. For `M = ∏ p^e`,

```
μ(M) = ∏_p p^{e-1}(p+1).
```

This is the standard coset count `[SL_2(ℤ) : Γ_0(M)]`. Now `M = 16 N_0`, so the power of `2` dividing `M` is at least `4`. The factor `2^{e-1}(2+1)` in `μ(M)` is therefore divisible by `2^3 = 8`. Hence `B(N_0) = μ(M)/8` is an integer.

Writing `N_0 = 2^{v} t` with `t` odd gives `M = 2^{v+4} t` and

```
μ(M) = 2^{v+3} · 3 · ∏_{p^a || t} p^{a-1}(p+1),
B(N_0) = 2^{v} · 3 · ∏_{p^a || t} p^{a-1}(p+1).
```

The empty product, when `N_0` is a power of `2`, equals `1`, and then `B(N_0) = 3 N_0`.

Pseudocode, using exact integer and rational arithmetic:

```
index_gamma0(M):
    n, result, p = M, 1, 2
    while p*p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0: n //= p; e += 1
            result *= p**(e-1) * (p+1)
        p += 1 if p == 2 else 2
    if n > 1: result *= n + 1
    return result

level(G0):
    return the least positive N such that N * G0^{-1} has integer entries
    (the lcm of the denominators of G0^{-1} in lowest terms)

cutoff(N0):
    M = 16 * N0
    mu = index_gamma0(M)
    assert mu % 8 == 0
    return mu // 8
```

The kernel basis for a binary column `b` with a coordinate `b_i` odd may be taken as the columns `2 e_i` and `e_j + (b_j mod 2) e_i` for `j ≠ i`. Its determinant is `±2`, and `G_0 = K^t G K` satisfies `det G_0 = 4 det G`.

## 6. Worked levels and shells

Checked computation, run in this pass with exact rational LDL enumeration from `exact_theta.complete_shells` and an independent Kronecker implementation. No census file was rewritten. The three requested levels match the engine controls in the brief.

For `b = (1,0,0)` the basis is `K = diag(2,1,1)` in every row below.

| Gram `G` | `det G` | `G_0` | `d_0` | `G_0^{-1}` denominators | `N_0` | `M` | `μ` | `B` | `4 d_0` |
|---|---:|---|---:|---|---:|---:|---:|---:|---:|
| `I_3` | 1 | `diag(4,1,1)` | 4 | `4,1,1` | 4 | 64 | 96 | 12 | 16 |
| `A_3` Cartan | 4 | `[[8,-2,0],[-2,2,-1],[0,-1,2]]` | 16 | lcm `16` | 16 | 256 | 384 | 48 | 64 |
| `2 I_3` | 8 | `diag(8,2,2)` | 32 | `8,2,2` | 8 | 128 | 192 | 24 | 128 |

The `A_3` inverse is

```
[[3/16, 1/4, 1/8], [1/4, 1, 1/2], [1/8, 1/2, 3/4]].
```

Its level is `16`. Substituting `det G = 4` for the level would be wrong. For `2 I_3` the character discriminant is `128`, and `(128/3) = (2/3) = -1`, so this sample has nontrivial nebentypus. The fourth-power character is still trivial. For `I_3`, `(16/v) = 1` for the odd `v` tested.

Shells, all complete within the listed bound, with the signed coefficient after the global phase has been removed as in `code/exact_theta.py`:

| Gram | `a` | `b` | First nonzero norm4 | Coefficient | A cancelled earlier shell |
|---|---|---|---:|---:|---|
| `I_3` | `(0,0,0)` | `(1,0,0)` | 0 | 1 | none through 12 |
| `I_3` | `(0,1,0)` | `(1,0,0)` | 1 | 2 | norm4 `5`, 8 vectors, coefficient 0 |
| `I_3` | `(0,1,1)` | `(1,0,0)` | 2 | 4 | none through 12 |
| `A_3` | `(0,1,0)` | `(1,0,0)` | 2 | 2 | norm4 `6`, coefficient 0 |
| `2 I_3` | `(0,1,0)` | `(1,0,0)` | 2 | 2 | norm4 `10`, coefficient 0 |
| `[[2,-1,-1],[-1,4,-1],[-1,-1,4]]` | `(0,1,1)` | `(1,0,0)` | 10 | 2 | norm4 `6`, 4 vectors, coefficient 0 |

For the last Gram, `det G = 20`, `d_0 = 80`, `N_0 = 80`, `M = 1280`, `μ = 2304`, and `B = 288`. The stored pack records the same `K`, `G_0`, `d_0`, `N_0`, unsigned level `1280`, and character discriminant `320`, with `cutoff` null. The coefficient `2` at norm4 `10` shows that the series is nonzero, while the zero coefficient at norm4 `6` shows that a cancelled first shell is not an identity. Both numbers sit below `B = 288`. The same pattern holds for the `I_3`, `A_3`, and `2 I_3` rows: each displayed nonzero coefficient occurs at an index `≤ B(N_0)`.

The pack file `data/census-rank3/det24-bound32.json` has SHA-256 `1083147c76d6b66aeed3491d35bcb1f83ebed65c9cc0c630e95eaac90b32816a`, matching `RANK3-CENSUS.md`. Its saved summary records 303 generated matrices, 120 classes, 4,320 even pairs, 723 symmetry zeros, 3,597 nonzero shells, and 0 unresolved pairs. Those are stored fields. This pass did not replay the census.

Across that pack, `N_0` runs from 2 to 96. The formula gives `B(2) = 6` and `B(96) = 384`. The census note records that the largest norm4 actually used by a nonzero certificate is 26. A bound of 384 is a sufficient identity test; the completed census did not need it.

## 7. Limitations

The criterion certifies identical vanishing of `F` from finitely many exact coefficients. It does not certify that the stabilizer character is trivial, and a symmetry-converse counterexample still requires both ingredients named in the brief. Supplied symmetry witnesses remain what they were: proofs of vanishing by an epsilon-one automorphism, not applications of this bound.

The level `16 N_0` was not proved to be the minimal modular level. A future proof that `F` lies on a proper subgroup, or on `Γ_0` of a proper divisor of `16 N_0`, could shrink `μ` and therefore `B`. No such reduction is claimed, and none is inferred from small coefficients. Elliptic-point corrections that improve valence bounds for special levels are not used.

The argument uses holomorphy at every cusp and refuses cuspidality. Weight `3/2` enters through the fourth power and the denominator `8`. The same pattern applies to no other rank in this note. Rank four, determinants above 24, CM provenance, and a further manuscript review stay outside the assignment.

Brunault’s theorem is used as a complex identity theorem. A mod-`m` Sturm congruence is a different statement and is not claimed for `F` or for `F^4`.

## 8. Certificate checklist

A machine certificate that quotes this cutoff should contain:

1. `G`, `a`, and `b`, with `b ≠ 0` and `p = a^t b` even.
2. The basis `K` of `{x : x^t b even}`, the matrix `G_0 = K^t G K`, the integer `d_0 = det G_0`, and the check `d_0 = 4 det G`.
3. `N_0`, together with the exact matrix `N_0 G_0^{-1}` and the fact that no smaller positive integer clears the denominators.
4. `M = 16 N_0`, the factored index `μ(M)`, and `B = μ(M)/8`.
5. The character label `χ_{4 d_0}` and the weight label `3/2`. The vanishing proof is the identity for `F^4` in `M_6(Γ_0(M))`, not a congruence and not a cusp-form hypothesis.
6. Either one complete shell with `c_N ≠ 0`, or complete shells for every norm4 from `0` through `B` with every signed coefficient `0`, and an explicit statement that the enumerator did not hit a resource limit. A zero coefficient at the first occupied shell is not item 6.
7. If the certificate is also offered as a symmetry-converse counterexample, the full characteristic stabilizer and the value epsilon `= 0` on every element. This cutoff does not compute that group.

## 9. What this pass ran and read

Ran: the exact index, inverse, Kronecker, and level arithmetic for `I_3`, the `A_3` Cartan matrix, `2 I_3`, and the determinant-20 Gram in `RANK3-CENSUS.md`; `exact_theta.complete_shells` on the ten shells listed in §6; SHA-256 of `data/census-rank3/det24-bound32.json`; a read of that pack’s summary and of the three stored `modular_inputs` records for `I_3`, `2 I_3`, and the determinant-20 Gram.

Read: `Grok-rank3.md`; `RANK3-CENSUS.md`; `docs/coset-modularity.md` §§1–3; the opening of `docs/theta-transformations.md`; Kane–Kim arXiv:2211.03987v2, §§2.1–2.2 and Proposition 2.3, including the slash operator and identity (2.9); Brunault’s note, all three pages; Kumar–Purkait, arXiv:1304.6586, Lemma 3.1. The lemma states the same threshold `k μ / 24` for a form of weight `k/2`. Their writeup treats the passage to an integral-weight power as immediate and writes the leading exponent as `B+1` for a real `B`. The proof above replaces that appeal by the slash calculation and the inclusive integer bound. Sturm’s 1987 chapter was not opened.

Unrun: `code/verify_rank3_census.py` and `tests/test_rank3_census.py`. The saved census summary was not reproduced by a new classification.
