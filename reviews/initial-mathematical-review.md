# Source audit and the rank-two parity theorem

Prepared 2026-09-27 from the brief in `Grok.md`. This note does not change
`docs/foundations-and-exact-enumeration.md`, the certificates, or the census artifacts. Codex's norm4=16 run
and any later norm4=25 rerun are left untouched.

Three kinds of statement are separated below.

- **Verified source finding.** A page that was actually read, with the
  normalization comparison written out.
- **Proved here.** The rank-two theorem in §2. Its hypotheses are the
  project's: positive definite integral lattice, half-integral characteristics,
  identical vanishing in `tau`.
- **Unchecked lead.** A plausible reference that was not opened, or was opened
  only through a secondary account. An unchecked lead is not a novelty claim
  and not a citation to put in a manuscript.

No substantive error turned up in the root-lattice classification, the E6
antipodal sign lemma, or the rank-at-most-two determinant-at-most-24 census
conclusions. The review did not rerun the E7, E8, or Dn certificate streams.

## 1. Literature audit

### 1.1. What "vanishing" means in this project

Throughout, vanishing means that the series

```
Theta_L[xi, delta](tau) = sum_{x in L} exp(pi i tau ||x+xi||^2) exp(2 pi i <x+xi, delta>)
```

is the zero function of `tau` in the upper half-plane. One nonzero complete
shell coefficient rules that out. A long initial string of zero coefficients
does not prove it. Vanishing of a Siegel theta constant at one period matrix,
or vanishing of `theta(z, tau)` at one point `z`, is a different statement and
is marked as such whenever a source only proves that.

### 1.2. Comparison table

| Our claim | Source actually read | Where | What that source assumes and concludes | Comparison |
|---|---|---|---|---|
| `K_k(s;2k) = ((-2)^k / k!) * product_{j=1}^k (s-2j+1)`, as the coefficient `[t^k](1-t)^s (1+t)^{2k-s}` | Chihara and Stanton, *Zeros of generalized Krawtchouk polynomials*, J. Approx. Theory **60** (1990), 43–57. Read in the author's PDF, not the Elsevier typeset file | Remark immediately after Theorem 3.4, author-PDF p. 5; generating function (2.2) and explicit sum (2.3), author-PDF p. 2; Proposition 6.2, author-PDF p. 10 | `k_n(x,q,N)` has `sum_n k_n z^n = (1+(q-1)z)^{N-x}(1-z)^x` and `k_n(x,q,N) = sum_j (-1)^j (q-1)^{n-j} C(N-x, n-j) C(x, j)`. The remark states that the zeros of `k_N(x,2,2N)` are `1,3,...,2N-1`. Proposition 6.2 is the `q`-analogue of "the trivial zeros of `x` odd for `k_N(x,2,2N)`" | **Same polynomial, not the same displayed formula.** See §1.3. No lattice and no `A_n` shell appears |
| The `A_n` minimal shell equals that Krawtchouk value, and even `p` makes it nonzero | Same paper; also the public MathOverflow question 515148, read in full on 2026-09-27, including both comments by Walter Speziali | MO update of 13 Sep 2026 and comment of 14 Sep 2026 | The MO post itself derives the shell evaluation and then records the Chihara–Stanton zero statement. The 14 Sep comment says Nebe replied that *Self-Dual Codes and Invariant Theory* gives coset theta series, but not as code weight enumerators, and that the looked-for identification is not on p. 261 | **Relevant background for the polynomial. The lattice interpretation was not found in the paper.** The book named in the comment was not opened here, so that negative report stays an unchecked lead |
| Odd `p` forces `Theta` to vanish identically in `tau` | Osorio and Vázquez-Mozo, arXiv:hep-th/9511149v2, Appendix A, local file `9511149v2.pdf`, PDF pp. 33–34 | Paragraph after (A.4); also (A.5) | Running hypothesis for the transformation laws and the parity paragraph: `Lambda` even and self-dual. Characteristics satisfy `2a, 2b in Lambda`. The series in the elliptic variable `v` is odd iff `4 a·b` is odd, hence vanishes at `v=0`. Formula (A.5) is a zero at the single point `v = tau delta_1 + delta_2` | **Relevant background, narrower hypotheses.** For even unimodular lattices, `4 a·b` is the project's `p`, and vanishing at `v=0` for every `tau` is identical vanishing in `tau`. The appendix does not treat `2 delta in L*` when `L*` is larger than `L`, odd lattices, or the converse. (A.5) is not identical vanishing of a `tau`-series |
| Odd half-integral characteristics give vanishing theta-nullwerte | Secondary expositions only: Poor, notes citing Igusa and Mumford; a lecture note in Mumford's notation recording `theta_{11}(0)=0` because `theta_{11}` is odd. Igusa, *Theta Functions*, and Mumford, *Tata Lectures on Theta* I, were not opened | Poor, Lemma 1.1.6 and the following paragraph on thetanullwerte; lecture note, parity bullet for `k,l in {0,1}` | For the Riemann theta function on the Siegel space, an odd characteristic makes `theta(z)` odd in `z`, so the nullwert `theta(0, Omega)` is zero for every period matrix. Even nullwerte are not asserted to be nowhere zero; their vanishing on special `Omega` is the Schottky/hyperelliptic theory | **Same mechanism for `Z^g`, different function.** Identical vanishing of `theta(0, Omega)` on all of Siegel space is stronger than, and not the same as, identical vanishing of the project's one-variable series at `tau` times a fixed Gram matrix. A product of two odd Jacobi thetas is an even genus-2 characteristic and is identically zero on diagonal period matrices, which is the project's orthogonal-sum counterexample for `Z^2`. The classical literature states odd implies zero at `z=0`. It does not state the parity converse |
| Stabilizer sign character: an automorphism fixing both characteristic classes and carrying sign `1` forces identical vanishing | Searches for this statement on lattice theta series, plus the sources above | — | The `-I` case is the odd-character argument. The extension to a general stabilizer element was not found | **Still unverified as a prior theorem.** Absence from the sources opened here is not a novelty claim |
| Parity converse for root lattices: exactly `A_n` and `E_6`, with `C_3=A_3` and `G_2=A_2` under short-root length 2 | Conway–Sloane was not opened. Elkies's 2009 AWS notes, *Theta functions and weighted theta functions*, were read at the `D_n` product formulas (his (36)–(38)) | Elkies notes, pp. 11–12 of the draft PDF | Ordinary, unphased theta series of `D_n`, `D_n^*`, and `D_n^+` as combinations of Jacobi `theta_2, theta_3, theta_4`. No half-integral characteristic and no converse | **Relevant background for the unphased root-lattice series. Inapplicable to the classification** |
| Rank-two theorem of §2 | Binary theta series were checked through secondary and expository sources: Zagier's weight-1 theta series of a positive definite binary form, as presented in an Oregon State thesis quoting that development; Kani, *The Space of Binary Theta Series* (abstract and basis theorem); Ehlen, *Vector valued theta functions associated with binary quadratic forms* (abstract) | Those pages | These are generating functions `sum r(n) q^n` for representation numbers. The constant term is 1, so the series does not vanish. They are modular of weight 1 | **Inapplicable to the phased series.** They do not discuss which even-`p` characteristics vanish. No prior statement of the indecomposability criterion was found; that negative search is not a novelty claim |
| Index-two coset difference `Theta = i^p (theta_{xi+L_0} - theta_{xi+t+L_0})` | Kane and Kim, *Theta series of ternary quadratic lattice cosets*, arXiv:2211.03987v2, read on ar5iv (abstract, §§1–2, Theorem 1.2, Proposition 2.3). Published version: Selecta Math. New Ser. **32** (2026), article 5, DOI `10.1007/s00029-025-01110-0`. The published PDF was not read line by line; the Springer abstract and the arXiv text agree on the theorem | arXiv HTML, Theorem 1.2 and (2.5)–(2.9) | Positive definite lattices, cosets `aL+nu` of conductor `a`, scaled to be integral. Theta series `sum_{x in aL+nu} q^{Q(x)}` with **no sign character**. For rank 3, this series splits into an Eisenstein series, a unary theta series, and a cusp form orthogonal to unary thetas, equal respectively to the genus average, the spinor-genus defect, and the class defect. Weight `3/2`, level `4 N_L a^2` | **Relevant background, inapplicable to the difference.** Their cosets are ordinary norm-counting series of `aL+nu`. Our identity is an equality of two index-two cosets of a single lattice, shifted by a half-integral `xi`, with opposite signs. Agreement of representation numbers inside a genus does not produce a pair `(L,a,b)` whose difference vanishes while epsilon is trivial. The paper supplies no finite-determination theorem for our phased series |
| Equal ordinary theta series need not come from an isometry, so coset equality need not come from an automorphism of `L` | Not read in the primary papers. Standard secondary accounts (Schiemann's dimension theorem as stated in expository theses and on MathStackExchange) | — | Schiemann: positive definite lattices in dimensions 1, 2, 3 are determined by their theta series; dimension 4 has isospectral non-isometric pairs. Milnor, Kneser, Kitaoka give higher-dimensional pairs. Witt: `E_8 perp E_8` and `D_16^+` have the same unphased theta series | **Relevant background, not a counterexample.** These are equal theta series of two lattices, not of the two cosets `xi+L_0` and `xi+t+L_0` inside one lattice. Primary pages were not read, so this row stays a lead for the manuscript's discussion, not a checked citation |
| Assaf–Nebe–Rains–Sloane and the orthogonal modular-form computations suggested in `NEXT-STEPS.md` | Not read | — | — | **Unchecked lead** |

Stable links for the rows that were read:

- Chihara–Stanton author PDF: <https://www-users.cse.umn.edu/~stant001/PAPERS/krlaur.pdf>
  DOI of the journal article: <https://doi.org/10.1016/0021-9045(90)90072-X>
- MathOverflow question: <https://mathoverflow.net/questions/515148/when-is-a-lattice-theta-function-with-half-integral-characteristic-identically-z>
- Kane–Kim arXiv: <https://arxiv.org/abs/2211.03987>
  DOI: <https://doi.org/10.1007/s00029-025-01110-0>
- Osorio–Vázquez-Mozo: <https://arxiv.org/abs/hep-th/9511149> (v2, the file already in this directory)
- Elkies notes: <https://people.math.harvard.edu/~elkies/aws09.pdf>

The Elsevier PDF of Chihara–Stanton was not opened, so the remark's journal page inside 43–57 is not renumbered here. In the author PDF the remark is the sentence after the proof of Theorem 3.4 on page 5, and it is immediately followed by Theorem 3.5. That placement is what the MathOverflow update describes.

### 1.3. Chihara–Stanton against the identity in `docs/foundations-and-exact-enumeration.md` §6

Their generating function (2.2) at `q=2` and length `N=2k` is

```
sum_{n=0}^{2k} k_n(s, 2, 2k) z^n = (1+z)^{2k-s} (1-z)^s.
```

The coefficient of `z^k` is therefore `[t^k](1-t)^s (1+t)^{2k-s}`, which is the project's `K_k(s;2k)`.

Their explicit formula (2.3) has leading term in `s`, as a polynomial of degree `n`, equal to `(-1)^n q^n / n!`. At `q=2` and `n=k` this is `(-2)^k / k!`. The same leading coefficient is computed directly in `docs/foundations-and-exact-enumeration.md` §6. The two polynomials agree.

The remark after Theorem 3.4 does not display the product. It says: the hypothesis `r < N/2` cannot be dropped, because the zeros of `k_N(x, 2, 2N)` are `1, 3, ..., 2N-1`. Substituting `N = k`, those are the zeros of `k_k(s, 2, 2k)`. A polynomial of degree `k` with those zeros and leading coefficient `(-2)^k / k!` is exactly

```
K_k(s;2k) = ((-2)^k / k!) * product_{j=1}^k (s - 2j + 1).
```

So `docs/foundations-and-exact-enumeration.md` §6 rederives a factorization whose roots are stated by Chihara–Stanton and whose leading coefficient is determined by their (2.3). The paper treats the roots as the trivial family: the sentence introducing Proposition 6.2 calls Proposition 6.2 the `q`-analogue of "the trivial zeros of `x` odd for `k_N(x, 2, 2N)`". Section 4's numbered "trivial" examples are narrower (the center zero, and the degree-1 family `(1, k, 2k)`); the full odd family is the one named next to Proposition 6.2.

Nothing in the paper evaluates a lattice shell. The `A_n` sentence in the manuscript should cite Chihara–Stanton for the zeros, or for the polynomial `k_k(s,2,2k)`, and should attribute the shell interpretation to the argument in `docs/foundations-and-exact-enumeration.md` §6.

### 1.4. Recommended manuscript wording

Parity, for every positive definite integral lattice. If `p = <2 xi, 2 delta>` is odd, the series vanishes for all `tau` in the upper half-plane. The proof is the substitution `x |-> -x - 2 xi`. When the lattice is even and unimodular and both characteristics lie in `(1/2)L`, this is the vanishing at `v=0` of an odd lattice theta function in the sense of Osorio–Vázquez-Mozo, Appendix A. The argument used here does not need those restrictions.

The polynomial factorization. In the normalization

```
sum_n k_n(s,2,2k) z^n = (1+z)^{2k-s}(1-z)^s
```

of Chihara–Stanton, `k_k(s,2,2k)` has zeros `1,3,...,2k-1` and leading coefficient `(-2)^k/k!`. The resulting product is the minimal-shell coefficient of `Theta_{A_n}` computed in §6. That shell evaluation is not claimed as a property stated by Chihara–Stanton.

Root lattices. For the lattice generated by a reduced irreducible crystallographic root system, with short roots of squared length 2, the parity property holds exactly for `A_n` and `E_6`. Under that normalization the root lattices of `C_3` and `G_2` are `A_3` and `A_2`. Orthogonal sums of positive-rank lattices fail in every rank.

Rank two. A positive definite integral lattice of rank two has the parity property if and only if it is orthogonally indecomposable. Rank one is positive for every determinant. Both directions are uniform; determinant 24 is not a hypothesis.

What is left open. The symmetry converse — identical vanishing of an even-`p` series implies a nontrivial stabilizer sign on the full characteristic stabilizer — is open. Orders 4, 8, and 16 are orders of supplied witnesses. They are not certified minimum orders, and they do not prove that minimal certifying order is unbounded. The involution-glue criterion, the algebraic origin of the stored cyclotomic Gram matrices, and any claim of an exhaustive ternary or quaternary census stay out of the paper until they have proofs of their own.

### 1.5. Literature gaps that matter before submission

1. Read the Elsevier typeset pages, or a library scan, only if a journal page number for the remark after Theorem 3.4 is required. The author PDF already contains the sentence.
2. Page-check Igusa and Mumford for the odd-nullwert statement rather than citing them from the MathOverflow post's memory. The mechanism is standard; the page is not yet a checked citation.
3. The stabilizer sign character, beyond `-I`, still needs a genuine literature pass through Nebe–Rains–Sloane and the Hecke–Schoeneberg–Venkov weighted theta series. The pass done here did not find the statement and did not open those books.
4. The `A_n` shell interpretation was not found. The MathOverflow author's report of Nebe's reply is a lead, not a substitute for the book.
5. Do not cite Schiemann, Witt, or Conway–Sloane for isospectral lattices until the primary page has been read. The distinction those examples support — equal theta series of two lattices is not an isometry — is the right warning next to the coset reformulation, and it is not a counterexample to the symmetry converse.
6. Kane–Kim does not need to be forced into the argument. If it is mentioned, the sentence is that their coset theta series are unsigned representation series of ternary cosets `aL+nu`, decomposed by genus and spinor genus.

## 2. Rank two: the parity property is indecomposability

**Theorem.** Let `L` be a positive definite integral lattice of rank two. The following are equivalent.

1. Every even-`p` half-integral characteristic has `Theta_L[xi, delta]` not identically zero.
2. `L` is orthogonally indecomposable: it is not the orthogonal sum of two rank-one lattices.

Rank one is already settled and is not part of this equivalence. Every rank-one positive definite integral lattice has the parity property, by the scaling argument recorded in the census note together with the `A_1` case of `docs/foundations-and-exact-enumeration.md` §6. The negative direction of the theorem is the orthogonal-sum obstruction in `docs/foundations-and-exact-enumeration.md` §7, specialized below so that this section can be read alone. The positive direction is new.

Determinant at most 24 is evidence, previously recorded in `docs/rank-one-and-two-census.md`, and is not a step in the proof. The same proof shows that the 33 pairs left unresolved at norm4 bound 16 are not unresolved mathematically: each of them is nonzero. That does not edit the saved census. A later norm4=25 run is a certificate rerun, not a missing lemma.

### 2.1. Reduced bases and indecomposability

A `Z`-basis may be chosen so that the Gram matrix is

```
G = [[A, B], [B, C]],    1 <= A <= C,    0 <= 2B <= A,
```

with `A, B, C` integers and `Delta = AC - B^2 >= 1`. Call such a basis reduced. Existence is the usual binary reduction: subtract multiples of the first basis vector from the second until the off-diagonal entry `B` satisfies `|2B| <= A`, swap if needed so that `A <= C`, and change the sign of the second vector so that `B >= 0`. Each successful swap strictly decreases the positive integer `A`, so the process stops.

**Lemma.** In a reduced basis, `N(m,n) = A m^2 + 2 B m n + C n^2` is at least `A` for every nonzero integer vector, so the minimal norm is `A`.

**Proof.** The cross term `2 B m n` is nonnegative when `mn >= 0`. When `mn <= 0`, `2B <= A` gives `N >= A m^2 - A |m| |n| + C n^2`. If either coordinate is zero, `N` is `A m^2` or `C n^2`, both at least `A`. If both are nonzero and of the same sign, `N >= A + C > A`. If they have opposite signs, write `u = |m|` and `v = |n|`. Then `N >= A u(u-v) + C v^2`. For `u >= v` this is at least `C v^2 >= A`. For `u < v` one has `v >= 2` and `C >= A`, so `N >= A(u^2 - u v + v^2) >= 3A`.

**Lemma.** `L` is orthogonally decomposable if and only if `B = 0`.

If `B = 0`, the reduced basis is orthogonal and both diagonal entries are positive, so `L` splits. Conversely, suppose `L = Z e perp Z f` with `d = ||e||^2 <= ||f||^2 = h`. Squared norms in this basis are `d p^2 + h q^2`, so the minimal norm is `d`, achieved only by `±e` when `d < h` and only by `±e, ±f` when `d = h`. The orthogonal basis is already reduced. In any reduced basis the first vector has norm `A`, and the previous lemma says `A` is the minimal norm, so `A = d`.

If `d < h`, that first vector is `±e`. The second basis vector is then `x e ± f`, because the pair is a `Z`-basis. Its inner product with `±e` is `|x| d`. The inequality `2B <= A = d` forces `x = 0`. If `d = h`, a `Z`-basis drawn from `{±e, ±f}` uses one vector from each line, and those pairs are orthogonal.

Thus a decomposable lattice has `B = 0` in every reduced basis, and an indecomposable lattice has `B >= 1`. Being off-diagonal in a basis that has not been reduced does not prove indecomposability. The matrix `[[2, 3], [3, 6]]` is `GL(2,Z)`-equivalent to `A_2` and is indecomposable, but it is not reduced.

### 2.2. Which even pairs can vanish

Coordinates `a` of `2 xi` and `b` of `2 delta` may be taken in `{0,1}^2`. Even `p = a·b` leaves ten pairs.

- If `a = (0,0)`, then `xi` lies in `L`. The zero vector is the unique vector of norm 0 and contributes signed coefficient `1`. All four pairs with `a = 0` are nonvanishing.
- If `b = (0,0)`, then `2 delta` lies in `2 L^*`. Every sign `(-1)^{x·b}` equals `+1`. The vector `y = a` lies in the coset, so some shell is a positive integer. The three pairs with `a != 0` and `b = 0` are nonvanishing.

The remaining even pairs are

```
(a, b) = ((1,0), (0,1)),   ((0,1), (1,0)),   ((1,1), (1,1)).
```

Here `y = 2x + a`, the shell label is `N(y) = y^t G y = 4 ||x+xi||^2`, and the signed coefficient of a shell is `sum (-1)^{x·b}` over that shell, as in `docs/foundations-and-exact-enumeration.md` §4. One nonzero value of this coefficient implies that `Theta` is not identically zero: if `N` is the least label with a nonzero coefficient `c`, the function `exp(pi t N/4) Theta(it)` tends to `i^p c` as `t -> +infinity`.

### 2.3. Two shells that never cancel

**Lemma.** For every reduced positive definite binary form, the coset `y ≡ (1,0) (mod 2)` has minimal norm `A`, achieved exactly by `y = ±(1,0)`. Both vectors have sign `+1` against `b = (0,1)`. The signed coefficient is `+2`.

The same statement with the axes exchanged: the coset `y ≡ (0,1) (mod 2)` has minimal norm `C`, achieved exactly by `y = ±(0,1)`, both of sign `+1` against `b = (1,0)`, with coefficient `+2`.

In particular both series are nonvanishing whether or not `B = 0`. A diagonal lattice has exactly one even-`p` vanishing characteristic, not more.

**Proof.** Write `N(m,n) = A m^2 + 2 B m n + C n^2`. Since `0 <= 2B <= A <= C` and `B >= 0`,

```
N(m,n) >= A m^2 - A |m| |n| + C n^2
```

whenever `mn <= 0`, and the cross term is nonnegative when `mn >= 0`.

Consider `y_1` odd and `y_2` even. The vectors `±(1,0)` have norm `A`. If `|y_2| >= 2` and `y_1`, `y_2` have the same sign, or `y_1 = 0`, then `N >= 4C >= 4A`. If they have opposite signs, set `m = |y_1|` and `n = |y_2| >= 2`. Then

```
N >= A m (m - n) + C n^2.
```

When `m >= n` this is at least `C n^2 >= 4A`. When `m < n`,

```
N >= A (m^2 - m n + n^2) = A ((m - n/2)^2 + (3/4) n^2) >= 3A,
```

because `n >= 2`. So no other vector of the coset has norm `A`, and nothing in the coset is smaller.

The signs for `b = (0,1)` depend only on `y_2/2`. Both `±(1,0)` have `y_2 = 0`, so both signs are `+1`.

For `y_1` even and `y_2` odd the same estimate, with `C >= A`, shows that any vector with `|y_1| >= 2` has norm strictly greater than `C`. The only minimal vectors are `±(0,1)`. Against `b = (1,0)` the sign depends only on `y_1/2`, which is `0` for both, so the coefficient is `+2`.

### 2.4. The remaining pair

**Lemma.** On the coset `y ≡ (1,1) (mod 2)`, the value `M = A - 2B + C` is achieved by `y = (1,-1)` and `y = (-1,1)`.

- If `B >= 1`, these are the only two vectors of norm `M`, every other coset vector has strictly larger norm, and both signs against `b = (1,1)` equal `-1`. The signed coefficient is `-2`.
- If `B = 0`, the four vectors `±(1,±1)` all have norm `A+C`, their signs are `+1, +1, -1, -1`, and the minimal coefficient is `0`. This cancellation does not, by itself, prove identical vanishing.

**Proof.** For odd integers `m, n`,

```
D(m,n) = N(m,n) - M = A(m^2 - 1) + 2B(mn + 1) + C(n^2 - 1).
```

Direct evaluation gives `D(1,-1) = D(-1,1) = 0` and `D(1,1) = D(-1,-1) = 4B`.

Both coordinates are odd, so `mn != 0`. Global sign does not change `D`.

If `mn >= 1` and `B >= 1`, the cross term is at least `4B >= 4` whenever `|m| = |n| = 1`, and if `|m|` or `|n|` is at least 3 then `A(m^2-1)` or `C(n^2-1)` is at least `8A` or `8C`. Thus `D > 0` on the same-sign pairs.

If `mn <= -1`, the inequality `2B <= A` reverses on multiplication by `mn+1 <= 0`, so

```
D(m,n) >= A m (m + n) + C (n^2 - 1).
```

Set `u = |m|` and `v = |n|`. In either opposite-sign orientation the bound becomes

```
D >= A u (u - v) + C (v^2 - 1).
```

If `u >= v`, this is at least `C(v^2 - 1)`, hence at least `0`, and it equals `0` only for `u = v = 1`. Those two vectors are `(1,-1)` and `(-1,1)`. If `v >= 3`, the same lower bound is at least `8C`. If `u > v = 1`, it is at least `6A`.

If `u < v`, then `v >= u + 2` because both are odd, and `C >= A` gives

```
D >= A (u^2 - u v + v^2 - 1).
```

The quadratic `u^2 - uv + v^2 - 1` is at least `6` for these odd pairs (the minimum occurs at `u = 1`, `v = 3`). Thus `D >= 6A > 0`.

The signs against `b = (1,1)` are independent of `G`. For `y = (1,-1)`, `x = (0,-1)` and `x·b = -1`. For `y = (-1,1)`, `x = (-1,0)` and `x·b = -1`. For `y = (1,1)`, `x = 0`. For `y = (-1,-1)`, `x = (-1,-1)` and `x·b = -2`. When `B >= 1` only the first two vectors lie on the minimal shell, and the coefficient is `-2`.

### 2.5. End of the proof

If `L` is indecomposable, then `B >= 1` in a reduced basis. The three pairs of §2.2 have signed coefficients `+2`, `+2`, and `-2` on their minimal shells. Together with the seven pairs settled by `a = 0` or `b = 0`, every even class is nonvanishing.

If `L` is decomposable, then `B = 0` and `L = Z e_1  perp  Z e_2` in the reduced basis, with norms `A` and `C`. Take `2 xi = e_1 + e_2` and `2 delta = e_1^* + e_2^*`. Each rank-one factor has odd characteristic pairing `1`, so each factor series vanishes identically by the substitution `x |-> -x - 2 xi`. Absolute convergence gives

```
Theta_L = Theta_{Z e_1} Theta_{Z e_2} = 0,
```

while the total pairing is `p = 2`, which is even. This is one even vanishing characteristic, so the parity property fails.

The boundaries asked for in the brief are included. `A = C` and `2B = A` is allowed: `A_2` has `A = C = 2`, `B = 1`, minimal coefficient `-2` on the third pair, and the property holds. Imprimitive forms such as `[[4, 2], [2, 4]]` have `B > 0` and the same coefficient `-2`. Odd diagonal entries are allowed; integrality of `G` was the only arithmetic hypothesis. A first-shell cancellation occurs for `B = 0` on the third pair and is not used as a proof of vanishing.

### 2.6. Finite check, labeled as evidence

The shell predictions of §2.3 and §2.4 were compared with `exact_theta.complete_shells` on every reduced form with `1 <= A <= 20` and `A <= C <= A+25`, for all three pairs, through the predicted minimal norm. That is 9360 shells. Every shell matched the predicted norm, vector count, and signed coefficient. Separate spots outside that box, including `[[100, 50], [50, 100]]`, `[[100, 1], [1, 10000]]`, and the imprimitive `[[4, 2], [2, 4]]`, matched as well. The `GL(2,Z)` image `[[2, 3], [3, 6]]` of `A_2` is not reduced; its transported critical class still has minimal coefficient `-2`.

This check is evidence that the inequality did not drop a vector inside that range. The proof is the inequality, and it does not stop at `A = 20`.

Every reduced rank-two form of determinant at most 24 has `A <= 5`, by the estimate `Delta >= 3 A^2 / 4` in `docs/rank-one-and-two-census.md`. Those forms sit inside the checked range. The census count — 26 positive non-diagonal classes and 44 negative diagonal classes, each negative class with one even zero — agrees with the theorem. The eight unresolved pairs `diag(1,c)`, `a = (0,1)`, `b = (1,0)`, `c = 17,...,24`, have minimal norm `c` and signed coefficient `2`, which is the second lemma of §2.3. The other 25 unresolved pairs have `b = 0` and are the positive-series case of §2.2.

## 3. Review of the established proofs

Checked, and consistent with the sources and with the theorem above.

- `docs/foundations-and-exact-enumeration.md` §6 and Chihara–Stanton (2.2)–(2.3) use the same polynomial, including the leading coefficient `(-2)^k / k!`. The project's decision not to claim the factorization as new is the right one. The shell interpretation remains the lattice-theoretic step.
- The E6 count `64 + 63*32 = 2080` even pairs is the correct count of even vectors in `F_2^6`, and it does not depend on `det G = 3`. The antipodal sign computation in `docs/e6-classification.md` matches `docs/foundations-and-exact-enumeration.md` §4: if `p` is even, `y` and `-y` carry the same sign, so an odd number of antipodal pairs cannot sum to zero. The shell sizes 1, 2, and 10 are certified computations; they were not rerun.
- The uniform `D_n` formula at `n = 4` equals 9, in agreement with the sequence quoted from the older draft. The algebraic identity `2^n + 2^{n-1} = 3 * 2^{n-1}` reconciles the two ways of writing the formula. This is an arithmetic check, not a re-proof of `docs/dn-vanishing-formula.md`.
- `docs/root-lattice-classification.md` distinguishes the root lattice from the Cartan matrix and records the scaling bijection `(xi, delta) -> (xi, delta/c)`. That bijection preserves `p` and identical vanishing. The positive list `A_n`, `E_6`, `C_3`, `G_2`, with `C_3` and `G_2` contributing `A_3` and `A_2`, is stated inside that normalization. Nothing in the sources read here contradicts it or anticipates it.
- `docs/rank-one-and-two-census.md` is right to keep the 33 pairs unresolved inside the norm4=16 artifact, and right that they do not affect the lattice-level verdicts. Section 2 strengthens the lattice-level statement from "determinant at most 24" to every determinant, and it also resolves the individual series. The saved JSON should continue to say unresolved at bound 16.

No correction to a proof file is proposed. The following are not promoted by this review: full automorphism groups of the higher-rank examples, minimality of the witness orders 8 and 16, unboundedness of minimal certifying order, the glue criterion, algebraic provenance of the stored Gram matrices, and any exhaustive census in rank greater than 2.

## 4. What remains open

The symmetry converse is untouched by §2. In rank two every even vanishing characteristic that occurs is the orthogonal-sum characteristic, and `-I` on each factor, or the explicit swap of the two halves when the summands have equal rank, supplies a symmetry witness. A counterexample needs an identically zero even-`p` series whose full characteristic stabilizer has trivial sign. Rank two does not contain one.

The useful caution from the isospectral-lattice literature is the one already in `docs/foundations-and-exact-enumeration.md` §8: equality of the two coset theta series is not, by itself, an automorphism of `L`. Kane–Kim does not produce such a pair. A search for one remains a search in rank at least 3, with the full stabilizer proved and with a proof of identical vanishing rather than a shell bound.
