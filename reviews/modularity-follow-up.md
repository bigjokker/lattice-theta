# Follow-up review

Prepared 2026-09-27 from the updated brief in `Grok.md`. `reviews/initial-mathematical-review.md` is
left as the first report. No proof, certificate, census input, or census
report was edited.

Three labels are used. **Checked** means the argument or page was read and
the stated comparison was made. **Inference** means a calculation from those
pages, not a quotation. **Unresolved** means a book or page that was not
opened. An unresolved item is not a novelty claim.

## 1. Adversarial review

No incorrect statement turned up in the rank-two theorem, in the diagonal
witness, or in the manuscript's account of what the closed census claims.
The issues below are missing lines, ambiguous wording, or sources still
unread. The exceptional packs and the norm4=25 shell enumeration were not
rerun.

### 1.1. What was checked directly

The three critical shells in `docs/rank-two-classification.md` were rederived, including the
step that replaces the cross term `2B` by `A`. The bounds are strict on
every vector outside the listed minimal shell:

- odd/even coset: any other vector has norm at least `min(3A, 4C) > A`;
- even/odd coset: the difference `N - C` is at least `A u(u - v) + C(v^2 - 1)`,
  and unequal parities make this positive both for `u > v` and for `u < v`
  (where then `v >= 3`);
- odd/odd coset: `N - M = A(u^2 - 1) - 2B(uv - 1) + C(v^2 - 1)`, and
  `2B <= A` produces a lower bound that is zero only at `u = v = 1`.

These cover `C` much larger than `A`, `A = C`, `2B = A`, `B = 0`, odd
diagonal entries, and nonprimitive forms. The same shells were recomputed
with `exact_theta.complete_shells` for `(A,B,C)` equal to `(2,1,10000)`,
`(5,2,5)`, `(4,2,4)`, `(9,4,100)`, and `(1,0,17)`. Counts and signs matched
the table, including coefficient `0` from four vectors when `B = 0`, and
coefficient `-2` from two vectors when `B > 0`.

The decomposition step also checks. If `L = Ze` perpendicular `Zf` with
norms `d <= h`, every reduced first vector is a shortest axis vector. After
signs and axis labels it is `e`, and the second basis vector is `xe + sf`
with `s = ±1` because the change-of-basis determinant is `±1`. Then
`2|x|d <= d`, so `x = 0`. This uses neither `d < h` nor shortness of the
second vector. So `B = 0` in every reduced basis if and only if `L` splits,
for an arbitrary unimodular basis of the splitting.

The even-zero witness is `T = diag(-1, 1)`, not a swap of the axes.
On a diagonal Gram matrix, for `a = b = (1,1)`:

- `T` preserves `G`;
- `(T - I)a = (-2, 0)` and `(T^{-t} - I)b = (-2, 0)`;
- the dual shift is `(-1, 0)` and `epsilon = -1`, which is `1` in `F_2`.

The same matrix works for unequal axes. The coordinate swap on `diag(7,7)`
fixes both characteristics with shifts `(0,0)` and has `epsilon = 0`. Both
values were returned by `exact_theta.symmetry_data`. Vanishing of that one
class is the product of the two odd rank-one series, equivalently this
witness. The other nine even classes have a nonzero complete coefficient, so
a decomposable rank-two lattice has exactly one even zero. Rank one has none:
the even pairs are `(0,0)`, `(0,1)`, and `(1,0)`, with coefficients `1`, `1`,
and `2`. Odd classes in both ranks have witness `-I`. Existence of one
epsilon-one element is the symmetry converse; the full stabilizer is not
required. That is the argument in `docs/rank-two-classification.md` and `docs/foundations-and-exact-enumeration.md` §15, and it
has no determinant restriction.

The closed pack `data/census/rank12-det24-closed.json` was read, not
re-enumerated. It contains 94 lattices and 772 even pairs: 728 recorded
nonzero, 44 recorded zero, no other verdict. Those zeros are exactly one per
diagonal rank-two Gram matrix (44 diagonal, 26 non-diagonal). Lattice-level
counts are 50 positive and 44 negative. This matches Theorem 3 and the
numbers in `MANUSCRIPT.md` §5. It is a check of the saved labels, not a new
shell computation.

`docs/foundations-and-exact-enumeration.md` §2 and §6 were reread against the manuscript. The sign-character
proof, the inverse-transpose formula, the `A_n` leading coefficient
`(-2)^k/k!`, and the anti-palindromic vanishing for odd `s` agree. Section 6
correctly separates that polynomial from the shell evaluation. The E6, E7,
E8, Dn, and root-lattice certificate runs were not repeated. The E6 antipodal
sign identity in `docs/e6-classification.md` is the same phase computation as
`docs/foundations-and-exact-enumeration.md` §4; the shell sizes 1, 2, and 10 remain certified data, not
something rechecked here. E7 and E8 are reproducible from the packs and the
milestone notes, not from a closed formula written out in the manuscript.
That is an appropriate division for this draft.

### 1.2. Issue table

| Location | Statement | Severity | Reasoning | Correction |
|---|---|---|---|---|
| `MANUSCRIPT.md` §3, proof of Theorem 3, the sentence beginning "Proposition 1 gives" | "Proposition 1 gives the converse symmetry implication." | wording | Proposition 1 is epsilon nontrivial implies vanishing. The symmetry converse is the opposite direction, and in ranks one and two it is the list of witnesses already given. The sentence can be read as attributing the converse to Proposition 1. | Replace it by two sentences: Proposition 1 is the sufficiency, and the witnesses `-I` and `diag(-1,1)` are the symmetry converse in these ranks. |
| `MANUSCRIPT.md` §3, the sentence on `T = diag(-1,1)` | "T=diag(-1,1) has epsilon one for that pair." | missing justification | The matrix is the correct witness, and `docs/rank-two-classification.md` computes the shifts. The manuscript itself does not. The axis swap is the nearby false witness: on equal axes it lies in the stabilizer and has epsilon zero. | Copy the three congruences into the manuscript: `T` preserves `G`, both stabilizer conditions hold, and the dual shift is `(-1,0)`. |
| `MANUSCRIPT.md` §3, the two displayed lower bounds | "opposite signs give `N >= Au(u-v)+Cv^2`" and "`N-M >= Au(u-v)+C(v^2-1)`" | missing justification | Both are true. Each uses `2B <= A` on a negative cross term, and the manuscript does not write that substitution. `docs/rank-two-classification.md` writes it for `N - M` and uses the same estimate for `N`. A reader checking only the manuscript has to reconstruct the step that fails if `2B > A`. | Insert `N - M = A(u^2-1) - 2B(uv-1) + C(v^2-1)` and the comparison `2B <= A` before the lower bound. Do the same one-line comparison for the shortest-vector estimate. |
| `docs/foundations-and-exact-enumeration.md` §2, "Separate open question" | The symmetry converse is stated as open, with no rank restriction and no pointer to §15. | wording | The general question remains open in rank at least 3. As written, the paragraph also reads as if ranks one and two were still open. | Add one sentence: §15 proves the converse for every characteristic in ranks one and two. |
| `docs/foundations-and-exact-enumeration.md` §14, last paragraph | "It supplies no all-determinant binary classification or unrestricted symmetry converse." | wording | True of the finite certification in §14. False if read as the status of the project after §15. | End the sentence with "Section 15 removes the determinant restriction in ranks one and two." |
| `MANUSCRIPT.md` abstract | "Exact finite certificates resolve all exceptional root-lattice classes" | wording | In the body, E6, E7, and E8 are finite certificates, while F4 and G2 are identifications with D4 and A2 plus their own packs. "Exceptional" can be read as only the E series. | Name the five lattices, or say "the exceptional root lattices E6, E7, E8, F4, and G2". |

Nothing in the manuscript silently treats a witness order as minimal, asserts
unbounded minimal order, claims the glue criterion, treats the norm4=16 file
as complete, or states the symmetry converse in rank at least 3. The census
paragraph says the uniform theorem does not need the census, and it keeps the
bound-16 file and the bound-25 file separate. The root-lattice statement
names the short-root normalization and the list `An, E6, C3, G2`, with `C3`
and `G2` contributing `A3` and `A2`. That matches `docs/root-lattice-classification.md`. The
unimodular basis matrices in that note were not rechecked line by line.

The three corrections in the brief are right, and the current draft already
uses them except for the two manuscript gaps in the table. The swap is not a
witness. Odd-character vanishing along `Omega = tau G` is a specialization of
the classical nullwert vanishing, so the narrower hypotheses in
Osorio–Vázquez-Mozo do not make the general odd case new; the direct
substitution is still a correct self-contained proof. Unsigned coset theorems
are not wholly inapplicable; §3 records the extra hypotheses the difference
needs.

## 2. Primary-source gaps

### 2.1. Half-characteristic parity, checked as an identification

**Checked algebraic comparison, not a book citation.** With `a` and `b` the
integer coordinates of `2 xi` and `2 delta`, and `Omega = tau G`,

```
sum_x exp( pi i (x + a/2)^t Omega (x + a/2) + 2 pi i (x + a/2)^t (b/2) )
```

is exactly `Theta_L[xi, delta](tau)`. Positive definite `G` puts `Omega` in
the Siegel upper half-space. The characteristic is the half-integral pair
`(a/2, b/2)`. Odd `p = a·b` is the classical oddness condition
`4 (a/2)·(b/2)` odd.

Substituting `x |-> -x - a` multiplies the series by `(-1)^p` and preserves
norms. For odd `p` the nullwert is zero at every such `Omega`, hence along
the ray `tau G`. This is identical vanishing in `tau`. It is not the
vanishing of an even nullwert at one period matrix.

**Unresolved as a printed citation.** Igusa, *Theta Functions*, Springer,
1972, DOI `10.1007/978-3-642-65315-5`, was not opened. The Springer table of
contents shows the geometric chapter ending on page 85 and the next chapter
beginning on page 86; a secondary citation of "p. 85" for the transformation
law is therefore plausible and unchecked. Mumford, *Tata Lectures on Theta* I,
was not opened. Archive.org lists it as access-restricted. No clean page image
of the general-genus parity law was read in this pass. Genus-one even/odd
statements are abundant and are not a substitute for that page.

### 2.2. The `GL(n,Z)` block, and what epsilon is

**Inference from the standard block, with the multiplier left unresolved.**
For `U` in `GL(n,Z)`, the period matrix `Omega = tau G` transforms by

```
[ U^t , 0 ; 0 , U^{-1} ]
```

inside `Sp(2n,Z)`, because

```
(U^t Omega) (U^{-1})^{-1} = U^t Omega U = tau (U^t G U).
```

The usual diagonal correction in the characteristic transformation is built
from the products of the off-diagonal blocks. Both of those blocks are zero
here, so the correction vanishes. The characteristics move by

```
a/2 |-> U^{-1}(a/2),    b/2 |-> U^t (b/2),
```

which is the project's coordinate rule `a' = U^{-1} a`, `b' = U^t b`. That
part of the transformation matches.

If `T` also preserves `G`, this symplectic matrix fixes the ray `Omega = tau G`.
When `T` stabilizes both characteristics, the new characteristic differs from
the old by an integral vector. Reducing that vector multiplies the nullwert by
the sign whose exponent is `a^t ((T^{-t} - I)b)/2`. That sign is epsilon.
The remaining automorphy factor, a square root of `det(C Omega + D)` times a
root of unity depending on the symplectic matrix, was not identified with `1`
from a printed formula. Epsilon is the characteristic-reduction sign in this
special case. It should not be described as the whole Igusa multiplier until
that page is read. The series proof in `docs/foundations-and-exact-enumeration.md` §2 does not need the
multiplier; it is the right proof to print.

### 2.3. Nebe–Rains–Sloane and the three missing prior statements

**Read.** Nebe, Rains, and Sloane, *Codes and Invariant Theory*,
arXiv:math/0311046, introduction and §§2–3. The main theorem is about
complete weight enumerators of self-dual isotropic codes over finite rings.
Section 2 says the explicit Clifford–Weil construction is for finite alphabets
and that self-dual lattices are left to the forthcoming book. The paper has
no `A_n` shell, no Krawtchouk evaluation, no stabilizer sign for a lattice
theta series, and no binary indecomposability criterion.

**Unresolved.** Nebe, Rains, and Sloane, *Self-Dual Codes and Invariant
Theory*, Springer, 2006. Google Books will display the existence of page 261
and withholds the text. Chapter 9 is the lattice chapter, pages 249–284 in
the publisher's table of contents. Those pages were not read. The
MathOverflow comment that Nebe did not find the Krawtchouk weld in the coset
and Jacobi sections remains a second-hand lead.

**Negative search, not a novelty claim.** In the sources opened for this
follow-up and the first report, there is still no printed statement of:

- the `A_n` minimal shell as the Krawtchouk value `k_k(s, 2, 2k)`;
- a stabilizer sign character for a general lattice automorphism, beyond the
  central inversion that gives odd parity;
- the theorem that a rank-two lattice has the parity property exactly when it
  is orthogonally indecomposable;
- the root-lattice list under short-root length 2.

Binary theta series in the unphased sense, `sum r(n) q^n`, are classical and
do not vanish. They do not contain this criterion.

## 3. Unsigned coset series and the signed difference

**Checked translation.** For even `p` and `b = 2 delta` not in `2 L^*`, let
`L_0 = {x in L : <x, b> even}` and choose `t` in `L` with `<t, b>` odd. Then
`a = 2 xi` lies in `L_0`, and `2t` lies in `L_0` but not in `2 L_0`, because
`<t, b>` is odd. The cosets `a + 2 L_0` and `a + 2t + 2 L_0` of the lattice
`2 L_0` are therefore distinct. With

```
Phi_C(z) = sum_{y in C} exp(2 pi i z ||y||^2),
```

the identity `||x + xi||^2 = ||2x + a||^2 / 4` gives

```
Theta_L[xi, delta] = i^p ( Phi_{a+2L_0}(tau/8) - Phi_{a+2t+2L_0}(tau/8) ).
```

This matches `docs/foundations-and-exact-enumeration.md` §8. Every `y` in either coset lies in `L`, so
`||y||^2` is an integer when `L` is integral.

**What Kane–Kim apply to, from arXiv:2211.03987v2, Proposition 2.3 and
Theorem 1.2, as read for the first report.** Their theta series of a coset
`c M + nu`, with `nu` in an integral lattice `M` and conductor `c`, is a
modular form of weight `rank(M)/2` on

```
Gamma_0(4 N_M c^2) intersect Gamma_1(c),
```

with character `chi_{4 det M}` or `chi_{(-1)^{k/2} 4 det M}` according to the
parity of the rank. The group action can send the coset `c M + nu` to
`c M + p nu`. It fixes the series when `p ≡ 1 mod c`, which is why `Gamma_1(c)`
appears.

Each doubled coset can be written in that shape with `M = L_0`. The conductor
is the least positive integer `c` such that `c` times the shift lies in
`2 L_0`. It need not be the same for the two shifts: if `a` lies in `2 L_0`
and `a + 2t` does not, the conductors are 1 and 2. The level also depends on
`N_{L_0}`, the level of `L_0`, which is not the level of `L`.

**What is still required before a finite-determination bound.** Both series
must be placed in one common space of modular forms. That is the forms for
the intersection of the two groups, using the least common multiple of the
conductors and the level of `L_0`. The substitution `z = tau/8` then has to
be pushed through the slash operator; a modular form in Kane–Kim's variable
`q = exp(2 pi i z)` does not automatically give a form in `tau` of the same
level. Theorem 1.2, the genus and spinor-genus splitting, is stated for rank
3 only and does not apply to the difference in another rank. No Sturm bound
is available until the weight, the level after the substitution `tau/8`, the
character, the Fourier normalization, and the cusps are computed for this
difference. Nothing in that paper proves the difference vanishes.

This is background for a later rank-three search. It is not a reason to
enlarge the census.
