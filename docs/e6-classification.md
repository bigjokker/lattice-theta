# E6: complete verification of the parity converse

Completed 2026-09-26. Every even-p half-integral characteristic on E6 has a nonzero minimal-shell coefficient. This is a finite exact proof with complete class coverage, checked by two different enumeration algorithms. Another person or model was not needed to execute it; both verifiers are available for anyone to rerun or review.

## Result and finite coverage

There are 4^6=4096 characteristic pairs, of which 2080 have even p and 2016 have odd p. The even classes all have a complete nonzero Fourier coefficient at squared norm at most 1:

| Squared norm | Four times norm (`norm4`) | Even characteristic pairs certified |
|---|---|---|
| 0 | 0 | 64 |
| 1/2 | 2 | 1152 |
| 1 | 4 | 864 |
| Total | | 2080 |

The odd classes vanish identically by the -I parity theorem in [docs/foundations-and-exact-enumeration.md §2](foundations-and-exact-enumeration.md#2-the-stabilizer-sign-character-theorem). The independent enumeration also verifies their coefficients through squared norm 1 are zero; that bounded observation supplements the theorem and is not the reason for claiming identical vanishing. Consequently E6 has no even-p vanishing classes, and its identically vanishing characteristics are exactly the odd ones.

The artifacts are:

- [e6-even.json](../data/e6-even.json): one nonzero-shell certificate for each of the 2080 even pairs, with a shared Gram matrix.
- [e6-even-results.json](../reports/e6-even-results.json): per-class checks through the existing exact verifier, including shell size and enumeration-node counts.
- [e6-ambient-results.json](../reports/e6-ambient-results.json): the independent Euclidean enumeration report and minimal-shell classification.
- [code/verify_even_pack.py](../code/verify_even_pack.py): independently checks exact coverage of all binary representatives and passes every entry to `exact_theta.verify`.
- [code/verify_e6_ambient.py](../code/verify_e6_ambient.py): verifies the same pack with integer sum-of-squares enumeration; imports neither the first verifier nor its builder.
- [code/build_e6.py](../code/build_e6.py): constructs the certificate table from bounded exact searches.

Both saved reports identify the same certificate pack by SHA-256. Neither verifier trusts a claimed class count or treats duplicate entries as coverage. The primary verifier's successful result requires complete, nonzero shell coefficients for all classes; any interrupted enumeration leaves the overall classification unresolved.

## Lattice and characteristic conventions

Use a simple-root basis ordered as Bourbaki roots (1,3,4,5,6,2), with Gram matrix

\[
G=\begin{pmatrix}
2&-1&0&0&0&0\\
-1&2&-1&0&0&0\\
0&-1&2&-1&0&-1\\
0&0&-1&2&-1&0\\
0&0&0&-1&2&0\\
0&0&-1&0&0&2
\end{pmatrix}.
\]

Thus nodes 1–5 form a chain and node 6 joins node 3. This is the E6 Dynkin diagram with roots of squared length 2; the numbering correspondence and ambient roots agree with the [Sage reference documentation for type E](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/root_system/type_E.html#sage.combinat.root_system.type_E.CartanType.dynkin_diagram). The primary verifier checks positive definiteness and determinant 3 exactly.

Let a be the coordinates of 2xi in this basis and b the coordinates of 2delta in the dual basis. Each lies in Z^6 and is defined modulo 2Z^6. Hence a,b in {0,1}^6 represent every characteristic pair exactly once, independently of det(G). The parity is a^t b modulo 2; b is not being represented in the primal basis.

For a=0 all 64 b are even. For each of the 63 nonzero a, the map b -> a^t b is a nonzero linear functional on F2^6, so exactly 32 b are even and 32 odd. This proves the counts 64+63*32=2080 and 63*32=2016. Both programs reconstruct and compare the full expected set rather than only these counts. No quotient by an automorphism group is used, so no orbit-completeness or full-group proof is needed.

## Minimal shells give a short proof

The two complete enumerations establish the following finite shell table, for y congruent to a modulo 2:

| Classes a | Minimum y^t G y | Number of vectors in that shell | Antipodal pairs |
|---|---|---|---|
| a=0 (one class) | 0 | 1 | zero vector |
| 36 nonzero classes | 2 | 2 | 1 |
| 27 nonzero classes | 4 | 10 | 5 |

These numbers account for every one of the 64 a classes. Enumeration through y^t G y<=4 is complete, so the asserted shell minima are justified, as well as their sizes. The sphere contains 1+36*2+27*10=343 vectors, including 72 at norm 2 and 270 at norm 4 in the unshifted y coordinates.

**Lemma.** If a nonzero coset a/2+L has a minimal shell comprising an odd number of antipodal pairs, every even-p characteristic with that a has nonzero minimal coefficient.

**Proof.** Write y=2x+a. The antipodal partner -y corresponds to x'=-x-a. Its signed phase satisfies

\[
(-1)^{x'^tb}=(-1)^{x^tb+a^tb}=(-1)^{x^tb}
\]

when p=a^tb is even. Thus every pair contributes +2 or -2. An odd number of such contributions sums to 2 times an odd integer and cannot be zero. For a=0 the norm-zero vector contributes 1. The global nonzero phase i^p changes neither conclusion.

The table has one or five antipodal pairs in every nonzero a class, so the lemma proves the E6 parity converse. The general cancellation argument is symbolic; its E6 shell sizes and exhaustive coverage are certified finite computations. The 2080 explicit coefficient certificates provide an additional direct verification of the same conclusion.

## Why the second enumeration is complete

Let C be the 8 by 6 matrix whose columns are twice the simple roots in the stated order:

\[
C=\begin{pmatrix}
1&-2&0&0&0&2\\
-1&2&-2&0&0&2\\
-1&0&2&-2&0&0\\
-1&0&0&2&-2&0\\
-1&0&0&0&2&0\\
-1&0&0&0&0&0\\
-1&0&0&0&0&0\\
1&0&0&0&0&0
\end{pmatrix}.
\]

Its Gram identity is C^t C=4G; the ambient verifier computes every entry and compares it to the pack. These roots match the [documented ambient E6 simple roots](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/root_system/type_E.html#sage.combinat.root_system.type_E.AmbientSpace.simple_root) under the same reordering. The explicit matrix, rather than a network lookup, is sufficient to run the check.

Using zero-based indices, write Cy=(z0,z1,z2,z3,z4,-h,-h,h) with h=y0. Every zi is congruent to h modulo 2, and

\[
4y^tGy=z_0^2+z_1^2+z_2^2+z_3^2+z_4^2+3h^2.
\]

For y^t G y<=4, the right side is at most 16, forcing |h|<=2. The program enumerates all h in this range and recursively all five zi of the required parity under the remaining sum-of-squares budget, using integer square roots. Every possible Cy is therefore considered. This uses no LDL decomposition or rational ellipsoid recursion.

Recover the unique candidate y in the following order:

\[
\begin{aligned}
y_0&=h,& y_4&=(z_4+h)/2,\\
y_3&=(z_3+h)/2+y_4,&y_2&=(z_2+h)/2+y_3,\\
y_5&=(z_0+z_1+2y_2)/4,&y_1&=(h+2y_5-z_0)/2.
\end{aligned}
\]

Reject a tuple if any division is not integral, then directly check Cy against the original ambient tuple. Every integral y produces a tuple passing these tests, and the inverse formulas show that no two y produce the same tuple. Thus the enumeration includes all lattice vectors in the sphere exactly once. It considers 627 ambient tuples and accepts 343 lattice vectors.

For each accepted y, compute a=y modulo 2, x=(y-a)/2 and its parity vector. The program sums (-1)^(x^t b) directly for the requested coefficients, separately checking all 2080 claims and the stronger minimal-shell statement. Its only dependencies are Python standard-library modules; it shares no enumeration code with `code/exact_theta.py`.

## Reproduce the milestone

From the workspace directory:

```powershell
python code/verify_even_pack.py data/e6-even.json --output reports/e6-even-results.json
python code/verify_e6_ambient.py data/e6-even.json --output reports/e6-ambient-results.json
python scripts/run_tests.py --all
```

Both verifier commands exit 0. The test suite has 19 passing tests, covering complete E6 verification, agreement on every coset's shell counts, invalid coefficients, missing and duplicate coverage, odd or invalid representatives, wrong E6 input, interrupted enumeration and the earlier examples.

To reconstruct the table:

```powershell
python code/build_e6.py
```

To verify any one entry with the original verifier, combine it with the shared G and `"schema": "work7-theta-v1"`; [tests/test_e6.py](../tests/test_e6.py) includes that check. The pack is a compact collection of ordinary nonzero-shell certificates, not a new source of mathematical assumptions.

## What this completes

E6 moves from an earlier reported result to a complete exact computational theorem with a minimal-shell proof and a separate enumeration cross-check. This settles this lattice's parity-converse classification. It does not establish the unrestricted classification, literature priority, or a general theorem about symmetry explanations. The subsequent E7/E8 milestone is complete in [docs/e7-and-e8-classification.md](e7-and-e8-classification.md), the uniform Dn counting audit in [docs/dn-vanishing-formula.md](dn-vanishing-formula.md), and the nonsimply-laced classification in [docs/root-lattice-classification.md](root-lattice-classification.md). The completed companion manuscripts present the subsequent low-rank and rank-four results.
