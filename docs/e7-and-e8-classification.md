# E7 and E8: complete exact classifications with batched shells

Completed 2026-09-26. E7 and E8 now have complete characteristic classifications, including exact proofs for every even-parity zero class. Two different shell enumerations give identical histograms. Both lattices fail the parity converse.

## Verified results

| Lattice | All pairs | Odd pairs (zero by parity) | Even pairs | Even zero, proved by symmetry | Even nonzero, proved by a complete shell | Unresolved |
|---|---|---|---|---|---|---|
| E7 | 16384 | 8128 | 8256 | 1260 | 6996 | 0 |
| E8 | 65536 | 32640 | 32896 | 9450 | 23446 | 0 |

The nonzero certificates use these squared norms, with the global phase i^p removed from their coefficients:

| Squared norm | E7 nonzero classes | E8 nonzero classes |
|---|---|---|
| 0 | 128 | 256 |
| 1/2 | 4032 | 15360 |
| 1 | 2772 | 7830 |
| 3/2 | 64 | 0 |

Every nonzero even class is detected by its minimal shell. For these two lattices, the even classes with a cancelling minimal shell are exactly the identically zero classes, as proved by their symmetry certificates. The ternary regression example remains a counterexample to using that stopping rule on arbitrary lattices. No finite run of zero coefficients, Sturm bound or complete automorphism-group enumeration is used to prove vanishing.

## Artifacts and reproduction

- [E7 pack](../data/e7-classification.json) and [E8 pack](../data/e8-classification.json.gz): complete even-class tables with shell or symmetry evidence.
- E7 reports: [LDL](../reports/e7-ldl-results.json) and [ambient](../reports/e7-ambient-results.json).
- E8 reports: [LDL](../reports/e8-ldl-results.json) and [ambient](../reports/e8-ambient-results.json).
- [code/exceptional_batch.py](../code/exceptional_batch.py): batch coefficients, exact action tables, witness construction and coverage validation.
- [code/exceptional_ambient.py](../code/exceptional_ambient.py): separate Euclidean shell enumeration and direct sign summation, without LDL imports.
- [code/build_exceptional.py](../code/build_exceptional.py), [code/verify_exceptional.py](../code/verify_exceptional.py), and [tests/test_exceptional.py](../tests/test_exceptional.py): construction, verification and corruption checks.

```powershell
python code/verify_exceptional.py data/e7-classification.json --output reports/e7-ldl-results.json
python code/verify_exceptional.py data/e7-classification.json --backend ambient --output reports/e7-ambient-results.json
python code/verify_exceptional.py data/e8-classification.json --output reports/e8-ldl-results.json
python code/verify_exceptional.py data/e8-classification.json --backend ambient --output reports/e8-ambient-results.json
python scripts/run_tests.py --all
```

Reconstruct the packs with `python code/build_exceptional.py 7 8`. Small self-contained examples also work with the original verifier:

```powershell
python code/exact_theta.py data/e7-zero.json data/e7-nonzero.json data/e8-zero.json data/e8-nonzero.json
```

Both displayed zero examples have primal a=e1+e3 and dual b=e4, in the bases below, with p=0. Their files contain actual isometry matrices. `export_certificate` expands any transported zero entry into the same ordinary G,a,b,T format. Tests check roots, terminal nodes and the deepest nodes.

The combined suite has 30 passing tests. It includes all prior checks, every E7/E8 class in both backends, transformed/direct coefficient agreement, coverage rejection, corrupted witnesses and transport paths, bad coefficients and enumeration limits.

## Lattice conventions and coverage

Simple roots are ordered as Bourbaki (1,3,4,...,n,2), for n=7 or 8, with squared length 2. The first n-1 nodes form a chain and the last node attaches to the third. G has diagonal 2 and entry -1 on these edges. The [Sage type-E documentation](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/root_system/type_E.html#sage.combinat.root_system.type_E.CartanType.dynkin_diagram) records these diagrams and numbering. The complete matrices are in the packs and constructed explicitly by `root_gram`.

Positive definiteness and determinant (2 for E7, 1 for E8) are checked exactly. Coordinates a represent 2xi in the primal basis; coordinates b represent 2delta in the **dual** basis. Binary vectors represent every class once, independent of det(G). Mask m encodes coordinate i as bit i, with zero-based indices.

Parity is a^t b modulo two. For a=0 all 2^n choices of b are even; for each nonzero a, half are even. Thus there are 2^(n-1)(2^n+1) even pairs and 2^(n-1)(2^n-1) odd pairs. The verifier compares the entire expected even set with the supplied entries, rejecting missing or duplicate classes. Odd classes vanish by the -I theorem in PROOFS.md.

## Batching and measured scalability

For each a, enumerate every y=2x+a satisfying y^t G y<=8 once. Store counts by shell N=y^t G y and x modulo two. All dual coefficients are

\[
c_N(a,b)=\sum_{m\in\mathbb F_2^n}h_{a,N}(m)(-1)^{m\cdot b}.
\]

The primary backend computes this binary Fourier transform using integer butterflies: at each bit replace pairs (u,v) by (u+v,u-v). Induction over processed bits gives the displayed sign sum, without normalization or division. Each shell takes O(n 2^n) additions to obtain all b coefficients.

Geometric enumeration therefore runs 128 times for E7 and 256 times for E8. It is not repeated for every characteristic pair. The LDL completeness proof in docs/foundations-and-exact-enumeration.md still applies. An optional observer records the vectors for batching; private histograms are published only after enumeration completes, so interruption cannot produce a partial certificate.

The independent backend enumerates one unshifted Euclidean sphere, sorts vectors by a=y modulo 2, and sums signs directly. It uses neither LDL nor the binary Fourier transform. Both backends agree on the hash of every count indexed by (a,N,x modulo two).

| Work through squared norm 2 | E7 | E8 |
|---|---|---|
| Lattice vectors | 7113 | 26641 |
| LDL nodes across all a | 62726 | 262858 |
| Ambient tuples considered | 13045 | 47553 |

Observed verification times on this machine were approximately 0.4 seconds for E7 and 2 seconds for E8 with the primary backend, and 0.2/1.3 seconds with the ambient backend. These are measured runs, not hardware-independent guarantees; reports save the actual timings.

## Adapted Euclidean enumeration: completeness

Let C have columns twice the simple roots. The first two are (1,-1,-1,-1,-1,-1,-1,1) and (-2,2,0,0,0,0,0,0). The subsequent chain columns are 2(e_j-e_(j-1)), using ambient indices j=2,...,n-2; the last branch column is (2,2,0,0,0,0,0,0). Direct multiplication checks C^t C=4G against the pack.

Writing h=y0 gives

\[
Cy=(z_0,\ldots,z_{n-2},\underbrace{-h,\ldots,-h}_{8-n\text{ entries}},h),
\quad4y^tGy=\sum_{j=0}^{n-2}z_j^2+(9-n)h^2.
\]

Each zj is congruent to h modulo two. The bound y^t G y<=8 gives |h|<=floor(sqrt(32/(9-n))). For every such h, integer square bounds enumerate all z tuples of the required parity under the remaining budget.

Recover y0=h, then for j=n-2 down to 2 recover yj=(zj+h)/2+y_(j+1), omitting the last term at j=n-2. Finally

\[
y_{n-1}=(z_0+z_1+2y_2)/4,\qquad y_1=(h+2y_{n-1}-z_0)/2.
\]

Reject nonintegral divisions and check reconstructed Cy directly. Every integral lattice vector in the sphere produces a considered tuple and passes this inverse test; the inverse is unique, so no vector is counted twice. This adapts the E6 proof to 6/7 free z coordinates and weights 2/1 for h.

## Vanishing proofs and compact transport certificates

The simple reflection in root i is S_i=I-e_i(e_i^t G). Every supplied generator is checked to satisfy S_i^t G S_i=G and S_i^2=I. Its dual action is S_i^t. No claim about generating the full automorphism group is assumed.

For witness discovery, write a'=S_i a modulo two and b'=S_i^t b modulo two, with S_i^t b=b'+2mu. The shift laws give

\[
\Theta[a/2,b/2]=(-1)^{a'^t\mu}\Theta[a'/2,b'/2].
\]

A closed path with sign product -1 gives a stabilizing isometry with epsilon=1. Breadth-first exploration finds such a path in each component of the zero-coefficient candidates. The builder multiplies the exact word into T and checks the witness with the original verifier. Failure to find a witness would leave a component unresolved.

Each lattice's even zero set is one connected component under these generators. Its pack records one root witness, an isometry tree reaching every other zero class, and a matching node reference for every zero entry. The verifier checks the root witness with `exact_theta.verify`, each edge by exact primal/dual actions, and every class reference. Parents must precede children, preventing circular proofs.

Inductively, each tree vertex has zero theta because an isometry and representative shifts relate it to its proved-zero parent by a nonzero scalar. Equivalently, if U transports the root to a vertex, U T U^-1 is a stabilizer witness with epsilon=1 there; this matrix can be exported to the original verifier. These trees give exact proofs for all 1260/9450 zeros. Initial shell cancellations are checked only for consistency.

## E7 glue and the E8 quadratic criterion

For E7, primal coordinates of the dual vector b are G^-1 b. Every even zero has a nonintegral coordinate, so 2delta lies in E7* outside E7. Also 2G^-1 is integral; changing b by twice a dual coordinate vector therefore preserves primal-lattice membership, making this condition well defined on classes. It is necessary for these zeros, not sufficient for every characteristic on the glue to vanish.

For E8, unimodularity identifies b with primal c=G^-1 b. On V=E8/2E8 define

\[
q(v)=\|v\|^2/2\pmod2,\qquad B(v,w)=\langle v,w\rangle\pmod2.
\]

Evenness makes q well defined; norm expansion gives B as its polar form. For the characteristic B(a,c)=a^t b=p **modulo two**. Complete classification verifies

\[
\Theta_{E8}[\xi,\delta]\equiv0\iff
p\text{ odd},\quad\text{or}\quad a,c\text{ span a totally singular two-plane}.
\]

For an even pair, the second condition means a,c are nonzero and distinct with q(a)=q(c)=0; B(a,c)=0 already follows from parity. The checker compares this condition against every proved even verdict and checks q(a+c)=0 for every plane.

There are 135 nonzero singular vectors and 1575 totally singular planes. Each plane has 6 ordered independent pairs, giving 9450 zeros. The transport tree puts all ordered pairs in one orbit of the supplied isometries, hence one orbit under Aut(E8).

Nondegeneracy follows from det(G)=1. Masks 5,17,34,65 form an independent totally singular four-space: all 16 vectors in their span have q=0. Totally isotropic dimension in a nondegenerate eight-dimensional polar space is at most four; this exhibits Witt index four and certifies plus type. Both reports record the basis and plane count.

## Status and next work

Combining this with the general An proof, the uniform Dn counterexample and the E6 theorem settles the parity converse for irreducible ADE root lattices: exactly An and E6 satisfy it. The exceptional counts and E8 criterion now have complete exact certificates.

The subsequent uniform Dn counting audit is complete in [docs/dn-vanishing-formula.md](dn-vanishing-formula.md), with a general proof and two exact checks of every even class in D2-D10. The nonsimply-laced milestone is complete in [docs/root-lattice-classification.md](root-lattice-classification.md); that note also qualifies the D8 two-coset explanation for E8, which needs a separate quarter-coordinate sector. The unrestricted classification, literature priority and minimum higher-order witness claims remain separate tasks.
