# Proof audit: fan coverage and the step from Theorem S to minus type

2026-10-23. Grok. Astra's task 3 in `RANK4-SIX-CELLS-QUICK-REVIEW-ASTRA.md`,
carried out from `BRIEF_20261023_Grok_literature_and_proof_audit.md`.

Read first: `RANK4-SIX-CELLS-QUICK-REVIEW-ASTRA.md` and
`RANK4-SIX-CELLS-CLASSIFICATION.md`. This audit does not replay the 502,740
face/marking cases. That replay is Grok Bot's brief. No determinant search
and no Magma were run. Nothing under `results/` was written.
`results/rank4-secondary-fan/regenerate_v1.py` was not run: it writes
`results/rank4-secondary-fan/fan.json`.

**Verdict.** The three cones used for Theorem S are Voronoi's three primitive
secondary cones. The constancy of minimal vectors on each face, and the chain
from Theorem S to the minus-type statement, check out. Theorem S itself stays
a candidate until the face/marking census is signed off.

Checking scripts, all under this directory:

| Script | What it checks |
|---|---|
| `check_proof_audit.py` | Recomputes the three cones, the v1 walk, the principal-domain rays, the mod-4 determinant theorem, Gamma, the det-21 witness, and the 210 markings. Output: `check_proof_audit.json` (112.6 s). |
| `match_primitive_rays.py` | Maps each recomputed cone onto Δ, Δ′, Δ″. Output: `primitive_match.json`. |
| `adjacent_cones.py` | Crosses two facets of Δ and names the new ray. |

`fan.json` is not on disk. The cones below were recomputed from the three
central forms recorded in the classification note, by the current
`secondary_cone` in `rank4_secondary_fan.py`. That routine still proposes
Delaunay simplices with a floating triangulation. Each proposal is then
checked exactly (empty circumsphere, and total absolute determinant 120). A
call that returns has passed those checks.

## Sources

Quadratic-form coordinates throughout are the project's order
`(Q11, Q22, Q33, Q44, Q12, Q13, Q14, Q23, Q24, Q34)`, so that
`Q(v) = v^t Q v`.

| Source | What was read | What it states |
|---|---|---|
| Vallentin, Weißbach, Zimmermann, *The chromatic number of 4-dimensional lattices*, arXiv:2407.03513. PDF opened (17 pages). | §2 and §2.1–2.2, Tables 2.2 and 2.3 captions, bibliography | For `n ≤ 3` one GL(n,Z) orbit of full-dimensional secondary cones; for `n = 4` exactly three. In these dimensions every full-dimensional secondary cone is simplicial. Rays `R1,…,R12` and the three cones Δ, Δ′, Δ″ are printed in §2.1. The 52 subdivisions are the GL4(Z) orbits of those faces whose relative interior meets the positive definite cone, at most `3·2^10 = 3072` candidates. Section 2 also states 51 iso-edge domains and 52 Delaunay subdivisions. Tables: 17 zonotopal and 35 non-zonotopal, summing to 52. |
| F. Vallentin, *Sphere coverings, lattices, and tilings (in low dimensions)*, TU München, oral exam 26 Nov 2003. PDF from https://mediatum.ub.tum.de/doc/602017/602017.pdf (128 pages). | §4.4 (printed p. 50) and §4.4.1 (printed p. 51); §8.4.2 | Delone 1929 claimed 51 combinatorial types and missed one, found by Štogrin 1973; there are 52, of which 17 are zonotopes. The thesis quotes Delone, p. 161, attributing to Voronoi exactly three primitive Dirichlet domains in dimension 4. §4.4.1: Voronoi's algorithm produces three non-equivalent Delone triangulations `D1, D2, D3`. `D1` is the principal form of the first type and has ten bistellar neighbours, all equivalent to `D2`. The same section introduces extreme rays `R1,…,R12`. The printed matrices are set as lower-triangular displays; the text extraction scrambled them, so the matrices used below are the ones from arXiv:2407.03513, checked by the scripts. |
| Dutour Sikirić, Schürmann, Vallentin, *A generalization of Voronoi's reduction theory and its application*, arXiv:math/0601084v4. Journal: Duke Math. J. 142 (2008), 127–164, DOI 10.1215/00127094-2008-003. The journal pages were not opened. | Theorem 2.3, arXiv PDF p. 6 | Quoted under Task A.3. |
| Schürmann, Vallentin, *Computational Approaches to Lattice Packing and Covering Problems*, arXiv:math/0403272v3 (27 July 2005). | Example 5.2 and the following paragraph, arXiv PDF pp. 10–11; Proposition 5.1 | The secondary cone of the principal triangulation is `{Q : q_ij < 0 for i ≠ j, and row sums > 0}`. Its closure is Voronoi's principal domain of the first type. Voronoi 1908 §97, as cited there, gives a facet-to-facet tessellation of the positive semidefinite cone; Gruber–Ryshkov (European J. Combin. 10 (1989), 83–84) upgrades facet-to-facet to face-to-face. Proposition 5.1 cites Voronoi 1908 §77 for the regulator description of a secondary cone. Those Crelle pages were not opened. |
| Conway, Sloane, *Low-dimensional lattices VI*, preprint http://neilsloane.com/doc/fedorov.pdf. The preprint's first page says the journal version is Proc. R. Soc. Lond. A 436 (1992), 55–68. | PDF p. 18 | "Delone (1929, 1937–1938), as corrected by Stogrin (1973), enumerated the Voronoi cells of four-dimensional lattices. There are three primary ones and 52 in all." The journal page of that sentence was not pinned. Theorem 9 and Table I of the same preprint are the five 3-dimensional cells. |
| Dutour Sikirić, Kummer, *Iso Edge Domains*, arXiv:2102.11139. Journal form, from the chromatic bibliography: Expositiones Mathematicae 40 (2022), 302–314. | Lemma 5.2, arXiv PDF p. 8 | Quoted under Task B. |
| DLMF §21.5, equation 21.5.9, https://dlmf.nist.gov/21.5.E9 | The transformation law and the symplectic conditions 21.5.2–21.5.3, read in this audit | Used under Task C.1. |
| `RANK4-MINUS-TYPE.md` §2, `RANK4-PARITY-CONJECTURE.md` Theorems 1, 3, 4, `RANK4-SIX-MINIMA-COMMON-LATTICE.md` §§1–2, `RANK4-STAGE2-STRATEGY.md` §2 | The project arguments cited by the classification note | Audited under Task C. |

Unread primaries are listed in `bot_work/NEEDS_BROWSER.md`: Delone 1929, Štogrin 1973, Voronoi's Crelle memoir, Schürmann's AMS book, and Igusa, *Theta Functions* (1972), pp. 78–85.

## Task A. The fan against the literature

### A.1. Three simplicial classes

The count of three GL4(Z) classes of primitive L-types in dimension 4 is stated
in four places that were opened:

- arXiv:2407.03513 §2.1, attributing the computation to Voronoi's second memoir, and stating that every full-dimensional secondary cone in dimensions `n ≤ 4` is simplicial.
- Vallentin's thesis §4.4.1, which reruns Voronoi's algorithm and describes the adjacency: every neighbour of the principal cone is equivalent to the second cone.
- Conway–Sloane VI, preprint PDF p. 18: three primary Voronoi cells and 52 in all.
- The same thesis §4.4, quoting Delone's attribution of the three primitive domains to Voronoi, and recording Štogrin's correction of 51 types to 52.

Delone 1929 and Štogrin 1973 were not opened. Their role in the count is the one given by the thesis and by Conway–Sloane.

Simpliciality of the cones actually used was also checked directly. Each recomputed cone has 10 rays and 10 facets. The 10×10 matrix of ray coordinates has determinant `1`, `-1`, and `1` on the three cones, so the rays are a basis of the 10-dimensional space of quadratic forms. A cone spanned by a basis is simplicial: every subset of the rays is a face.

### A.2. Explicit match

The published cones, arXiv:2407.03513 §2.1, are

- Δ = cone{R1,…,R10},
- Δ′ = cone{R1,…,R4, R6,…,R11},
- Δ″ = cone{R1,…,R4, R6,…,R9, R11, R12}.

`R1,…,R10` and `R12` are primitive. `R11` as printed is twice a primitive matrix. Its primitive representative, in project coordinates, is

```
(2, 2, 2, 2, 1, -1, -1, -1, -1, 0),
```

the matrix `[[2,1,-1,-1],[1,2,-1,-1],[-1,-1,2,0],[-1,-1,0,2]]`. Replacing `R11` by this primitive vector does not change the cone. Summing the printed matrices without that reduction doubles the weight of one ray and produces a different interior form. The first matching attempt, against those weighted sums, correctly found no GL4(Z) isometry. The match below uses primitive rays.

`check_proof_audit.py` called `secondary_cone` on the three recorded centrals. Each call returned 10 rays, the recorded central, and a positive definite form. The v1 walk (the depth-first algorithm in `regenerate_v1.py`, reimplemented here and writing nothing under `results/`) found the same three orbits, with the same ray sets, in 105.2 s. `regenerate_v1.log`, which was read and not re-run, already reports the same three centrals and reports that all 5,256 saved cells have their generators on the stated face. That cell incidence was not replayed.

**Cone 0 is Δ.** The identity `U = I` carries its ten rays onto `R1,…,R10` in that order. The sum of `R1,…,R10` is the recorded principal central

```
[[4,-1,-1,-1],[-1,4,-1,-1],[-1,-1,4,-1],[-1,-1,-1,4]].
```

The closed principal domain `{q_ij ≤ 0 for i ≠ j, row sums ≥ 0}` has exactly these ten extreme rays (`inequality_cone` in the JSON). That is the closure described in arXiv:math/0403272, Example 5.2.

**Cone 1 is Δ′.** With our rays in the order stored in `check_proof_audit.json`,

```
U = [[-1, 1, 0, 0],
     [-1, 0, 0, 0],
     [-1, 0, 0, 1],
     [-1, 0, 1, 0]],    det U = -1,
```

and `U^t R U` runs through `R11, R1, R7, R6, R2, R9, R8, R4, R3, R10`, where `R11` means the primitive ray. `box_isometries` found 24 isometries of the primitive centrals; the script keeps the first one whose ray image is the whole published set. The action is `U^t (our form) U`.

**Cone 2 is Δ″.**

```
U = [[-1, 0, 1, 1],
     [-1,-1, 1, 1],
     [-1, 0, 0, 1],
     [-1, 0, 1, 0]],    det U = 1,
```

and `U^t R U` runs through `R11, R12, R7, R6, R2, R3, R4, R8, R9, R1`. There were 72 isometries of the centrals.

`adjacent_cones.py` is an independent view of the same fact. The form just outside Δ across the facet opposite `R5` has secondary cone with rays `{R1,…,R4, R6,…,R10}` together with the primitive `R11`. The form just outside across `R1` has the ray set of recorded cone 1. Both adjacent centrals have determinant 392, the determinant of the primitive central of cone 1.

`fan-exact.json` is a different set of representatives. Its three cones match the recorded cones by

```
cone 0: U = [[-1,0,0,0],[-1,0,0,1],[-1,0,1,0],[-1,1,0,0]],  det 1
cone 1: U = [[-1,0,0,1],[0,0,0,1],[0,-1,0,1],[0,0,-1,1]],   det -1
cone 2: U = [[-1,0,0,1],[0,0,0,1],[-1,0,1,1],[0,-1,0,1]],   det -1
```

again with `U^t (fan-exact form) U` equal to the recorded central, and with the rays carried onto the recorded rays. The sign pattern of the fan-exact centrals is the one in `fan-exact.log`.

### A.3. Coverage

Theorem 2.3 of arXiv:math/0601084 (Voronoi's reduction theory, as stated there) says: the secondary cone of a Delone triangulation with vertex set Z^d is a full-dimensional open polyhedral cone in the positive definite cone. The topological closures of these cones form a polyhedral subdivision of the rational closure of the positive semidefinite cone. Two closures share a facet if and only if the triangulations are bistellar neighbours. GL_d(Z) acts by `U ↦ U^t Δ U`, and there are finitely many orbits.

Combined with the match in A.2: every positive definite quaternary form lies in the closure of some GL4(Z) image of Δ, Δ′, or Δ″, and the form's L-type domain is a face of one of those images. A coarser Delaunay subdivision is the subdivision attached to a face of a triangulation cone. arXiv:2407.03513 §2 records the same relationship: every Delaunay subdivision is a coarsening of some triangulation, and the 52 types are read off the faces of these three cones.

The Crelle pages of Voronoi 1908 were not opened. The theorem number used here is Theorem 2.3 of the arXiv paper, which attributes the tiling to Voronoi.

### A.4. Faces and the 52 types

Because each cone is simplicial, every subset of its rays is a face. A face meets the positive definite cone if and only if the sum of its rays is positive definite: the relative interior is connected, the positive definite cone is open, and the face lies in the positive semidefinite cone. The script counted those sums.

| Cone | Nonempty ray subsets | Positive definite faces |
|---|---:|---:|
| Δ (recorded 0) | 1,023 | 728 |
| Δ′ (recorded 1) | 1,023 | 826 |
| Δ″ (recorded 2) | 1,023 | 840 |
| Total, faces counted once per cone | 3,069 | 2,394 |

This is the classification note's pair of counts, 3,069 and 2,394. Shared faces of the arrangement are counted once for each cone that contains them. The number 2,394 is a count with those repetitions. It is the input to the marking census. It is not the number of GL4(Z) orbits.

The orbit count 52 was not recomputed. It is the published quotient of these faces: arXiv:2407.03513 §2.2, with 17 orbits in Table 2.2 and 35 in Table 2.3, and the same 17 + 35 in Vallentin's thesis §4.4. Since the three cones have been identified with Δ, Δ′, and Δ″, that published quotient is the quotient of these faces. A fresh pairwise isometry search of the 2,394 forms was not run.

The 51 iso-edge domains are a coarser partition. One of them contains two secondary cones. Task B uses that.

## Task B. Minimal vectors are constant on a face

The classification note, section 2(b), is right, by the following checks. Let `Q` be positive definite on R^4, and let `v` be a nonzero class representative in Z^4.

**(1) Proved.** `v` is minimal in the coset `v + 2Z^4` if and only if `v/2` lies in the Voronoi cell of 0. Expanding,

```
Q(v - 2y) - Q(v) = 4 Q(y) - 4 B(v, y),
Q(v/2 - y) - Q(v/2) = Q(y) - B(v, y).
```

The two differences have the same sign. The vectors that arise are Conway and Sloane's relevant vectors: the minimal shells of the nonzero cosets of `2L`. The preprint of *Low-dimensional lattices VI* develops that language for dimension 3; the identity above is the same comparison in any dimension.

**(2) Proved.** The Voronoi cell is closed, convex, and centrally symmetric. If `v/2` lies in `Vor(0)`, then so does `-v/2`, hence `v/2` lies in `Vor(v) = v + Vor(0)`. Conversely, a point of the intersection and its reflection through `v/2` are both in the intersection, so the midpoint `v/2` is in it.

**(3) Proved.** A point of `Vor(0) ∩ Vor(v)` is the centre of an empty ball touching 0 and `v`, so 0 and `v` are vertices of a common Delaunay cell. A centre of an empty sphere through 0 and `v` lies in that intersection. The Delaunay subdivision is constant on the relative interior of a secondary cone, by definition of the cone, and constant on the relative interior of each face, where the subdivision is the coarsening cut out by that face. Therefore the set of coset-minimal vectors is constant on each such relative interior. On a proper face the subdivision coarsens, so the minimal set can grow. The classifier recomputes the set on each face. The numerical check "0 exceptions in 2,394 faces" is a second look at the same fact. It is not the proof, and it was not re-run here.

**Lemma 5.2 of arXiv:2102.11139** says something narrower, and says it correctly. The map Θ from positive definite forms to R^{2^n − 1}, recording one closest-vector value for each of the `2^n − 1` nonzero parity classes, is linear and injective on every iso-edge domain. The proof shows that if `Del(A)` is a triangulation, a chosen closest lattice point to `c` remains a closest lattice point throughout the iso-edge domain, so the value `A'[v_c − c]` is linear, and by continuity linear on the closure. Injectivity uses their Lemma 3.1.

That is a statement about the minimum values. An iso-edge domain can contain more than one secondary cone (51 domains against 52 subdivisions, arXiv:2407.03513 §2). Across that extra wall the value stays linear while the set of minimizers can grow. Multiplicities are constant on secondary-cone faces, by (1)–(3). The common-lattice note cites Lemma 5.2 for linearity of the parity-coset minima on an iso-edge domain, including the closure. That citation matches the values. The classification is right to test multiplicities on secondary cones.

## Task C. From Theorem S to the minus-type statement

### C.1. Plane invariance — proved

`RANK4-STAGE2-STRATEGY.md` §2 builds, for an even rank-4 Gram `G` of odd determinant `D`,

```
H = D G^{-1} = adj(G),    δ = (3D + 1)/2,
Gamma = [[2I, 3G], [H, δ I]].
```

**Determinant mod 4.** `check_proof_audit.py` enumerated every symmetric 4×4 matrix over Z/4Z with even diagonal: `2^4 · 4^6 = 65,536` matrices. Of these, 28,672 have odd determinant, and every one of those determinants is `1 (mod 4)`. None is `3 (mod 4)`. The determinant modulo 4 depends only on the matrix modulo 4, so every even rank-4 Gram of odd determinant satisfies `D ≡ 1 (mod 4)`. Then `δ` is an even integer.

**Symplectic check.** On the Cartan Gram of `A4` (determinant 5) and on the determinant-21 control, the script found `det Gamma = 1` and `Gamma J Gamma^t = J` for both

```
J = [[0, I], [-I, 0]]    and    J = [[0, -I], [I, 0]].
```

The second matrix is the negative of the first, and the identity is homogeneous of degree 1 in `J`, so the two checks are one fact. The block identity `A^t D − C^t B = I` holds on both examples, and holds in general because `2δ − 3D = 1` and `GH = DG = D I`. This is the condition DLMF 21.5.2–21.5.3 asks of a symplectic matrix of determinant 1. The script also checked `A^t C = 2H` and `B^t D = 3δ G`, both symmetric.

**Characteristic map.** DLMF 21.5.9, with integer characteristics `α = a/2`, `β = b/2`, gives

```
a_new = δ a − H b + diag(δ H),
b_new = −3 G a + 2 b + diag(6 G),
```

read modulo 2. Both diagonal corrections are even: `δ` and the diagonal of `G` are even, and `H` is an integer matrix. Since `D` is odd, `H b ≡ c (mod 2)` where `c = G^{-1} b`. Since `δ` is even and `−1 ≡ 1 (mod 2)`,

```
a_new ≡ c,    b_new ≡ G a    (mod 2).
```

The script checked all 256 pairs `(a, b)` in `(Z/2Z)^4 × (Z/2Z)^4` on both sample Grams. All 256 map as this rule says.

Translation `τ ↦ τ + 1` generates the rest of the ordered bases of the plane `P = span(a, c)`. It is the integer characteristic shift in the same DLMF section. Ordinary inversion `τ ↦ −1/τ` changes the lattice and is not the argument.

**The factor κ.** DLMF 21.5.9 writes the transformed series as `κ(α, β, Γ)` times a non-vanishing holomorphic factor (`det(CΩ+D)` on the Siegel upper half space, and the exponential) times the original series. κ is independent of `z` and `Ω`. The pages of Igusa that prove the formula were not opened. Vanishing transfers anyway: the identity sends the zero series to the zero series, and the inverse modular substitution is another instance of the same identity, so the six series of `P` vanish together. κ itself is nonzero because the formula at the zero characteristic, `z = 0`, `Ω = i I`, equates a classical theta value (a sum of positive terms) to `κ` times a non-zero factor.

### C.2. Six series and the marking — proved

`RANK4-SIX-MINIMA-COMMON-LATTICE.md` §§1–2 was read and the algebra agrees.

`K = {y ∈ L : y mod 2L lies in P}` has index 4 and contains `2L`. For `x, y ∈ K` the product `B(x, y)` is even because both classes lie in the isotropic plane `P`, so `2H`, the Gram of `K`, is even and `H` is integral. The index formula `det(2H) = 16 det H = 16 det G` gives `det H = det G`.

Each signed half is `θ_{v+2K} − θ_{w+2K}` at every coefficient: the sign is constant on each coset of `2K` and opposite on the two cosets, because `B(c, ·)` is non-zero on the class that separates them. Equal minimum and equal multiplicity of those two ordinary cosets is exactly a balanced half. In the coordinates of `H` this is the hypothesis of Theorem S.

The marking is `τ(v) = v + j(π(v))` on `V = K/2K`, with `W = 2L/2K` and `π : V → P`. Then `τ^2 = I`, `im(τ − I) = ker(τ − I) = W`, and `N = τ − I` has rank 2 and `N^2 = 0`.

The count of such markings is 210. The Gaussian count is `(15·14)/(3·2) = 35` planes `W`, times `|GL2(F2)| = 6`. The script enumerated all `2^{16} = 65,536` matrices over F2 and found exactly 210 with `N^2 = 0` and rank 2. For every one of those, `im N = ker N` follows from the rank.

The 19-matrix and 10,522-plane checks in that note are certificates for those inputs. They are not the general identity, which is the algebra above.

### C.3. Extension to a symmetry witness — proved

This fills the sentence in the classification note, section 1, step 4. Let `T ∈ GL(K)` satisfy `T^t H T = H` and `T ≡ τ (mod 2)` on `K/2K`.

`τ` fixes `W` pointwise and fixes every class of `P = K/2L`. For `x ∈ 2L` the class of `x` in `K/2K` lies in `W`, so `Tx = x + 2k` with `k ∈ K`. Both `x` and `2k` lie in `2L`, so `T` preserves `2L`. Define `h(y) = T(2y)/2` for `y ∈ L`. Then `h(y) ∈ L`. For `y ∈ K`, `T(2y) = 2 T(y)`, so `h` restricts to `T`. The same formula with `T^{-1}` is the inverse, and

```
B(h y1, h y2) = (1/4) B(T(2 y1), T(2 y2)) = (1/4) B(2 y1, 2 y2) = B(y1, y2),
```

because `T` preserves the Gram `2H` of `K`. Thus `h ∈ Aut(L)`. Since `τ` adds an element of `W ⊂ 2L`, `h` fixes every class of `P`, including `a`, `c`, and `a + c`.

The witness sign is the formula in `RANK4-PARITY-CONJECTURE.md`: `ε(g) = B((g a − a)/2, c)`. Take a representative `α ∈ K` of the class `a`. Then `Tα − α = 2λ + 2k'` where `2λ` represents `j(a) ∈ W` and `k' ∈ K`. The pairing `W × P → F2` is perfect, and `j(a)` is defined to be the unique nonzero element of `W` that annihilates `a`. The annihilator of the line `F2 a` is that line in `W`. It is distinct from the annihilator of `c`, because `j` is bijective on the three nonzero classes and `c` is independent of `a`. So the pairing of `j(a)` with `c` is 1. The correction `B(k', c)` is even because the class of `k'` lies in `P` and `P` is isotropic. Therefore `ε(h) ≡ 1 (mod 2)`.

This is the definition of a symmetry witness in `RANK4-MINUS-TYPE.md` §2: `h` fixes `a` and `c` modulo `2L`, and `ε(h) = 1`, on an even characteristic with `a` and `b` both nonzero.

The same `h` serves every ordered pair of distinct nonzero classes of `P`, because each such pair determines the same `τ`.

### C.4. Theorem 5 — proved by the text that was read

`RANK4-MINUS-TYPE.md` §2 was read. The argument is complete.

- Lemma A: an element of order 4 in `Sp4(F2)` has a single nonzero fixed vector. The written proof is the Jordan argument given there (nilpotence of `ḡ − I` in characteristic 2, exclusion of a `3+1` block by the alternating form). The finite enumeration of all 720 symplectic matrices is a separate check. It was not re-run. The written proof does not use it.
- Lemma B: an element of order 8 in `GL4(Z)` has minimal polynomial dividing `x^4 − 1` and divisible by the 4th cyclotomic polynomial `x^4 + 1`, hence characteristic polynomial `x^4 + 1`, hence `g^4 = −I`.
- Theorem 4 of `RANK4-PARITY-CONJECTURE.md`: `g^2 = −I` forces the fixed plane `(g+I)L/2L` to be totally singular, hence plus type (`det ≡ 1 (mod 8)`), and every ordered pair of distinct nonzero classes there is a symmetry zero with witness `g`.
- Theorem 3: an involution cannot be a witness on an even, indecomposable, odd-determinant rank-4 lattice. The corollary `ε(g^m) = m ε(g)` puts the 2-part of any witness in order 4 or 8.
- Theorem 5 then splits into image order 4 (Lemma A forces `a = c`, and Theorem 1 says that series does not vanish), order 8 (Lemma B produces a square root of `−I`, whose mod-2 image cannot be trivial), and order 4 with image order at most 2 (either `g^2 = −I`, the plus-type case, or a splitting involution, which contradicts indecomposability).

Minus type has no automorphism squaring to `−I`, by Theorem 4. So an even, orthogonally indecomposable, odd-determinant, minus-type lattice has no symmetry witness.

### C.5. Scope, and the det-21 control

Odd characteristics vanish for a separate reason (the odd nullwert; Theorem 1 of the parity note, with witness `−I`). They are not symmetry witnesses in the §2 definition, which asks for an even characteristic with both vectors nonzero. The printed Corollary 2 of the classification note already puts them outside the statement.

The determinant-21 lattice in `RANK4-STAGE2-STRATEGY.md` §3,

```
G = [[2,1,0,0],[1,2,0,0],[0,0,2,1],[0,0,1,4]],
a = (1,0,1,0),  c = (0,1,0,1),
```

has determinant 21 ≡ 5 (mod 8) and is decomposable. The script checked `T = diag(−1,−1,1,1)`: `T` preserves `G`, `T^2 = I`, `T` fixes `a` and `c` modulo 2, `B(a, c)` is even, and `ε(T) = −1`, hence `ε(T) ≡ 1 (mod 2)`. It is a symmetry witness. Theorem 5 does not apply, because the lattice is decomposable. The example shows why indecomposability stays in the final statement, and it agrees with the sentence that a vanishing nontrivial even characteristic is a symmetry zero whether or not the lattice is indecomposable.

The printed Corollary 2 is the right statement, once Theorem S is available. Written out:

> Let `L` be even, positive definite, of rank four and odd determinant, and let `(a, b)` be a nontrivial even characteristic: `a` and `b` both nonzero modulo `2L`, and `a·b` even. If the characteristic theta series vanishes identically, then the isotropic plane `P = span(a, c)`, `c = G^{-1} b (mod 2)`, admits a symmetry witness in `Aut(L)`.
>
> If in addition `L` is orthogonally indecomposable, the witness exists only in the situation of Theorem 5: plus type, with some `h` in `Aut(L)` satisfying `h^2 = −I` whose fixed plane contains `a` and `c` as distinct nonzero classes. In particular, an orthogonally indecomposable lattice with `det ≡ 5 (mod 8)` has no identically vanishing theta series of a nontrivial even characteristic.

Odd characteristics stay outside this statement. Decomposable lattices can have other witnesses, as the determinant-21 involution does. "Level 1" in the project's sense means a symmetry zero. It is a different invariant from the level of the even Gram.

## What this does to Theorem S

The hypotheses of Theorem S are GL4(Z)-equivariant. By A.2 and A.3 it is enough to test the positive definite faces of Δ, Δ′, and Δ″. By A.4 there are 2,394 such faces, counted with the repetitions of the arrangement. By Task B the minimal vectors of each parity class are constant on the relative interior of each face. The remaining step is the one this brief assigns to Grok Bot: for each of those faces and each of the 210 markings, either the six equalities of minimum and multiplicity cut out a cell that carries a matrix `T`, or the cell misses the relative interior and is handled on a smaller face. The classification note reports 5,256 cells, all certified, and Astra's quick review checked the saved witnesses algebraically. Neither that census nor Astra's recheck was repeated here. Until that census is signed off, Theorem S stays a computed candidate, and so does Corollary 2.

The old floating ellipsoid bound in an earlier `secondary_cone` is not the coverage argument used here. Coverage is the identification with Δ, Δ′, and Δ″ plus Theorem 2.3.

## Open gaps, by importance

1. **The 502,740 face/marking dispositions.** This is the remaining obstruction to Theorem S. It is Grok Bot's brief. The saved-witness recheck and the exact empty-sphere repair are Astra's, and were not repeated.

2. **Igusa, pp. 78–85, for the proof of DLMF 21.5.9.** The audit uses the formula as printed in DLMF, including the shape of `κ`. Vanishing of the six series does transfer from the identity, as written in C.1. A reader who wants the multiplier from Igusa rather than from DLMF still has those pages to open.

3. **Unread primaries for the classical count.** Delone 1929, Štogrin 1973, Voronoi's Crelle memoir, and Schürmann's AMS book. The statements used in Tasks A and B were read in the sources in the table above. Details are in `bot_work/NEEDS_BROWSER.md`.

4. **The orbit count 52 was not recomputed.** It is inherited from arXiv:2407.03513 and from Vallentin's thesis through the cone identification. A direct GL4(Z) quotient of the 2,394 faces was not run.

5. **The 5,256-cell incidence in `regenerate_v1.log` was not replayed.** The log says the generators lie on the stated faces of the recorded cones. This audit recomputed those cones and matched them to Δ, Δ′, and Δ″. It did not reopen `cells.json`.

6. **Two certificate-edge remarks from Astra's quick review, not re-tested here.** The certificate checker accepts a sign flip of the Gram, and `process_face` can emit a normal record when two samples of the minimal set disagree. Both sit in the census pipeline of gap 1.
