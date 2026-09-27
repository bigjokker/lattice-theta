# Adversarial review of rank-four feasibility

Prepared 2026-09-27 from `Grok-rank4-feasibility.md`. This note reviews the generation argument and the weight-two identity cutoff in `docs/rank-four-census-domain-and-cutoff.md`, together with the feasibility code, tests, and saved controls named there. The rank-three theorem, census, and cutoff were not reopened. No rank-four census was started, and no rank-five question was taken up.

Manuscripts, earlier reviews, inventories, certificate packs, and ZIPs were left unchanged. The feasibility packs were read and replayed; they were not rewritten.

**Checked** means a source statement read in this pass, or an exact integer, shell, group, or domain computation run in this pass. **Deduction** means a step proved from those statements and from `docs/coset-modularity.md` §§1–3, which this pass re-read. This review does not address novelty. The unrestricted rank-four symmetry converse remains open. Astra was not invoked.

## Verdicts

**Generation: accept.** Every positive definite integral rank-four lattice of determinant at most \(H\) has at least one admitted basis. The inequalities are closed at \(|\mu|=1/2\) and at \(d_{i+1}=(3/4)d_i\). Extra bases are allowed. For \(H=24\) the executable domain contains 1,510 such bases.

**Cutoff: accept.** Let \(c_N\) be the complete signed norm4 coefficient of \(F\), the coefficient of \(q^N\) with \(q=\exp(2\pi i z)\). If \(c_N=0\) for every integer \(N\) with \(0\le N\le B\), where \(B=\mu(16N_0)/6\) and \(\mu(M)=[\mathrm{SL}_2(\mathbb Z):\Gamma_0(M)]\), then \(F\) is identically zero. The constant term and the endpoint \(B\) belong to the test. The weight, the square, and the denominator 6 are the rank-four values.

**Evidence and implementation: accept.** The saved counts, full stabilizers, and cutoff verdicts describe what was computed. Candidate bases are ordered Gram matrices. The five named lattices are controls. A zero coefficient below the cutoff stays unresolved. A proper subgroup cannot stand in for the full stabilizer. These checks leave the symmetry converse open.

## 1. Generation

### 1.1. Saturated bases and successive projections

Checked construction, expanded here. Let \(L\) be a positive definite integral lattice of rank four. A shortest nonzero vector \(v\) is primitive: if \(v=mw\) with an integer \(|m|\ge 2\), then \(w\in L\) is shorter. A primitive vector is the first column of some matrix in \(\mathrm{GL}(4,\mathbb Z)\).

Let \(\pi\) be orthogonal projection along \(\mathbb R v\). Then \(\ker\pi\cap L=\mathbb Z v\). The images of the other three basis vectors are linearly independent over \(\mathbb R\) and generate \(\pi(L)\) as a group, so \(\pi(L)\) is a free \(\mathbb Z\)-module of rank three inside \(v^\perp\). It is discrete. Its shortest nonzero vector is primitive in \(\pi(L)\). Lifting a basis of \(\pi(L)\) and repeating produces a \(\mathbb Z\)-basis \(v_1,\ldots,v_4\) of \(L\) whose successive projections are shortest in the corresponding quotients.

Size reduction subtracts integer multiples of earlier basis vectors, in reverse index order. Each step is unimodular, so the basis stays saturated. Subtracting a vector from the span of \(v_1,\ldots,v_j\) does not change the projection orthogonal to that span, and the reverse order leaves the already reduced later coefficients fixed. After reduction every Gram–Schmidt coefficient satisfies \(|\mu_{ij}|\le 1/2\). When a coefficient is exactly half an odd integer, either adjacent integer leaves absolute value \(1/2\). Both choices satisfy the same closed bound.

### 1.2. The pivot inequality

Deduction. At stage \(i\), the projection of the size-reduced \(v_{i+1}\) onto the orthogonal complement of \(\mathrm{span}(v_1,\ldots,v_{i-1})\) is \(v_{i+1}^*+\mu_{i+1,i}v_i^*\). Its squared length is

\[
d_{i+1}+\mu_{i+1,i}^2 d_i,
\]

where \(d_j=\lVert v_j^*\rVert^2\). That projected vector lies in the lattice whose shortest squared length is \(d_i\). Hence \(d_{i+1}+\mu^2 d_i\ge d_i\). The bound \(|\mu|\le 1/2\) gives \(\mu^2\le 1/4\), so

\[
d_{i+1}\ge\frac34 d_i.
\]

Equality is attainable: \(|\mu|=1/2\) and the projected vector itself shortest. A strict inequality would omit those lattices. The note uses \(\ge\), and both programs keep equality. The builder’s diagonal threshold is \(\mathrm{ceil}(q+(3/4)d_{\mathrm{prev}})\), which preserves an integral endpoint. The scanner rejects a new pivot only when it is strictly less than \((3/4)\) times the previous pivot.

Checked examples inside the determinant-24 domain:

- Both Gram matrices of \(I_2\perp A_2\), with lower-right block \(\begin{pmatrix}2&\pm 1\\ \pm 1&2\end{pmatrix}\), occur. Their pivots are \(1,1,2,3/2\), and \(3/2=(3/4)\cdot 2\). The cross entry \(\pm 1\) is an endpoint of the interval of radius \(1\) about \(0\).
- One admitted basis of \(A_4\) is
  \[
  \begin{pmatrix}2&-1&-1&-1\\ -1&2&0&0\\ -1&0&2&1\\ -1&0&1&2\end{pmatrix},
  \]
  with pivots \(2,3/2,4/3,5/4\). Thus \(3/2=(3/4)\cdot 2\) and \(5/4=(3/4)\cdot(4/3)\). Every \(|\mu|\) is at most \(1/2\). The matrix \(U\) with columns \((-1,-1,-1,-1)\), \((0,0,0,1)\), \((1,0,0,0)\), \((1,1,0,0)\) satisfies \(U^tG_{\mathrm{Cartan}}U\) equal to this Gram. The Cartan basis itself has a coefficient \(3/4\) and is outside the admitted inequalities.
- One admitted basis of \(D_4\) is
  \[
  \begin{pmatrix}2&-1&-1&-1\\ -1&2&0&0\\ -1&0&2&0\\ -1&0&0&2\end{pmatrix},
  \]
  with pivots \(2,3/2,4/3,1\) and \(1=(3/4)\cdot(4/3)\). The control Gram of \(D_4\) has a coefficient \(2/3\) and is outside the admitted inequalities; the lattice is present through this equivalent basis. The same pattern holds for \(A_3\) perpendicular to a norm-one line: an admitted basis begins with the norm-one vector and has pivots \(1,2,3/2,4/3\).

Across all 1,510 bases, 1,402 have some \(|\mu|=1/2\) and 988 have some adjacent pivot equality. The closed bounds are the generic case of this domain, not an empty endpoint.

### 1.3. Determinant bounds and cross entries

Deduction. The shortest squared length \(A=d_1\) is a positive integer. Chaining the pivot inequality gives

\[
\Delta=d_1d_2d_3d_4\ge\Bigl(\frac34\Bigr)^{0+1+2+3}A^4=\Bigl(\frac34\Bigr)^6 A^4=\frac{729}{4096}A^4.
\]

For \(\Delta\le 24\),

\[
\frac{729}{4096}\cdot 81=\frac{59049}{4096}<24,\qquad \frac{729}{4096}\cdot 256=\frac{729}{16}>24.
\]

So \(A\le 3\), and \(A=3\) still satisfies the inequality. Both programs use the closed test \((3/4)^6 A^4\le H\). The recomputed domain has 328 bases with first diagonal 1, 1,166 with first diagonal 2, and 16 with first diagonal 3.

For a prefix of length \(k-1\), with \(P=\prod_{j<k}d_j\) and \(r=5-k\) later pivots including \(d_k\),

\[
\Delta\ge P\, d_k^r\Bigl(\frac34\Bigr)^{r(r-1)/2}.
\]

The exponent is \(0+1+\cdots+(r-1)\). This is an upper bound on the next pivot once the right-hand side exceeds \(H\). Equality is kept: the builder’s loop continues while the lower bound is at most \(H\), and the scanner breaks only after the bound has become strictly larger.

Every original basis vector has squared length at least \(A\), because \(v_1\) was shortest, so each later diagonal is at least \(A\). The Gram–Schmidt relation

\[
\mu_{kj}=\bigl(G_{jk}-\sum_{l<j}\mu_{jl}\mu_{kl}d_l\bigr)/d_j
\]

puts \(G_{jk}\) in the closed real interval of radius \(d_j/2\) about that sum. The integral points are exactly the integers from \(\mathrm{ceil}(\mathrm{center}-d_j/2)\) through \(\mathrm{floor}(\mathrm{center}+d_j/2)\). Checked on the half-integral endpoints: the interval of radius \(1\) about \(0\) contributes \(\{-1,0,1\}\), and \(\mathrm{ceil}(1/2+(3/4)\cdot 2)=2\).

The last diagonal is solved from \(P(G_{44}-q)=\Delta\) for every integer \(\Delta\) from 1 through \(H\). An integral value that meets the same pivot and diagonal bounds is kept. Positive pivots are the LDL certificate of positive definiteness. The determinant of an integral Gram equals the product of the pivots, so the solved \(\Delta\) is the lattice determinant. The recomputed domain contains the identity, of determinant 1, and 184 bases of determinant 24.

The scanner’s coarser boxes were checked as outer bounds, not as a second completeness proof. From \(d_{j+1}\ge(3/4)d_j\) and \(d_1\ge 1\),

\[
d_i\le H\Big/\Bigl(\frac34\Bigr)^{7-i}
\]

in the note’s one-based indexing. The scanner’s zero-based list \(H/(3/4)^{6-i}\) is the same sequence. Then

\[
G_{ii}=d_i+\sum_{j<i}\mu_{ij}^2 d_j\le U_i+\frac14\sum_{j<i}U_j,
\]

and \(|G_{jk}|\le\sum_{l\le j}|\mu_{jl}|d_l/2\). An integer upper bound contributes every integer up to its floor, including an integral endpoint. The Schur test on the scanner rejects a new coefficient only when its absolute value is strictly greater than \(1/2\). These boxes may admit still more matrices before the filter; they do not cut away a matrix that meets the closed inequalities. The two enumerations agree for every bound from 1 through 24.

Sixteen admitted bases have first diagonal 3. Each has minimum norm 2, so the first vector of that basis is not shortest: the note’s extra bases occur in this domain. Each of the sixteen is integrally equivalent to an admitted basis whose first diagonal equals that minimum. The \(A=3\) branch is exercised, and those extras do not displace a shortest-vector basis of the same lattices.

A generator node limit raises and discards the partial list. `candidates` does not catch the exception, and the resource-limit test requires the raise. The mathematical inequalities are not restricted to determinant 24. The executable feasibility domain is the interval of bounds implemented and compared here, \(1\le H\le 24\).

## 2. The weight-two cutoff

### 2.1. The common space

Checked source. Kane–Kim, arXiv:2211.03987v2, §2.1, define the discriminant as \(\det G\) and the level as the least positive integer \(N\) for which \(N G^{-1}\) is integral. Proposition 2.3, read in this pass with its proof, places the theta series of an integral coset of conductor \(a\) and rank \(k\) in weight \(k/2\) on \(\Gamma_0(4N_L a^2)\cap\Gamma_1(a)\), with character \(\chi_{4d_L}\) when \(k\) is odd and \(\chi_{(-1)^{k/2}4d_L}\) when \(k\) is even. Identity (2.9) in that proof is the slash relation before the conductor forces the shift to be trivial. Theorem 1.2 of the same paper is the ternary decomposition into genus, spinor genus, and class. The cutoff uses Proposition 2.3.

For rank \(k=4\) the weight is 2 and the character is \(\chi_{4d_L}\). The displayed Shimura calculation in the proof of Proposition 2.3 reduces, for even \(k\), as follows. The exponential phase is an integer, hence 1. The powers of the conductor and of 2 cancel in the Kronecker symbol, leaving \((d_L/s)(r/s)^k\varepsilon_s^{-k}\). For even \(k=2m\), \((r/s)^k=1\) and \(\varepsilon_s^{-k}=((-1)^m/s)\). Since the lower-right entry is odd, \((4/s)=1\). The product is \(\chi_{(-1)^{k/2}4d_L}(s)\). At \(k=4\) this is \(\chi_{4d_L}\). The integral-weight slash is then \((cz+d)^{-2}\).

Checked project note, re-read in this pass. `docs/coset-modularity.md` §§1–3 package the two doubled cosets. For \(b\ne 0\), \(L_0\) has index two and \(d_0=\det G_0=4\det G\). A shift already in \(2L_0\) is the lattice \(2L_0\), with base \(2L_0\) and conductor 1. A shift outside \(2L_0\) has base \(L_0\) and conductor 2. Section 2’s valuation shows that the level of \(2L_0\) is exactly \(4N_0\): if a prime divided every entry of \(H=N_0 G_0^{-1}\), then \(0=4 v_\ell(N_0)=v_\ell(d_0)+v_\ell(\det H)\ge 4\).

Conductor 2 on \(L_0\) gives \(\Gamma_0(16N_0)\cap\Gamma_1(2)\). Conductor 1 on \(2L_0\) gives \(\Gamma_0(4\cdot 4N_0)=\Gamma_0(16N_0)\), with discriminant multiplied by \(4^4\). On this group the lower-right entry \(s\) is odd, and \((4\cdot 4^4 d_0/s)=(4d_0/s)\). Section 3 records that \(\Gamma_1(2)\) is redundant, because a matrix of \(\Gamma_0(16N_0)\) already has odd diagonal entries. Both unsigned summands, and therefore their difference \(F\), lie in

\[
M_2\bigl(\Gamma_0(16N_0),\,\chi_{4d_0}\bigr).
\]

The Fourier index of \(F(z)=\sum c_N\exp(2\pi i Nz)\) is norm4, including \(N=0\). The cusp width of \(\Gamma_0(16N_0)\) at infinity is 1. This is the variable of the cutoff. Section 4 of the same note is the rescaling law, whose group has cusp width 8; that group is not the group in this cutoff.

Checked on the five control kernels and on the cancelled-shell Gram: the character of discriminant \(4\cdot 4^4 d_0\) agrees with \(\chi_{4d_0}\) at every unit modulo \(M=16N_0\).

### 2.2. The square

Deduction. For \(\gamma=\begin{pmatrix}*&*\\ *&s\end{pmatrix}\in\Gamma_0(M)\) the transformation is \(F(\gamma z)=\chi(s)(cz+d)^2 F(z)\). Every prime divisor of \(d_0\) divides \(N_0\), and \(2\) divides \(M\), so \(\gcd(s,4d_0)=1\) whenever \(\gcd(s,M)=1\). Thus \(\chi(s)=\pm 1\) and

\[
F(\gamma z)^2=(cz+d)^4 F(z)^2.
\]

The character of \(F^2\) on \(\Gamma_0(M)\) is trivial, and the weight is 4. Weight 4 is even, so \((-I)\) acts by \(+1\), which agrees with the trivial character at \(s=-1\). Brunault’s §0.1 records that an odd-weight form for a congruence subgroup containing \(-I\) is zero. That vanishing does not apply to \(F^2\).

On \(A_4\), with the default kernel \(b=(1,0,0,0)\), one has \(d_0=20\), character discriminant \(80\), and level \(M=160\). The matrix \(\begin{pmatrix}107&2\\ 160&3\end{pmatrix}\) lies in \(\Gamma_0(160)\), and \((80/3)=-1\). The square of that value is \(1\). The same square check holds at \(-1\) and at every unit modulo \(M\) for \(I_4\), \(2I_4\), \(A_4\), \(D_4\), \(A_3\) plus a line, and the cancelled-shell Gram.

Holomorphy. Each coset theta series lies in Kane–Kim’s space, so it is holomorphic on the upper half-plane and of polynomial growth at the cusps. The same is true of \(F\). The square remains holomorphic on the upper half-plane, and polynomial growth passes to it. At a cusp \(\sigma(\infty)\), trivial character and even weight make

\[
\phi(z)=(c_\sigma z+d_\sigma)^{-4}F(\sigma z)^2
\]

periodic under a positive translation, the width of that cusp. Its Fourier series in the local parameter therefore has only finitely many negative powers, and polynomial growth as the imaginary part tends to infinity forces those coefficients to vanish. Thus \(F^2\) is holomorphic at every cusp. It lies in \(M_4(\Gamma_0(M))\) in Brunault’s sense. No coefficient is asserted to vanish at any cusp.

### 2.3. Brunault’s bound and the inclusive endpoint

Checked source. François Brunault, *Sturm bounds for general congruence subgroups*, 28 May 2021, author’s PDF, all three pages. Theorem 1: if \(\Gamma\) is a congruence subgroup of \(\mathrm{SL}_2(\mathbb Z)\), \(m\) is the index of \(\pm\Gamma\) in \(\mathrm{SL}_2(\mathbb Z)\), and \(f\in M_k(\Gamma)\) satisfies \(\mathrm{ord}_\infty(f)>km/12\), equivalently \(a_n(f)=0\) for \(0\le n\le\lfloor km/12\rfloor\), then \(f=0\). The order is the least index \(n\) with \(a_n\neq 0\) in the expansion \(\sum a_n q^{n/w}\), where \(w\) is the width of the cusp infinity. Section 0.1 defines \(w\) as the least positive integer such that the translation by \(w\) lies in \(\pm\Gamma\).

For \(\Gamma_0(M)\) the translation by 1 lies in the group, so \(w=1\). Also \(-I\in\Gamma_0(M)\), hence \(\pm\Gamma=\Gamma\) and Brunault’s \(m\) equals \(\mu(M)\). The proof’s comparison of \([\mathrm{SL}_2(\mathbb Z):\Gamma]\) with \([\mathrm{SL}_2(\mathbb Z):\pm\Gamma]\) contributes the factor \([\pm\Gamma:\Gamma]=1\) in this case. For even weight the expansion is an ordinary integer power series in \(q\). Here the index \(n\) is the exponent of \(q\), which is the norm4 index of \(F\).

Deduction. Apply Theorem 1 to \(F^2\) in weight 4. Then \(km/12=\mu/3\). Since \(M\) is divisible by 16, the two-primary part of \(\mu\) is \(2^{e-1}\cdot 3\) with \(e\ge 4\), so \(2^{e-1}\) is a multiple of 8 and \(24\) divides \(\mu\). In particular \(6\) divides \(\mu\), \(B=\mu/6\) is an integer, and \(\mu/3=\lfloor 4\mu/12\rfloor\). A nonzero form has order at most \(\mu/3\).

Let \(m\) be the least index with \(c_m\neq 0\). The coefficient of \(q^{2m}\) in \(F^2\) is \(c_m^2\): any other pair of indices summing to \(2m\) has one index strictly smaller than \(m\). In \(\mathbb C\), \(c_m\neq 0\) implies \(c_m^2\neq 0\). Thus \(\mathrm{ord}_\infty(F^2)=2\,\mathrm{ord}_\infty(F)\) whenever \(F\neq 0\), and \(m\le\mu/6\). Vanishing of \(c_N\) for \(0\le N\le B\) forces \(m\ge B+1\), hence

\[
\mathrm{ord}_\infty(F^2)\ge 2(B+1)=\mu/3+2>\mu/3.
\]

Brunault’s theorem gives \(F^2=0\). The ring \(\mathbb C[[q]]\) is an integral domain, so \(F=0\). This is an identity of complex Fourier series.

The endpoint is required by that arithmetic. Vanishing only through \(B-1\) allows order exactly \(\mu/3\), with leading coefficient \(c_B^2\). For \(I_4\) and \(b=(1,0,0,0)\), one has \(N_0=4\), \(M=64\), \(\mu=96\), and \(B=16\). Then \(\mu/3=32\), while \(2\cdot 16=32\) and \(2\cdot 17=34\). The stored control whose bound is one less than the cutoff remains unresolved.

The same integrality was recomputed for every control below: each index is divisible by 24, and each cutoff equals that index divided by 6.

| Lattice, default kernel \(b=(1,0,0,0)\) | \(N_0\) | \(M\) | \(\mu\) | \(B\) | Character discriminant |
|---|---:|---:|---:|---:|---:|
| \(I_4\) | 4 | 64 | 96 | 16 | 16 |
| \(2I_4\) | 8 | 128 | 192 | 32 | 256 |
| \(A_4\) | 10 | 160 | 288 | 48 | 80 |
| \(D_4\) | 4 | 64 | 96 | 16 | 64 |
| \(A_3\) perpendicular to a norm-one line | 16 | 256 | 384 | 64 | 64 |

The note’s four cutoffs 16, 32, 48, and 16 match the first four rows. The fifth row is the same kernel convention on the fifth control.

The cancelled-shell Gram \(\begin{pmatrix}2&-1&-1&0\\ -1&4&-1&0\\ -1&-1&4&0\\ 0&0&0&1\end{pmatrix}\), with \(a=(0,1,1,0)\) and \(b=(1,0,0,0)\), has \(d_0=80\), \(N_0=80\), \(M=1280\), \(\mu=2304\), and \(B=384\). Recomputed shells: the norm4 coefficient at 6 is 0, and the norm4 coefficient at 10 is \(+2\). The run bounded by 6 is unresolved. The run bounded by 10 is proved nonzero by the coefficient 2. A zero at the first occupied shell is not an identity test.

On \(I_4\) with \(a=0\) and \(b=(1,0,0,0)\), the complete shell at norm4 0 has signed coefficient \(+1\). The range starts at 0 because that coefficient can be nonzero. One nonzero complete coefficient proves nonvanishing at once. The bound \(B\) is used in the vanishing direction.

An interrupted enumeration stores an unresolved verdict and omits the shell list. Replay of that certificate does not promote partial cancellation to a proof. The cutoff statement has no determinant cap and does not assert that \(16N_0\) is the minimal level.

## 3. Evidence and implementation

Checked against `code/rank4_feasibility.py`, `code/verify_rank4_feasibility.py`, `tests/test_rank4_feasibility.py`, and the packs in `data/rank4-feasibility/`.

The saved candidate pack has determinant bound 24, generation marked complete, and isometry deduplication marked false. Its 1,510 bases equal `independent_candidates(24)`, which this pass also found equal to `candidates(24)`. The same equality holds for every bound from 1 through 24; that is the first test. The saved sizing rows are the counts 1, 2, 47, 166, 406, 688, and 1,510 at bounds 1, 2, 4, 8, 12, 16, and 24. Those counts enumerate admitted bases. They do not enumerate isometry classes. A shear of \(I_4\), with Gram \(U^t U\) for \(U=I+E_{12}\), has full group order 384 and is absent from the candidate domain: it fails size reduction, while the identity basis of the same lattice is present.

The control file’s domain string is “five named controls; not a census.” Replay covers \(I_4\), \(2I_4\), \(A_4\), \(D_4\), and \(A_3\) plus a line, each with all 256 binary pairs. The five full groups and all 1,280 stabilizer sign lists were recomputed from transported characteristics. Odd pairs replay as symmetry zeros. The even counts match the note:

| Control | Group order | Even symmetry zeros | Even nonzeros | Odd zeros |
|---|---:|---:|---:|---:|
| \(I_4\) | 384 | 55 | 81 | 120 |
| \(2I_4\) | 384 | 55 | 81 | 120 |
| \(A_4\) | 240 | 0 | 136 | 120 |
| \(D_4\) | 1,152 | 9 | 127 | 120 |
| \(A_3\) perpendicular to a norm-one line | 96 | 28 | 108 | 120 |

The independent group search lists every integer vector whose squared length is at most the largest basis diagonal, inside the box \(|x_i|\le\sqrt{\mathrm{bound}\cdot(G^{-1})_{ii}}\). For a positive rational \(x\), \(\lfloor\sqrt{\lfloor x\rfloor}\rfloor=\lfloor\sqrt x\rfloor\), so the integer box contains the sphere. Compatible column pairs are then joined, and a matrix is kept when its determinant is \(\pm 1\) and it reproduces the Gram matrix. A finished run is the full automorphism group. The stored list is required to equal that group. Deleting one stored matrix, one characteristic pair, or one stabilizer sign makes the certificate invalid. A proper subgroup therefore cannot be used to omit an epsilon-one element of the full stabilizer.

\(D_4\) contributes nine even symmetry zeros. One stored witness is \(a=(0,0,1,1)\), \(b=(0,1,0,0)\), with epsilon 1 on the recomputed stabilizer. That witness is a symmetry proof for this pair. The unrestricted rank-four symmetry converse is not a consequence of five controls, and \(P\) if and only if indecomposable is not a consequence of the rank-three classification.

Ten cutoff certificates were replayed. Three remain unresolved: the run whose bound stops one short of \(B\), the interrupted run, and the cancelled first shell at norm4 6. Every cutoff replay in the test reports that the full stabilizer was not checked. A modular zero is an identity certificate for \(F\). It is not a stabilizer computation. The resource-limit test raises on a one-node candidate search and on a one-node group search, marks the one-node control replay unresolved, and keeps the enumeration object out of an interrupted cutoff certificate. Tampered level, index, cutoff, weight, character discriminant, shell list, or bound is rejected. A right shear of the index-two kernel, with the Gram of the kernel updated, still replays as the same modular zero. The index test compares \(\mu(M)\) with the count of bottom rows modulo \(M\) up to units, for the levels arising from \(M\in\{16,32,48,64\}\) in that diagonal family.

Pack hashes read in this pass, SHA-256:

| File | SHA-256 |
|---|---|
| `data/rank4-feasibility/candidates-det24.json` | `bdd2719c3f500d43459e47a8f6e4b6dcf137b6db736bffd453d1edd348f27d5d` |
| `data/rank4-feasibility/candidate-sizing.json` | `0f9f10497b1a7f48c9f92474164d61f43ae7ca2ea1bde4ce588e2a7193ed924d` |
| `data/rank4-feasibility/controls.json` | `7490e6e5a4e9a0e792bf677dcd30068fc59c6f037e66ad7fb8cba0b02e7a5e63` |
| `data/rank4-feasibility/cutoff-controls.json` | `e2107b28667ff08011f8521ab9b850e64f9cdf9993de57da778fe58e4fb9e31a` |

## 4. What this pass ran and read

Ran: `python scripts/run_tests.py --all`, eight tests, all passing in 7.189 seconds; `candidates(24)` against `independent_candidates(24)` and against the saved pack; the pivot, coefficient, and determinant census of those 1,510 bases; explicit isometries for the \(A_4\) and \(D_4\) witness bases above; shortest-basis coverage of the sixteen first-diagonal-3 extras; `modular_data` and Kronecker values for the five controls and the cancelled-shell Gram; `cutoff_build` for the norm4-0 shell of \(I_4\) and for the norm4-6 and norm4-10 shells of the cancelled-shell Gram.

Read: `Grok-rank4-feasibility.md`; `RESEARCH-CONTINUATION-PLAN.md`; `docs/rank-four-census-domain-and-cutoff.md`; `docs/coset-modularity.md` §§1–3; the generation, modular, replay, and test sources named above; Kane–Kim, arXiv:2211.03987v2, §§2.1–2.2 and Proposition 2.3, including identity (2.9) and the statement of Proposition 2.2 quoted there; Brunault’s note, Theorem 1 and §§0.1–0.2, from the author’s PDF.

Unread, and not used as proof: Shimura’s paper beyond the transformation formula quoted in Kane–Kim’s Proposition 2.2; Kane–Kim §§3–6; Brunault’s references Diamond–Shurman, Stein, and Sturm, Lecture Notes in Mathematics 1240 (1987), 275–280. Section 4 of `docs/coset-modularity.md` was not used as the cutoff group.
