# Prior art for the uniform rank-three theorem

Prepared 2026-09-27 from `Grok-rank3-prior-art.md`. This report looks for a predecessor of the signed-characteristic theorem in `docs/rank-three-classification.md`. It does not reopen the shell proof, the cutoff, or the determinant-at-most-24 census. The manuscript and the review bundle were not changed.

The target is the following statement. For every positive definite integral \(3\times 3\) Gram matrix \(G\), and for \(a,b\in\{0,1\}^3\),

\[
\theta[a,b](\tau;G)
=\sum_{x\in\mathbb Z^3}
\exp\bigl(\pi i\tau\,(x+a/2)^tG(x+a/2)
+\pi i\,(x+a/2)^tb\bigr).
\]

Parity is \(p=a^tb\bmod 2\). Property \(P\) means that every even characteristic gives a function of \(\tau\) that is not identically zero. The reviewed theorem says that \(P\) holds exactly for the orthogonally indecomposable lattices, with no determinant restriction, and that every identically zero characteristic, odd or even, has an integral isometry \(T\) with \(T^tGT=G\), \(Ta\equiv a\), \(T^{-t}b\equiv b\pmod 2\), and

\[
\varepsilon=a^t(T^{-t}b-b)/2\equiv 1\pmod 2.
\]

No predecessor of that pair of claims was located in the sources opened below. That is a statement about this search, not a novelty claim.

## 1. Applicable results

Nothing opened here has hypotheses that imply the indecomposable direction, the cycle, diamond, or \(K_4\) signed-shell evaluation, or the symmetry converse for even zeros.

Two classical ingredients do imply the odd-characteristic half of claim 2, and the decomposable failure of \(P\). Both are already the background of `docs/foundations-and-exact-enumeration.md` sections 2 and 7. They are recorded here with the convention translation.

### Odd nullwerte, for every period matrix

Checked source. NIST DLMF 21.2.5 and 21.3.6, read at `https://dlmf.nist.gov/21.2` and `https://dlmf.nist.gov/21.3`. The notes to §21.3(ii) cite Mumford, *Tata Lectures on Theta I* (1983), pp. 120–122. The formulas themselves are:

\[
\theta\begin{bmatrix}\alpha\\\beta\end{bmatrix}(z\mid\Omega)
=\sum_{n\in\mathbb Z^g}
\exp\Bigl(2\pi i\Bigl(\tfrac12(n+\alpha)\cdot\Omega\cdot(n+\alpha)
+(n+\alpha)\cdot(z+\beta)\Bigr)\Bigr),
\]

\[
\theta\begin{bmatrix}\alpha\\\beta\end{bmatrix}(-z\mid\Omega)
=(-1)^{4\alpha\cdot\beta}\,
\theta\begin{bmatrix}\alpha\\\beta\end{bmatrix}(z\mid\Omega),
\]

for half-period characteristics. Set \(z=0\), \(\Omega=\tau G\), \(\alpha=a/2\), and \(\beta=b/2\). The quadratic term is \(\pi i\tau\,(x+a/2)^tG(x+a/2)\). The linear term is \(2\pi i(x+a/2)\cdot(b/2)=\pi i(x+a/2)^tb\). This is the series in the brief. The sign exponent is \(4\alpha\cdot\beta=a^tb=p\).

At \(z=0\), an odd characteristic therefore satisfies \(\theta(0\mid\Omega)=-\theta(0\mid\Omega)\) for every period matrix \(\Omega\). The nullwert is the zero function on Siegel space, and its restriction to the ray \(\Omega=\tau G\) is the zero function of \(\tau\). On the lattice, \(x\mapsto -x-a\) is the map induced by \(T=-I\). It preserves \(G\), fixes every characteristic modulo 2, and has \(\varepsilon\equiv p\pmod 2\). For odd \(p\) this is the witness in claim 2. DLMF states the functional equation in \(z\). It does not mention \(\mathrm{Aut}(L)\) or indecomposability.

The same identity does not constrain an even characteristic. For even \(p\) the sign is \(+1\), and 21.3.6 is compatible with a nonzero nullwert.

### The product on a block-diagonal period matrix

Checked formulas, then a deduction. If \(G=\mathrm{diag}(G_1,G_2)\) in a \(\mathbb Z\)-basis, then \(\Omega=\tau G\) is block diagonal for every \(\tau\). The sum in DLMF 21.2.5 separates, and

\[
\theta\begin{bmatrix}\alpha_1,\alpha_2\\\beta_1,\beta_2\end{bmatrix}(0\mid\tau G)
=\theta\begin{bmatrix}\alpha_1\\\beta_1\end{bmatrix}(0\mid\tau G_1)\,
\theta\begin{bmatrix}\alpha_2\\\beta_2\end{bmatrix}(0\mid\tau G_2).
\]

If both factors are odd, 21.3.6 makes each factor identically zero, so the product is identically zero, while the total parity \(p\) is even. An orthogonal sum of three lines is the same argument with three factors. This is the decomposable failure of \(P\). It uses a splitting basis that already exhibits the orthogonal decomposition. It does not decide the indecomposable lattices, and it does not evaluate a shell.

The unsigned special case \(a=b=0\) is written out by Elkies, lecture notes *Lattices, Linear Codes, and Invariants*, `https://people.math.harvard.edu/~elkies/aws09.pdf`, displayed equation (15):

\[
\Theta_{L_1\oplus L_2}(q)=\Theta_{L_1}(q)\,\Theta_{L_2}(q),
\qquad
\Theta_L(q)=\sum_{v\in L}q^{\langle v,v\rangle/2}.
\]

That identity was read in the notes. It is the characteristic-zero product. It does not by itself produce an even vanishing characteristic.

## 2. Sources read that do not settle the target

### Selling, *Des formes quadratiques binaires et ternaires*

Checked. Édouard Selling, Journal de Mathématiques Pures et Appliquées (3) 3 (1877), 21–60, `https://www.numdam.org/item/JMPA_1877_3_3__21_0.pdf`. Section III, journal p. 43, begins “Formes ternaires définies” and “Nouvelles conditions de réduction.” Formula (5) on journal p. 44 writes a positive ternary form, after introducing a fourth variable \(t\) initially zero, as a sum of six homogeneous coefficients times squared differences of four variables. Adding the same integer to all four variables does not change the value. Formula (6) is the same expression in the six paired coefficients. The last page, journal p. 60, classifies the resulting symmetry by crystal systems (rectangular, rhombohedral, cubic, and the others listed there).

A text search of all 41 pages found no theta series and no half-characteristic. The hits on “theta” and “caractère” are the reduction and crystallographic vocabulary. This is the classical superbase expansion. It does not discuss \(\theta[a,b]\).

The German paper cited by Conway and Sloane, Selling, J. Reine Angew. Math. 77 (1874), 143–229, was not opened. See section 3.

### Conway and Sloane, Low-dimensional lattices VI

Checked. Author’s text at `http://neilsloane.com/doc/fedorov.pdf`, read in full. The file states that a slightly different version appeared in Proc. Roy. Soc. London Ser. A 436 (1992), 55–68. Citations below use the theorem numbers of the author’s text.

Theorem 3 gives Selling’s formula \(N(\sum m_iv_i)=\sum_{i<j}p_{ij}(m_i-m_j)^2\) for an obtuse superbase, with \(p_{ij}=-v_i\cdot v_j\), and identifies the minimal vectors in each class of \(L/2L\). Theorem 8 states that every 3-dimensional lattice has an obtuse superbase. The proof is the adjacent-superbase algorithm in §7: if a putative conorm equals \(-\varepsilon\) with \(\varepsilon>0\), one partial-sum vonorm drops by \(4\varepsilon\), and the process stops because those norms are positive and decrease. Theorem 9 and Table I list Fedorov’s five Voronoi cells by the positions of the zero conorms. Table I names the types “indecomposable tertiary” (rhombic dodecahedron, canonical example \(A_3\)), “decomposable tertiary or simply decomposable” (hexagonal prism, canonical example \(A_2I_1\)), and “quaternary or fully decomposable” (cuboid, \(I_3\)). Section 8.0.2 defines the Delone diagram: \(n+1\) nodes, with an edge labelled by the Selling parameter \(p_{ij}\) whenever \(p_{ij}\neq 0\), attributed to Delone, Uspekhi Mat. Nauk 3–4 (1937–1938), Fig. 37, p. 138.

The word “decomposable” here means an orthogonal splitting of the lattice, read off the zero pattern of that diagram. The diagram is the positive-edge graph used in `docs/rank-three-classification.md`. The paper contains no theta series, no characteristic, and no signed coefficient. The word “theta” does not occur. Theorem 8 supplies existence of the superbase. It does not evaluate shells.

### Kurlin, and Conway’s third lecture

Checked in the previous review, and not re-audited as a proof. Kurlin, arXiv:2201.10543, Theorem 2.8 and Appendix A, and Conway, *The Sensual (Quadratic) Form*, third lecture, printed pp. 69–76, give the same obtuse-superbase framework. `reviews/rank-three-proof-review.md` records that those portions contain no signed characteristic theorem. This pass did not find a later section of either source, within the pages already read, that adds one.

### Kane and Kim

Checked. Ben Kane and Daejun Kim, arXiv:2211.03987v2, `https://arxiv.org/html/2211.03987v2`. The abstract, §1 through Theorem 1.2, and §2 through Proposition 2.3 were read.

The coset series is unsigned:

\[
\Theta_{aL+\nu}(z)=\sum_{x\in aL+\nu}q^{Q(x)}=\sum_{n\ge 0}r(n,aL+\nu)\,q^n,
\qquad q=e^{2\pi iz}.
\]

Theorem 1.2 identifies the Eisenstein, unary-theta, and orthogonal-cuspidal pieces of this series with the proper genus, the proper spinor genus, and the proper class. Proposition 2.3 places \(\Theta_{aL+\nu}\) in weight \(k/2\) on \(\Gamma_0(4N_La^2)\cap\Gamma_1(a)\). The paper studies one unsigned coset at a time, and averages of such cosets. It does not state a criterion for two particular cosets inside one ternary lattice to have identical representation numbers, which is the form the signed difference takes in `docs/foundations-and-exact-enumeration.md` section 8. Indecomposability is not a hypothesis of Theorem 1.2.

### Schiemann’s audibility theorem, as stated in sources that were opened

The primary article, Alexander Schiemann, Math. Ann. 308 (1997), 507–517, was not opened; the publisher page is access-restricted. Two accounts were read.

Georg Hein, arXiv:1106.4895v1, first section, defines

\[
\Theta_A(z)=\sum_{\lambda\in\mathbb Z^n}\exp\bigl(2\pi i\,(\lambda^tA\lambda)\,z\bigr)
\]

and states that if the rank is at most three, then \(\Theta_A=\Theta_B\) if and only if the positive definite forms are equivalent. A thesis exposition, `https://gupea.ub.gu.se/server/api/core/bitstreams/27343d34-fb1f-4b4e-b381-70451449a542/content`, Chapter 5, Theorem 5.0.1, states the same theorem as equality of representation numbers \(R(q,t)=\#\{x\in\mathbb Z^n:q(x)=t\}\).

Both statements are about the unsigned theta series of the whole lattice. Equal representation numbers of two ternary lattices are not the equality of the two signed cosets \(\xi+L_0\) and \(\xi+t+L_0\) inside one lattice. The opened statements give no bridge from audibility to property \(P\).

### Mumford scan, transformation law

Checked pages of the local scan in `sources/pages`. Printed pp. 4–9 introduce the one-variable theta function, the heat equation, and the Heisenberg group. Printed pp. 132–133 discuss complex tori. Printed pp. 190–195 develop the characteristic series and the symplectic transformation. Equation (5.3') on p. 192 is the series with rational characteristics. Proposition 5.5 and Cases I–III on pp. 193–195 give the transformation under \(\mathrm{Sp}(2g,\mathbb Z)\), including the pure basis change \(\vartheta(Az,A\Omega A^t)=\vartheta(z,\Omega)\) up to an eighth root of unity when \(\det A=\pm 1\). Printed pp. 200–201 place half-integral-weight theta series in the Weil representation. Printed pp. 232–233 concern pluriharmonic polynomials.

These pages are transformation laws and definitions. They do not classify the ray \(\Omega=\tau G\) for a ternary Gram matrix, and they do not contain the cycle or diamond shell. The DLMF normalization above is the one matched term-by-term to the project’s series. Mumford’s printed pp. 120–122, the pages DLMF cites for §21.3(ii), were not the spreads opened in this pass.

### Cohen–Zagier, and the project’s MathOverflow question

Checked only for scope, as in the previous review. Cohen and Zagier, *Vanishing and non-vanishing theta values*, study \(\Theta(\chi)=\sum n^\epsilon\chi(n)e^{-\pi n^2/N}\) at one point of the upper half-plane. That is a Dirichlet theta value, not a lattice characteristic along \(\tau G\).

MathOverflow question 515148 asks whether even parity implies that the lattice series is not identically zero, and records that every decomposable lattice is a counterexample. It is the question this theorem answers. The answers visible on the page do not state the rank-three classification.

## 3. Unread or inaccessible leads

These were not used as evidence.

- Selling, J. Reine Angew. Math. 77 (1874), 143–229. This is the German paper named by Conway and Sloane. The 1877 French text above was the Selling file actually read. De Gruyter and EUDML did not return the Crelle pages.
- B. N. Delone, “Geometry of positive quadratic forms,” Uspekhi Mat. Nauk 3 (1937), 16–62, and 4 (1938), 102–164, especially Fig. 37 on p. 138. Known through Conway and Sloane §8.0.2. The Russian article was not opened.
- G. Voronoi, J. Reine Angew. Math. 133, 134, and 136 (1908–1909), beyond the citations in Conway and Sloane.
- Jun-ichi Igusa, *Theta Functions*, Springer, 1972. Named by DLMF as a source for Chapter 21. No copy was opened.
- Mumford, *Tata Lectures on Theta I*, printed pp. 120–122, the specific citation in DLMF §21.3(ii). Nearby scanned pages were opened, as listed above.
- Alexander Schiemann, Math. Ann. 308 (1997), 507–517, and the Bonn thesis, Bonner Mathematische Schriften 268 (1994). The theorem was read only in the two expositions above.
- Nebe, Rains, and Sloane, *Self-Dual Codes and Invariant Theory*, Springer, 2006, Chapter 9. Still not opened. The earlier negative report about p. 261 remains an unchecked lead.
- A systematic MathSciNet or Zentralblatt pass. The searches below are web and archive searches, plus the files named in the brief.

Search terms used, in combination: “theta characteristic identically zero indecomposable ternary”; “half-integral characteristic lattice theta”; “reducible principally polarized theta constant”; “Selling conorm theta”; “Conway Sloane Voronoi theta”; “Schiemann theta series ternary characteristic”; “Kane Kim coset theta difference”; “orthogonal sum theta series product characteristic.”

Coverage gap. The geometric reduction literature and the standard theta-constant identities were opened far enough to separate them from the target. A book-length reading of Igusa, of Mumford’s chapters on the parity of theta characteristics, or of Delone’s 1937–38 memoir could still contain a remark that was not on the pages read. No such remark was found in the portions that were opened, and this search does not close that possibility.

## 4. Attribution language

The following division matches what was read.

Established background, safe to cite for the framework and for the odd and decomposable directions:

- The series is the Riemann theta constant \(\theta[a/2,\,b/2](0\mid\tau G)\) in the normalization of DLMF 21.2.5. For odd \(a^tb\), DLMF 21.3.6 implies that this nullwert vanishes for every period matrix, hence as a function of \(\tau\). The lattice witness is \(T=-I\).
- If the Gram matrix is block diagonal over \(\mathbb Z\), DLMF 21.2.5 factors the series. Two or more odd factors give an even characteristic whose series is identically zero. Elkies’s equation (15) is the unsigned case of the same product.
- Obtuse superbases, Selling’s formula, and the four-vertex Delone diagram are classical. Cite Selling, J. Math. Pures Appl. (3) 3 (1877), formulas (5)–(6), for the four-variable expansion, and Conway–Sloane, Low-dimensional lattices VI, Theorems 3, 8, and 9 and §8.0.2, for the three-dimensional existence theorem and the diagram. Kurlin’s Theorem 2.8 is a modern existence proof with nonnegative conorms. Conway–Sloane’s Table I already separates orthogonally decomposable ternary lattices from the indecomposable conorm types. That separation is geometric.

Present argument, for which this search found no predecessor:

- The converse, that every even characteristic of an orthogonally indecomposable positive definite integral ternary lattice is not identically zero, at every determinant.
- The signed minimal-shell evaluation on the cycle, the diamond, and \(K_4\), including the diamond case in which one coordinate slice cancels and the next shell has complete coefficient of magnitude 2.
- The symmetry converse for the even zeros: each such zero comes from an odd factor, and the factor’s central inversion extends to an epsilon-one isometry of the rank-three lattice.

A manuscript sentence can say that the odd nullwert and the product on an orthogonal sum are classical, that the obtuse-superbase reduction is classical, and that the indecomposable signed-shell classification and the even symmetry converse are proved in `docs/rank-three-classification.md`. It should not say that the search establishes those last two statements as new.

No contradiction with the reviewed theorem appeared. The sources that use the words “indecomposable” and “decomposable” use them for the orthogonal geometry of the lattice or for a principally polarized abelian variety, not for identical vanishing of \(\theta[a,b](\tau;G)\).
