# Adversarial review of the uniform rank-three graph-shell proof

Prepared 2026-09-27 from `Grok-rank3-structure.md`. This note reviews `docs/rank-three-classification.md` as a proposed theorem for every positive definite integral lattice of rank three. The manuscript, cutoff report, census packs, and structure certificates were not changed.

**Checked** means a source read in this pass, or an exact integer, shell, or certificate computation run in this pass. **Deduction** means a step proved from those statements and from `docs/foundations-and-exact-enumeration.md` sections 1–4 and 7 and `docs/rank-two-classification.md`, which this pass re-read and treats as established.

## Verdict

I accept the uniform theorem.

For every positive definite integral lattice \(L\) of rank three, the parity property \(P(L)\) holds if and only if \(L\) is orthogonally indecomposable. Every identically zero half-integral characteristic has a stabilizer element of epsilon \(1\). Vanishing means identical vanishing as a function of \(\tau\). One nonzero complete norm4 coefficient is enough to prove nonvanishing.

The six steps in the brief are correct, including saturation, termination, the articulation direction, the sign normalization, ties between shells, and cancellation on the later diamond shell. The argument is uniform in the determinant. The determinant-\(\le 24\) certificates are a regression of that argument.

Four justifications are asserted in the note and are written out below: a component sum is nonzero, an articulation vertex produces a block-diagonal Gram matrix, antipodal signs agree exactly when \(p\) is even, and the norm-zero class has coefficient \(\pm 1\). The symmetry converse also needs one explicit reminder that odd \(p\) is witnessed by \(-I\). These are insertions into a correct argument. They do not change the statement, and they do not restrict the determinant, the conorms, or the scaling.

## 1. Orthogonal splitting

Checked against the note and against `split_basis`, `decomposition`, and `check_split`.

Let \(d=v^tGv\). If \(v\) is a column and \(d\) divides every entry of \(Gv\), the row \(\lambda=(Gv)^t/d\) is integral and \(\lambda v=1\). For every lattice vector \(x\),

\[
x=(\lambda x)\,v+\bigl(x-(\lambda x)\,v\bigr),
\]

the first summand lies in \(\mathbb Z v\), and the two summands are orthogonal because \(Gv=d\lambda^t\). Since \(\lambda v=1\), the summand \(\mathbb Z v\) is saturated and \(v\) is primitive. Conversely, a primitive generator of a rank-one orthogonal summand meets the same divisibility because its pairings with \(L\) are exactly \(d\mathbb Z\).

In rank three the only orthogonal splittings are \(1+2\) and \(1+1+1\). Every decomposable lattice therefore has a rank-one summand. If \(e_i=m_iv+w_i\) in that splitting, some \(m_i\) is nonzero because the basis spans \(L\), and \(G_{ii}=m_i^2d+\lVert w_i\rVert^2\ge d\). Every splitting generator satisfies \(d\le\max_i G_{ii}\) in the supplied basis. An enumeration of the full sphere through that bound proves both decomposability and indecomposability. An interrupted sphere proves neither.

Extending a primitive \(v\) to a matrix in \(\mathrm{GL}(3,\mathbb Z)\) and shearing the other columns by \(\lambda\) produces a unimodular basis whose Gram matrix is \(\mathrm{diag}(d,H_2)\). Binary reduction of \(H_2\), as in `docs/rank-two-classification.md`, separates three lines (\(B=0\)) from a line plus an indecomposable binary lattice (\(B>0\)).

Checked computations. `complete_primitive` recovered unimodular completions of \((1,0,0)\), \((-1,0,0)\), \((0,-1,0)\), \((0,0,1)\), \((2,3,5)\), and \((-7,11,13)\), and rejected \((0,0,0)\) and \((2,4,6)\). On \(I_3\), on \(\mathrm{diag}(1,A_2)\), and on the skewed three-line Gram \([[2,1,0],[1,1,0],[0,0,1]]\), the sphere found a saturated block basis of the expected type. Each obtuse superbase of those three decomposable lattices has an articulation vertex.

## 2. Selling reduction and the norm identity

Deduction, re-expanded in this pass. The move

\[
v_i'=-v_i,\qquad v_j'=v_j,\qquad v_k'=v_k+v_i,\qquad v_l'=v_l+v_i
\]

preserves the relation \(v_0+v_1+v_2+v_3=0\). Omitting \(v_j\), the new triple is the old triple times

\[
\begin{pmatrix}-1&1&1\\0&1&0\\0&0&1\end{pmatrix},
\]

of determinant \(-1\). One unimodular triple and the sum-zero relation imply that every triple is a \(\mathbb Z\)-basis.

Let \(S\) be the sum of the four squared lengths and let \(\varepsilon=v_i\cdot v_j\). The new sum equals \(S+2\lVert v_i\rVert^2+2v_i\cdot(v_k+v_l)\). The sum-zero relation gives \(v_k+v_l=-v_i-v_j\), so the added term is \(-2\varepsilon\). For an integral lattice a positive inner product is an integer at least \(1\), and \(S\) is a positive integer, so every reducing step drops \(S\) by at least \(2\). The process reaches a superbase with all conorms \(w_{ij}=-v_i\cdot v_j\ge 0\). Any acute pair may be chosen; the descent does not depend on the order.

The same move appears in Kurlin’s Lemma A.1. Conway’s printed pp. 74–76 and Kurlin’s Appendix A track a different quantity: one partial-sum vonorm drops by \(4\varepsilon\). The factor \(2\) in the note is the change in the sum of the four squared lengths. Replacing it by \(4\) would misstate this Lyapunov function. The verifier recomputes \(S_{\mathrm{new}}=S-2\varepsilon\) on every stored step.

Selling’s formula is the expansion

\[
\Bigl\lVert\sum_i x_iv_i\Bigr\rVert^2=\sum_{i<j}w_{ij}(x_i-x_j)^2.
\]

It uses \(\lVert v_i\rVert^2=\sum_{j\neq i}w_{ij}\), which is \(v_i=-\sum_{j\neq i}v_j\). In this note the squared length is the Gram value \(y^tGy\). That integer is the project’s norm4. Fixing any root coordinate \(x_r=0\) represents every lattice vector exactly once, because the other three vectors are a basis.

Integrality is used here and again to guarantee \(h\ge 1\) for a positive conorm. Positive definiteness is used to conclude that a vector of norm \(0\) is the zero vector.

## 3. The positive-edge graph

Deduction. Join \(i\) to \(j\) when \(w_{ij}>0\).

The graph is connected. If it split into components, an integer vector constant on each component and nonconstant overall would have norm \(0\). It remains only to see that the corresponding lattice vector \(y\) is nonzero. A nonempty proper subset of a superbase cannot sum to zero: a single vector is nonzero, the sum of three vectors is the negative of the remaining basis vector, and a two-element subset summing to zero forces the complementary pair to sum to zero as well, so those two vectors are negatives and some three are linearly dependent. Thus \(y\neq 0\), which is impossible. This is the missing line in the connectedness sentence. The claim itself is right, and it holds whether or not the lattice is indecomposable.

Suppose \(r\) is an articulation vertex. The three remaining superbase vectors are a basis of \(L\), and there is no positive edge between distinct components of the graph with \(r\) deleted. Distinct components are therefore orthogonal. Their Gram matrix, in that basis, is block diagonal with nonempty blocks, so \(L\) is an orthogonal direct sum of positive-rank saturated summands. An indecomposable lattice therefore has an obtuse-superbase graph with no articulation vertex.

A connected graph on four vertices with no articulation vertex has minimum degree \(2\). Three or fewer edges cannot meet that degree bound. Four edges force the \(2\)-regular graph \(C_4\). Five edges are \(K_4\) minus one edge: the complement is a single edge, and a degree-\(4\) vertex would leave a degree-\(1\) vertex among the other three. Six edges are \(K_4\).

The shell analysis assumes indecomposability only through this graph. A decomposable lattice whose obtuse superbase had no articulation vertex would satisfy \(P(L)\) by the same shell argument, contradicting the product obstruction in `docs/foundations-and-exact-enumeration.md` section 7. The note’s final paragraph records that consequence in the right direction. The three decomposable controls above have articulation vertices, as this predicts.

## 4. Signs, parity, and the lower bound

Deduction. Coordinates \(a,b\) and the signed coefficient are those of `docs/foundations-and-exact-enumeration.md` section 4:

\[
c_N=\sum_{\substack{y\equiv a\pmod 2\\ y^tGy=N}}(-1)^{((y-a)/2)^tb}.
\]

A superbase dual vector \(\beta_i=\langle v_i,b\rangle\) satisfies \(\sum\beta_i=0\). With root \(r\) normalized so that \(a_r=x_r=0\), the superbase exponent \(\sum\beta_i(x_i-a_i)/2\) equals the engine exponent. Changing root by subtracting \(a_r\) and \(x_r\) preserves the exponent because \(\sum\beta_i=0\). An even shift of the representative \(a\) multiplies every coefficient by one global sign \((-1)^{\langle\ell,b\rangle}\). Adding \(1\) to every coordinate does not change the lattice vector. The parity \(p\equiv\sum a_i\beta_i\pmod 2\) agrees with \(a^tb\) and is unchanged by complementation.

**Antipodal signs.** For a vector \(y\) with root-fixed coordinates \(x\), the engine sign of \(-y\) and the engine sign of \(y\) differ by \((-1)^{\sum\beta_ix_i}\). Since \(x_i\equiv a_i\pmod 2\), the exponent \(\sum\beta_ix_i\) has the same parity as \(p\). Even \(p\) makes the two signs equal. Odd \(p\) makes them opposite. The note asserts this and then uses it in every shell; the calculation above is the missing writeup. `docs/foundations-and-exact-enumeration.md` section 2 already gives the global form of the odd case: \(-I\) lies in the stabilizer and \(\varepsilon(-I)\equiv p\pmod 2\).

**Lower bound.** Inside a parity class, a cross difference is odd and a within-part difference is even. Hence

\[
\lVert y\rVert^2\ge m(S)=\sum_{i\in S,\,j\in T}w_{ij},
\]

with equality if and only if every positive within-part edge has difference \(0\) and every positive cross edge has difference \(\pm 1\). The empty set \(S\) forces every coordinate equal to the root, so the only vector is \(0\). Its signed coefficient is \(+1\) for the representative \(a=0\), and \(\pm 1\) after an even change of representative. Those eight even pairs, one for each \(b\), are nonzero. Section 4 should say this; section 3 already identifies the vector.

## 5. The four indecomposable shells

Deduction, then checked on complete shells. Complement \(S\) so that \(|S|\le 2\). No articulation vertex gives minimum degree \(2\).

**One odd vertex.** Its complement is connected, or else that vertex would be an articulation point. A root in the complement forces the three even coordinates to vanish. The odd vertex has degree at least \(2\), hence a neighbor at coordinate \(0\), so its coordinate is \(\pm 1\). Any larger odd value makes a cross difference of size at least \(3\). The shell is one antipodal pair. Even \(p\) gives coefficient \(\pm 2\).

**Both within-part edges positive.** Each part is constant. A root in \(T\) puts both even coordinates at \(0\). Each odd vertex has degree at least \(2\) and already uses its within-part edge, so each has a cross edge to \(0\). The common odd value is \(\pm 1\). One antipodal pair, coefficient \(\pm 2\).

**Both within-part edges zero.** Each vertex has degree at least \(2\) and no within-part edge, so every cross conorm is positive. The graph is \(C_4\). With root \(r\in T=\{r,t\}\), each odd coordinate is at distance \(1\) from both \(0\) and \(x_t\). Thus \(x_t\in\{0,\pm 2\}\). The value \(x_t=0\) gives four vectors, with the odd coordinates independent signs. The value \(x_t=2\) forces both odd coordinates to \(+1\), and \(x_t=-2\) forces both to \(-1\). These six vectors are three antipodal pairs, all of norm \(m(S)\), for arbitrary positive cross conorms. Even \(p\) makes each pair contribute \(\pm 2\). The sum of three such terms lies in \(\{\pm 6,\pm 2\}\).

**One within-part edge zero.** Complementation puts that edge inside \(S\) and the positive edge \(h=w_{rt}\) inside \(T\). A missing cross edge would drop the corresponding vertex of \(S\) to degree at most \(1\). All four cross conorms are positive, and the graph is the diamond. At norm \(m(S)\) one has \(x_r=x_t=0\) and the odd coordinates independently \(\pm 1\): four vectors.

Even \(p\) means the two values \(\beta\) on \(S\) have the same parity, because \(p\) is their sum once the root representative vanishes on \(T\). If both are even, flipping either odd coordinate from \(+1\) to \(-1\) preserves the sign, so all four signs agree and the coefficient is \(\pm 4\).

If both are odd, the same flip reverses the sign: the exponent changes by \(\beta_p x_p\), a product of two odd integers, and an odd coordinate is not fixed by negation. The whole slice \(x_t=0\) is partitioned into cancelling pairs at every norm. This is the identity that has to be proved. Agreement of the two antipodal pairs inside the minimal shell is compatible with those pairs having opposite signs; their contributions \(\pm 2\) and \(\mp 2\) sum to \(0\). The stronger parity statement, that both values of \(\beta\) are even or both are odd, separates coefficient \(\pm 4\) from this cancellation.

For \(x_t=2k\), the cross terms contribute at least \(m(S)\) and the edge \(h\) contributes \(h(2k)^2\). Thus \(|k|\ge 2\) gives norm at least \(m+16h\). Since \(h\ge 1\), this is strictly above \(m+4h\). For \(k=\pm 1\) the norm \(m+4h\) occurs exactly when every cross difference has size \(1\), hence at exactly one vector for each sign of \(k\): both odd coordinates equal \(+1\), or both equal \(-1\). Those vectors are \(\pm(v_t-v_r)\). Even \(p\) gives them equal signs. Directly, their exponents are \(\beta_t\) and \(\beta_r\), which agree because both values of \(\beta\) on \(S\) are odd and the four values of \(\beta\) sum to \(0\).

Any vectors of the slice \(x_t=0\) that happen to have norm \(m+4h\) cancel among themselves. The slices \(|k|\ge 2\) lie strictly higher. The complete coefficient at \(m+4h\) is therefore \(\pm 2\).

The support at \(m+4h\) may be larger than those two vectors. A \(k=0\) vector with an odd coordinate \(\pm 3\) raises the cross contribution by \(8\) times the sum of the two positive cross conorms at that vertex, hence by at least \(16\). Whenever that excess equals \(4h\), those cancelling vectors lie on the same shell. The coefficient is still \(\pm 2\). The certificate field `selected_vectors` lists the two vectors that determine the sign. The independent replay compares complete signed sums. Requiring the stored list to exhaust the shell would reject valid certificates.

Checked shells, each compared three ways: `shell_claim`, a direct walk in root-fixed superbase coordinates, and `box_shells`.

| Lattice | Later shell |
|---|---|
| Diamond, all cross conorms \(1\), \(h=1\) | Norm \(m+4\), exactly \(2\) vectors, coefficient \(\pm 2\) |
| Same diamond, \(h=4\), \(a=(0,1,1)\), \(b=(0,1,1)\) | Norm \(20\), \(10\) vectors, coefficient \(+2\) |
| Same diamond, \(h=8\), same \(a,b\) | Norm \(36\), \(6\) vectors, coefficient \(+2\) |
| \(h=4\) after the basis change \([[1,3,-2],[0,1,5],[0,0,1]]\) | Norm \(20\), \(10\) vectors, coefficient \(+2\) |
| Census class \(97\), \(G=[[2,0,-1],[0,2,-1],[-1,-1,6]]\), determinant \(20\) | Norm \(4\): \(4\) vectors, coefficient \(0\). Norm \(20\): \(10\) vectors, coefficients \(+2\) and \(-2\) for \(b=(1,1,0)\) and \(b=(1,1,1)\) |
| Stated Gram \([[2,-1,-1],[-1,4,-1],[-1,-1,4]]\), \(a=(0,1,1)\), \(b=(1,0,0)\) | Norm \(6\): \(4\) vectors, coefficient \(0\). Norm \(10\): \(2\) vectors, coefficient \(+2\). Here \(h=1\) |

The same three-way comparison covered all \(36\) even pairs on equal and unequal \(C_4\), on \(K_4\), on the unequal diamond, on the determinant-\(20\) Gram, and on the determinant-\(29\) Gram from `data/ternary.json`. Coefficients were \(+1\) on the norm-zero class, \(\pm 2\) on singletons and on connected pairs, \(\pm 4\) on the uncancelled diamond minimum, and \(\pm 2\) or \(\pm 6\) on every four-cycle. Every stored coefficient was nonzero and matched the complete shell. No proper subset sum of a superbase was zero.

On \(C_4\) the six minimal vectors really do lie on one shell when the four cross conorms are unequal. The three pair-sums cannot cancel.

## 6. Decomposable lattices and the symmetry converse

Deduction from `docs/foundations-and-exact-enumeration.md` section 7 and `docs/rank-two-classification.md`. The series of an orthogonal sum is the product of the series. If \(c_n\) and \(d_m\) are the least nonzero coefficients of two factors, the product has coefficient \(c_nd_m\) at degree \(n+m\). The same holds for three factors. A factor that is identically zero makes the product identically zero.

Rank one has three even characteristics, with complete coefficients \(1\), \(1\), and \(2\), and one odd characteristic. An indecomposable binary lattice has all ten even characteristics nonzero and has six odd characteristics. A decomposable binary lattice has exactly one even zero, at \(a=b=(1,1)\), with witness \(\mathrm{diag}(-1,1)\).

A sum of three lines therefore has \(3^3=27\) even characteristics with every factor even, all nonzero, and

\[
\binom{3}{2}\cdot 1\cdot 1\cdot 3=9
\]

even characteristics with exactly two odd factors, all zero. A line plus an indecomposable binary lattice has \(3\cdot 10=30\) even nonzero characteristics and \(1\cdot 6=6\) even zeros, the pairs in which both factors are odd. These are the counts in the note. Across the saved \(120\) classes they reproduce the census: \(51\cdot 9+44\cdot 6=723\) even zeros, and the remaining \(25\cdot 36=900\) even characteristics are the indecomposable shells.

Every such zero has an odd factor. Negating that factor fixes every characteristic class modulo \(2\), preserves a block-diagonal Gram matrix, and has epsilon \(1\) because the negated block has odd parity. In coordinates, if \(B\) is the saturated block basis and \(D\) is the corresponding diagonal sign matrix, \(T=BDB^{-1}\) is the witness on the original basis. Extending a binary witness by the identity is the special case in which the odd block is one line of a decomposable binary factor; \(\mathrm{diag}(-1,1)\) is that negation. Odd characteristics in every rank, decomposable or not, are witnessed by \(-I\). The indecomposable lattices contribute no even zeros, so their only zeros are these odd classes.

Checked on \(I_3\), on a line plus \(A_2\), and on the skewed three-line basis: the builder’s witnesses had epsilon \(1\), and the zero counts were \(9\), \(6\), and \(9\). The full structural replay also requires epsilon \(1\) for every stored zero witness.

The symmetry converse is the combination of the three preceding paragraphs. One sentence in the consequences section should cite \(\varepsilon(-I)\equiv p\pmod 2\) from `docs/foundations-and-exact-enumeration.md` section 2 so that the odd classes are visibly included.

## 7. Finite controls

Checked. `python scripts/run_tests.py --all` passed all six tests in this pass. The full replay reports \(7{,}680\) characteristics, \(4{,}320\) even pairs, \(51\) three-line lattices, \(44\) line-plus-indecomposable-binary lattices, \(25\) indecomposable lattices, graph counts \(12\) cycles, \(11\) diamonds and \(2\) copies of \(K_4\), and \(22\) diamond-later shells. The structure pack hash is

```
ae8288b4f1ae0bbddf3f6161200ea9e7e7c52eab28df772146b84aafb1a762f6
```

The replay binds the census pack hash recorded in `docs/rank-three-classification.md`. It recomputes spheres, unimodular block bases, Selling steps, conorms, coverage, witnesses, and complete box coefficients. It does not call the structural builder. `shell_claim` is used by the builder and by the boundary tests; the box comparison is separate.

These controls show that the saved domain matches the theorem. They are not the proof for arbitrary determinant. The determinant-\(29\) Gram and the scaled and unequal conorm graphs above are outside that domain and were checked directly.

## 8. Primary sources

Checked pages and their scope.

**Conway, *The Sensual (Quadratic) Form*, third lecture, printed pp. 69–76.** The ranick PDF `https://webhomes.maths.ed.ac.uk/~v1ranick/papers/conwaysens.pdf` was read at those pages. Conway states Selling’s formula \(N(v)=\sum_{i<j}P_{ij}(m_i-m_j)^2\), defines an obtuse superbase, and gives the adjacent-superbase move used in the note. The existence argument on printed pp. 71–76 deforms a known obtuse superbase and snaps across a wall when one conorm becomes negative. The worked example says the algorithm terminates because a putative vonorm drops. The definition sentence on printed p. 69 uses the strict inequality \(v_i\cdot v_j<0\). The same pages’ conorm diagrams and the Voronoi-cell discussion include the value \(0\). The note’s closed inequality \(w_{ij}\ge 0\) is the one the shell argument requires: the cycle and the diamond are the zero-conorm graphs. These pages contain no half-integral characteristic, no signed shell, and no parity property.

The index entries for theta functions point to the second lecture, printed pp. 36, 45, and 49, which were sampled. They concern the ordinary theta series, audibility, and Schiemann’s theorem that the representation numbers determine a three-dimensional lattice. That is a different statement.

**Kurlin, arXiv:2201.10543, Theorem 2.8 and Appendix A.** Read in full in that portion. Definition 2.5 makes a superbase obtuse when every conorm is nonnegative, and strict when every conorm is positive. Theorem 2.8 says every lattice in \(\mathbb R^3\) has an obtuse superbase. Lemma A.1 is the same vector move as the note, with the conorm update written explicitly, and the termination is the \(4\varepsilon\) drop in one vonorm. The paper’s subject is a continuous isometry classification by root invariants. It contains no theta characteristic.

**Selling (1874), Voronoi, Delone, and Conway–Sloane, *Low-dimensional lattices VI*.** These are the classical sources named by Conway and by Kurlin for superbases, Selling parameters, and Voronoi’s first kind. Kurlin records that the example in Conway–Sloane’s section 7 needed correction and supplies Appendix A in its place. The note’s integer descent for the sum of the four squared lengths is an independent proof for integral lattices. The shell argument does not treat a literature example as a theorem.

**Adjacent theta literature, checked only for scope.** Kane–Kim, arXiv:2211.03987, decomposes unsigned theta series of ternary cosets into genus, spinor, and class pieces. Cohen–Zagier’s vanishing theta values are Dirichlet theta constants at the point of symmetry. Schiemann’s theorem is the audibility result above. The project’s own MathOverflow question 515148 asks the parity converse; it is the question, and it does not contain this classification. I found no statement of the rank-three parity property, and no graph-shell evaluation of these signed coefficients, in the sources opened here. That search does not establish novelty.

## 9. Insertions to make in `docs/rank-three-classification.md`

The statement can stand while these sentences are added.

1. In the connectedness sentence, record that no nonempty proper subset of a superbase sums to zero, by the rank and basis check in section 3 above.
2. At an articulation vertex, say that the complementary triple is a basis and that its Gram matrix is block diagonal.
3. Record the antipodal calculation: the sign ratio of \(y\) and \(-y\) is \((-1)^p\).
4. In section 4, give the empty class coefficient \(\pm 1\), and cite \(\varepsilon(-I)\equiv p\pmod 2\) for every odd class.
5. Replace “degree one” by “degree at most one” in the diamond cross-edge sentence.
6. Record the arithmetic \(27+9\) and \(30+6\) next to the stated zero counts.
7. Say explicitly that \(y^tGy\) is the norm4 index, and that on the diamond the complete coefficient at \(m+4h\) remains \(\pm 2\) when cancelled \(x_t=0\) vectors lie on that shell.

## 10. Unreviewed

Rank four and higher, any larger determinant census, complex-multiplication or algebraic provenance, manuscript packaging, and another modular cutoff were outside this assignment. `MANUSCRIPT.md`, `docs/foundations-and-exact-enumeration.md` beyond the cited sections, and `reviews/rank-three-cutoff-review.md` were not reopened. The ordinary theta series of a ternary lattice may determine its isometry class and still say nothing about these signed characteristics; that distinction was checked only to the extent of the pages named above. Selling’s 1874 paper and Delone’s book were not read in the original.
