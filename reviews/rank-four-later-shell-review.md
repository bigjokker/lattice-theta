# Review of the later-shell lemma and the fibre identity

Prepared 2026-09-27 from `Grok-rank4-later-shell.md`. The census and the Walsh reduction stay accepted. This note checks the uniform coefficient argument in `docs/rank-four-later-shell-lemma.md` and the block identity in `docs/rank-four-decomposition-and-fibres.md`. It does not replay the 91 MB census.

Earlier reports, manuscripts, inputs, inventories, and ZIPs were left unchanged.

**Checked** means a calculation or test run in this pass, or a stated identity re-derived here. **Deduction** means a step obtained from that algebra. The unrestricted rank-four symmetry converse remains open. Astra is not indicated.

## Decisions

**Later-shell lemma: accept.** For every pair of integers \(s,t\ge 3\), the three characteristics in the statement have signed norm4 coefficient 0 at every index \(N<s+t+8\), and the coefficient at \(s+t+8\) is \(-2\), \(+2\), or \(-2\) according to the three columns \(b\). The 21 exact controls are a finite regression, including parameter pairs outside determinant 24. They are not the proof.

**Fibre identity: accept.** Completing the square in a positive definite \(2\times 2\) block gives the stated double sum. Positive definiteness supplies absolute convergence, so the rearrangement is legitimate. The inner shift \((-4/3,-2/3)\) is not a half-integral characteristic of an integral binary lattice, so the rank-at-most-three converse does not apply to that fibre.

**Next task.** Do not look for one cancellation rule that covers every indecomposable corner. On the indecomposable Gram with upper block \(\begin{pmatrix}2&-1\\-1&3\end{pmatrix}\) and lower diagonal \((4,4)\), the same outer labels \(a=(0,0,1,1)\), \(b=(0,1,1,1)\) already have first coefficient \(+2\) at norm4 8. Leading cancellation is special to the \(A_2\) block solved by the lemma. The bounded continuation is to relate a fibre-count that vanishes at every level to an epsilon-one automorphism, using the determinant-15 series (first coefficient \(-2\) at norm 14, trivial epsilon) as the regression that must stay nonzero.

## 1. The coefficient lemma

Let

\[
G(s,t)=\begin{pmatrix}2&-1&-1&-1\\-1&2&0&0\\-1&0&s&0\\-1&0&0&t\end{pmatrix},\qquad a=(0,0,1,1),\qquad m=s+t,
\]

with \(s,t\ge 3\). A supported vector is \(y=(x_0,x_1,u,v)\) with \(x_0,x_1\) even and \(u,v\) odd. Its squared length is

\[
Q(y)=2x_0^2+2x_1^2-2x_0x_1-2x_0(u+v)+su^2+tv^2.
\]

The binary block \(A=\begin{pmatrix}2&-1\\-1&2\end{pmatrix}\) has determinant 3 and inverse \(\frac13\begin{pmatrix}2&1\\1&2\end{pmatrix}\). Completing the square in \(w=(x_0,x_1)\) subtracts \(\frac23(u+v)^2\). Since \(A\) is positive definite,

\[
Q(y)\ge su^2+tv^2-\frac23(u+v)^2.
\]

The right-hand side is at least \(\frac53(u^2+v^2)\) for \(s,t\ge 3\), because \((u+v)^2\le 2(u^2+v^2)\). Thus \(G(s,t)\) is positive definite. The Schur complement also gives the determinant formula in the note:

\[
\det G(s,t)=3st-2s-2t.
\]

### 1.1. Large odd coordinates

For a supported vector, \(|u|\) and \(|v|\) are odd and at least 1, so \(u^2-1\) and \(v^2-1\) are nonnegative. Therefore

\[
Q(y)-m\ge 3(u^2+v^2-2)-\frac43(u^2+v^2)=\frac53(u^2+v^2)-6.
\]

If \(|u|\ge 3\) or \(|v|\ge 3\), then \(u^2+v^2\ge 10\) and \(Q(y)-m\ge 32/3>8\). Every such vector lies strictly above \(m+8\). The sphere \(Q\le m+8\) contains only \(u,v\in\{+1,-1\}\).

### 1.2. The two levels

Put \(x_0=2k\) and \(x_1=2l\). Direct expansion gives

\[
\begin{aligned}
u=v=+1&\Longrightarrow Q=m+8R,&R&=k^2+l^2-kl-k,\\
u=+1,\ v=-1&\Longrightarrow Q=m+8H,&H&=k^2+l^2-kl.
\end{aligned}
\]

The form \(H=(l-k/2)^2+3k^2/4\) is positive definite, hence nonnegative, and zero only at \((0,0)\). The identity \(R=H(k-2/3,l-1/3)-1/3\) gives \(R\ge -1/3\). For integer \(k,l\), \(R\) is an integer, so \(R\ge 0\). Checked: the minimum of \(R\) on \(|k|,|l|\le 20\) is 0. Consequently every supported vector with coordinates \(\pm 1\) has \(Q\ge m\), and \(Q-m\) is a multiple of 8. There is no supported norm strictly between \(m\) and \(m+8\), and none below \(m\).

The antipode \(-y\) stays in the support because \(-1\equiv 1\pmod 2\). If \(y=a+2x\), then \((-y-a)/2=-a-x\). The extra pairing \(-a\cdot b\) is even for each of the three columns \(b\) in the statement, so \(f(-y)=f(y)\). No supported vector is its own antipode. The four sign pairs are therefore two antipodal copies of the charts \(u=v=+1\) and \(u=+1,\,v=-1\).

### 1.3. The solution lists

From \(3k^2/4\le H\le 1\) one gets \(k\in\{-1,0,1\}\). Substitution produces exactly

\[
H=0:\ (0,0);\qquad H=1:\ (1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1).
\]

From \(3k^2/4-k\le R\le 1\) one gets \(k\in\{0,1,2\}\). In particular \(k=-1\) forces \(R\ge 2\). Substitution produces exactly

\[
R=0:\ (0,0),(1,0),(1,1);\qquad R=1:\ (0,1),(0,-1),(2,1).
\]

A direct search over \(|k|,|l|\le 30\) returned these lists and no others. Each pair determines one vector in its chart.

### 1.4. The three sign patterns

For \(b=(0,1,1,1)\) the chart signs are \((-1)^l\) and \(-(-1)^l\). The sums of \((-1)^l\) are \(1,1\) on \(R=0,H=0\) and \(-3,-2\) on \(R=1,H=1\). Two antipodal copies give

\[
c_m=2(1-1)=0,\qquad c_{m+8}=2(-3+2)=-2.
\]

For \(b=(1,0,0,0)\) both charts have sign \((-1)^k\). The sums are \(-1,1\) at level 0 and \(3,-2\) at level 1, so \(c_m=0\) and \(c_{m+8}=2(3-2)=+2\).

For \(b=(1,1,1,1)\) the signs are \((-1)^{k+l}\) and \(-(-1)^{k+l}\). The sums of \((-1)^{k+l}\) are again \(1,1\) and \(-3,-2\), so \(c_m=0\) and \(c_{m+8}=-2\).

These are the complete coefficients through \(m+8\). The value \(\pm 2\) is nonzero, so it is the first nonzero coefficient. An epsilon-one stabilizer element would make every coefficient zero, so these three characteristics have trivial epsilon for every \(s,t\ge 3\). The argument does not need the automorphism group, and it does not identify the group average’s support. Support of the average at norm \(m\), for the two census rows, remains the separate finite scan.

Checked again by an independent double loop, not by the control builder, on \((s,t)=(3,3),\ (3,4),\ (6,9),\ (4,100)\). In each of the twelve series the only occupied norm below \(m+8\) is \(m\), with coefficient 0, and the coefficient at \(m+8\) is the predicted value. No vector with an odd coordinate of absolute value at least 3 occurred at norm at most \(m+8\) inside the searched box. The box enumerator agrees on \((6,9)\) and \(b=(0,1,1,1)\): coefficient 0 at norm 15 and first nonzero coefficient \(-2\) at norm 23. The pairs \((6,9)\) and \((4,100)\) are outside the 21 stored controls.

The stored controls are the seven parameter pairs \((3,3),\ (3,4),\ (3,5),\ (4,4),\ (5,7),\ (8,11),\ (20,31)\), each with the three columns \(b\). Their box replay reports 21 cases and `uniform_proof_performed: false`. `python scripts/run_tests.py --all` passed all six tests in 0.024 seconds. Those tests reject \(s=2\), a zero claimed coefficient, a bound that stops before \(m+8\), a changed shell entry, and a missing case. Representative 131 is \(G(3,3)\), determinant \(15\), with first support 6 and first nonzero 14. Representative 132 is \(G(3,4)\), determinant \(22\), with first support 7 and first nonzero 15. The six averaged-support cancellations in `data/rank4-remainder/det24.json` are exactly these two Gram matrices and the three columns \(b\). The source hash on that pack is the accepted census hash.

## 2. Orthogonal sums

No correction to the corollary is needed. Factorization of the signed series on a block-diagonal Gram is the identity in `docs/foundations-and-exact-enumeration.md` §7. A holomorphic product that vanishes on the upper half-plane has a vanishing factor, the rank-at-most-three converse supplies an epsilon-one isometry of that factor, and extension by the identity keeps epsilon equal to 1. Transport by \(U\in\mathrm{GL}(4,\mathbb Z)\) is the same rule accepted with the Walsh review.

The plane-norm bound is correct. A Gauss-reduced basis \(\begin{pmatrix}a&b\\b&c\end{pmatrix}\) of a summand of determinant \(d\) satisfies \(1\le a\le c\) and \(|b|\le a/2\). Then \(d\ge 3a^2/4\). If \(a=1\), then \(b=0\) and \(c=d\). If \(a\ge 2\), then \(c\le d/a+a/4\le 4d/(3a)\le d\). Both basis vectors therefore lie in the sphere of squared norm \(d\le\Delta\). The saved search uses the uniform sphere bound 24. A negative answer is recorded only after that sphere and every saturated pair have been scanned: representative 131 records 374 primitive vectors up to sign and \(374\cdot 373/2=69751\) pairs. The tests compare the line, plane, and indecomposable outcomes of the two vector backends on \(I_4\), \(A_2\oplus A_2\), and the determinant-4 Gram of the \(D_4\) class, and they reject a projector that fails idempotence or \(G\)-adjointness.

The remainder pack’s summary is 138 line summands, 4 plane summands, 27 indecomposable classes, 73 even symmetry zeros, and 6 cancelled first-support pairs. That agrees with the census review’s decomposition count. This pass did not repeat the census replay.

## 3. The fibre identity

Write \(G=\begin{pmatrix}A&C\\C^t&D\end{pmatrix}\) with \(A\) positive definite, and split \(y=(y_1,y_2)\), \(a=(a_1,a_2)\), \(b=(b_1,b_2)\). Then

\[
Q(y)=y_1^tAy_1+2y_1^tCy_2+y_2^tDy_2.
\]

Completing the square in \(y_1\) produces

\[
Q(y)=(y_1+A^{-1}Cy_2)^tA(y_1+A^{-1}Cy_2)+y_2^tSy_2,
\]

where \(S=D-C^tA^{-1}C\). On the support, \(y_2=w\equiv a_2\pmod 2\) and \(y_1=a_1+2x\). The shifted vector is \(2x+\delta(w)\) with \(\delta(w)=a_1+A^{-1}Cw\). The sign \((-1)^{((y-a)/2)\cdot b}\) factors as the product of the outer sign in \(w\) and \((-1)^{x\cdot b_1}\). Summing first in \(x\) and then in \(w\) is the identity in the note.

Both sides are sums of absolute values bounded by a positive definite theta series, which converges on the upper half-plane. Absolute convergence justifies the rearrangement.

For \(G(3,3)\) and \(w=(1,1)\), the same algebra gives \(\delta=(-4/3,-2/3)\) and Schur matrix \(\begin{pmatrix}7/3&-2/3\\-2/3&7/3\end{pmatrix}\). The fibre contribution of \(x=0\) is \(8/3\), and \(w^tSw=10/3\), so the two exponents sum to 6. Direct evaluation of \(Q(0,0,1,1)\) is also 6. The shift \(-4/3\) is not an element of \(\frac12\mathbb Z\). The set \(\{2x+\delta\}\) is not an integral coset of a binary integral lattice. The half-characteristic theorems in ranks one through three therefore do not apply to this inner series. When \(u=-v\), one has \(\delta=0\), and that one fibre is an integral binary series; the neighbouring fibres are not.

The stored regression `reports/rank4-remainder/fibre-identity-check.json` records 13,689 rational norm identities and 219,024 sign factorizations on the 169 census Grams. Its own scope line says that this is a transcription check. The derivation above is the proof of the identity.

## 4. What the identity does and does not isolate

The coefficient of norm \(N\) is a finite sum, over pairs \((w,x)\) with combined squared length \(N\), of the outer sign of \(w\) times \((-1)^{x\cdot b_1}\). Finiteness is the positive definite ellipsoid. This bookkeeping is not a new nonvanishing theorem.

In the \(A_2\) family the two charts at level 0 land on the same norm \(m\). Their signed counts are \(+1\) and \(-1\) before the antipodal factor, for \(b=(0,1,1,1)\), and the total is 0. The same counts at level 1 total \(-2\) after doubling. The later-shell lemma is exactly that evaluation. Each fibre is nonzero at its own minimum; the rank-four coefficient vanishes because those minima collide and the counts cancel. Fibrewise appeal to a lower-rank theorem would not see the cancellation.

The collision is not the general indecomposable pattern. For

\[
G=\begin{pmatrix}2&-1&-1&-1\\-1&3&0&0\\-1&0&4&0\\-1&0&0&4\end{pmatrix},
\]

the upper block is an indecomposable binary form of determinant 5. This rank-four Gram has determinant 56 and no orthogonal summand through the sphere of squared norm 24. For \(a=(0,0,1,1)\) and \(b=(0,1,1,1)\), the complete shell at norm 8 contains the six vectors

\[
(0,0,1,1),\ (2,0,1,1),\ (0,0,-1,-1),\ (-2,0,-1,-1),\ (0,0,1,-1),\ (0,0,-1,1),
\]

with signs \(+1,+1,+1,+1,-1,-1\). The coefficient is \(+2\), and it is the first nonzero coefficient. Several fibres meet at that norm, and their signed count does not vanish.

A proposed lemma saying that the first combined fibre shell always cancels is therefore false. A proposed lemma saying that some combined shell is nonzero, for an arbitrary indecomposable Gram and an epsilon-trivial pair, is the full remaining converse written in fibre coordinates. The identity does not supply a reduction that sits strictly between those two statements and still covers a new infinite family. The one infinite family in which the counts are evaluated for every parameter is the \(A_2\) family already proved.

The bounded continuation is to connect the fibre count to epsilon. If every combined level had signed count 0, the series would vanish and the accepted Walsh criterion would require an epsilon-one automorphism. The missing step is the converse direction for a collision: a trivial epsilon character should force some level’s signed count to be nonzero. The determinant-15 series is the regression. Its epsilon character is trivial, its count at norm 6 is 0, and its count at norm 14 is \(-2\). An argument that turns a trivial character into a nonzero count has to return those two numbers on \(G(3,3)\) and \(b=(0,1,1,1)\), and it has to return a nonzero first count on the determinant-56 Gram above.

## 5. What this pass ran and read

Ran: `python scripts/run_tests.py --all`; an independent solution search for \(H\) and \(R\); direct signed shells for twelve series, including \((6,9)\) and \((4,100)\); `box_shells` on \(G(6,9)\); the determinant, summand test, and norm-8 shell of the determinant-56 Gram; the fibre arithmetic for \(w=(1,1)\) on \(G(3,3)\).

Read: `Grok-rank4-later-shell.md`; `docs/rank-four-later-shell-lemma.md`; `docs/rank-four-decomposition-and-fibres.md`; `code/rank4_cancellation.py`; `tests/test_rank4_remainder.py`; the summaries of `data/rank4-remainder/det24.json`, `reports/rank4-remainder/later-shell-box.json`, and `reports/rank4-remainder/fibre-identity-check.json`; the summand and average routines in `code/rank4_remainder.py` used by the tests. The census certificate was not reloaded.
