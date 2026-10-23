# Six signed minima in one common lattice

Astra review, 2026-10-22. The reduction below is proved algebraically and checked by an isolated program. It is **not a proof of C1 or C5**, and it does not assert that the sought indecomposable minus-type pairing exists.

## 1. Replace three index-two lattices by one index-four lattice

Let L be even, positive definite, of rank four and odd determinant. Write B for its bilinear form and N(x)=B(x,x). Let P be a two-dimensional B-isotropic subspace of L/2L. Since B modulo 2 is nondegenerate, P=P-perpendicular.

Define

\[
 K=\{y\in L:y\bmod 2L\in P\}.
\]

Then [L:K]=4 and 2L is contained in K. For any basis a,c of P, represented in L,

\[
 K=\{x:B(a,x)\equiv B(c,x)\equiv0\pmod2\}.
\]

For a nonzero a in P put L_a={x:B(a,x) is even}. Thus K is contained in L_a with index two.

The shift with character c, where c is any other nonzero element of P, sums over y=a+2x with sign (-1)^{B(c,x)}. Its two halves are indexed by epsilon=B(a,x) modulo 2. Their norms occupy disjoint residues, because

\[
 N(a+2x)\equiv N(a)+4B(a,x)\pmod8.
\]

Each half a+2x_0+2L_a lies in K and is the disjoint union of **two ordinary cosets of 2K**. The sign is constant on each: translating x by k in K does not change B(c,x) modulo 2. The signs on the two cosets are opposite, since B(c,.) is nonzero on L_a/K.

Consequently, with the sign orientation chosen appropriately, every signed half has the exact identity

\[
 S_{a,\epsilon}(q)=\theta_{v+2K}(q)-\theta_{w+2K}(q).
\]

This holds at every coefficient, not merely at the minimum. In particular, its minimal shell balances exactly when the two ordinary cosets have **equal minimum norm and equal minimal multiplicity**. This formulation includes shells with more than four vectors.

If B_K is a matrix of basis columns for K in L, then every entry of B_K^t G B_K is even. Set

\[
 H=\tfrac12 B_K^t G B_K.
\]

H is integral and det H=det G. On a {1,1,0} plane H is odd; on a totally singular plane H is even. Scaling all norms by two does not affect the comparisons or automorphisms.

## 2. The six pairs form one linear involution

Let V=K/2K and W=2L/2K. V has 16 elements, W has four, and V/W=P. There is a perfect pairing

\[
 W\times P\longrightarrow\mathbf F_2,
 \qquad(2x\bmod2K,a\bmod2L)\longmapsto B(x,a)\bmod2.
\]

For each nonzero a in P, let j(a) be the unique nonzero element of W annihilating a under this pairing. Set j(0)=0. In dimension two over F2, j is linear and invertible: for a basis a,c, j(a+c)=j(a)+j(c).

If pi:V->P is the quotient, define

\[
 \tau(v)=v+j(\pi(v)).
\]

Then tau is linear, tau squared is the identity, and both the image and kernel of tau-I are W. It fixes the four W-classes and swaps the remaining twelve classes in exactly the six pairs above.

Thus all six minima balance if and only if the ordinary parity-coset minimum function **and its multiplicity function** are invariant under tau. Full vanishing is equivalent to invariance of the entire sixteen-component coset-theta vector under tau; the four fixed components impose no conditions.

An integral automorphism T of K whose reduction modulo 2 is tau is an exact certificate of all six identities. It preserves W, hence extends to L. There is no need for a Sturm expansion when such a witness is supplied. **The existence of that lift from minimum data alone remains to be proved.** A theorem saying that minimum data determine an *unmarked* lattice does not by itself prove this prescribed lift exists.

## 3. A smaller, finite classification target

Parity-coset minima are linear functions of the Gram matrix on each iso-edge domain, including its closure. The relevant result is Lemma 5.2 of [Dutour Sikiric and Kummer, *Iso Edge Domains*](https://arxiv.org/pdf/2102.11139). This makes the six minimum-norm equalities linear on a known cell.

In dimension four the published classification has 52 Delaunay subdivision types; the corresponding coarser decomposition has 51 iso-edge types. There are three full-dimensional secondary cones. See [Vallentin, Weissbach and Zimmermann, *The chromatic number of 4-dimensional lattices*](https://arxiv.org/pdf/2407.03513). Their unmarked classification is a starting point, not a ready-made classification of this pairing problem.

There are only **210** possible linear involutions tau on a fixed F2^4 with rank(tau-I)=2: choose W=ker(tau-I)=im(tau-I) in 35 ways, then an isomorphism V/W->W in 6 ways. Enumerating these on each classified cell covers the possible markings; stabilizers can reduce duplication. The crude product 52*210=10,920 is a count of cell/marking combinations before refinement, not a runtime estimate or a claim that one representative point settles each combination.

The next tool should:

1. Import and validate the exact ray/inequality descriptions of the three maximal cones and their faces, with source hashes and a coverage record. Keep faces containing positive definite forms, including cases with extra minimal vectors.
2. Attach complete sets of minimizers for all 15 nonzero parity classes. Verify their parity labels and constant multiplicities on each relative interior.
3. Enumerate the 210 markings, impose the six equal-minimum equations, and reject incompatible multiplicities. Use exact rational feasibility certificates. If the intersection lies on a face, process that face with its enlarged minimal sets.
4. Enforce the arithmetic that reconstructs the original even L. In K coordinates, with Gram 2H, let L_W consist of the half-vectors whose doubled classes lie in W. It is even exactly when H[w] is 0 modulo 4 for every w in W. Require det H odd, the correct induced pairing tau, and the target determinant residue. Real feasibility alone does not satisfy these conditions.
5. For each surviving feasible family, either certify the required symmetry over the whole family or return an integral Gram matrix without that certificate. A symmetry found at one sampled interior point does not certify the entire family.
6. Send any surviving indecomposable minus-type matrix through the existing exact three-shift/Sturm/isometry checks. Six balanced minima alone are not a hit.

This replaces an unbounded determinant search with a proposed finite geometric and arithmetic classification. Completing and certifying that classification is still work to do. It can potentially find a counterexample or prove the restricted search universe empty.

## 4. Tests and reusable certificate format

[common_lattice.py](results/rank4-night-review-20261022/common_lattice.py) constructs K by integer HNF, checks its index and determinant, constructs all sixteen classes and tau, and compares two independently enumerated descriptions: the builder's mod-4 histogram of L and the verifier's rectangular-box parity histogram of K.

The check covered **19 matrices, 271 marked planes, and 1,626 signed sub-coset series through norm 200**. Every coefficient agreed. It includes the decomposable determinant-21 control, the determinant-9373 five-of-six example, twelve sampled production matrices, and all five indecomposable plus-type controls from the handoff. Both ordinary minima were reached in every tested pair.

For all eleven six-of-six planes in this test, the program found an integral K-automorphism inducing tau and recorded its extension to L. This upgrades those particular zeros to exact all-order symmetry certificates. The five plus controls are included. Results: [common-lattice-checks.json](results/rank4-night-review-20261022/common-lattice-checks.json).

A separate certificate pass applied the same criterion to every saved box-sweep hit: **10,522 matrices and 63,132 marked planes**. All have an integral splitting witness and all marked planes have an automorphism inducing tau. These witnesses are saved in [six-hit-certificates-all.json.gz](results/rank4-night-review-20261022/six-hit-certificates-all.json.gz); generating them took 125.88 seconds on one worker. This settles these particular finite examples exactly, without proving the general minimum-to-symmetry implication.

For a compact certificate, store G, P, the basis B_K, and T. An independent checker only needs to verify integral invertibility, T^t(B_K^t G B_K)T=B_K^t G B_K, and T modulo 2=tau. For a splitting certificate on L, also check a nontrivial h with h squared equal to I and h congruent to I modulo 2; then (I+h)/2 and (I-h)/2 are integral orthogonal projections. This closes the gluing gap that a rational orthogonal span leaves open.
