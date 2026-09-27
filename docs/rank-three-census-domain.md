# Rank-three determinant <=24 symmetry-converse census

Complete determinant-24 domain-generation proof and exact coverage policy. The completed census has 120 isometry classes and all 4,320 even pairs resolved. See the lower-rank manuscript and reproduction guide for final results and replays.

## Domain and certificates

Enumerate all positive definite symmetric integral rank-three Gram matrices
up to GL3(Z) isometry, with determinant 1 through 24. Include odd, even,
nonprimitive and decomposable lattices. Different positive scalings remain
distinct. Squared norm is x^t G x. Each representative has all 36 even and
28 odd binary primal/dual characteristic pairs; odd pairs vanish by -I.

An even zero needs an exact symmetry witness or an independently proved
identity. A nonzero needs a complete nonzero signed coefficient. Complete
cancellation through a finite norm4 bound remains unresolved. A symmetry-converse counterexample would additionally require epsilon zero on the entire characteristic stabilizer. No such example occurs in this domain.
Full groups, rather than selected ambient generators, are required.

## Elementary finite-generation proof

Choose a shortest nonzero vector v1, of squared norm A. It is primitive:
an integer multiple of another vector would not be shortest. A primitive
integer column extends to a unimodular basis (Bezout / integer elementary
operations). Orthogonal projection along v1 identifies L/Zv1 with a rank-two
Euclidean lattice. The projection is discrete: the projections of either
completion basis are linearly independent and generate it.

Gauss-reduce a basis w2,w3 of this projected lattice. Subtract a nearest
integer multiple to reduce the pairing, and interchange when the second
vector is shorter. Every interchange strictly decreases the first squared
norm; projected norms lie in (1/A)Z, so reduction terminates. Thus its Gram
entries S,Q,R satisfy S<=R and |2Q|<=S. Lift this basis to v2,v3 in L and
subtract integer multiples of v1 independently to reduce both pairings.
These operations preserve a basis. Write

    G = [[A,B,D],[B,C,F],[D,F,E]],
    |2B|<=A, |2D|<=A,
    S=C-B^2/A, Q=F-BD/A, R=E-D^2/A,
    |2Q|<=S<=R.

Since v2 is nonzero, C>=A. Hence S>=3A/4. For determinant Delta,

    Delta=A(SR-Q^2)>=3 A S^2/4>=27 A^3/64.

Consequently A is a positive integer with 27 A^3<=64*24, hence A<=3.
For an arbitrary declared bound H<=24, impose 27 A^3<=64H and
3(AC-B^2)^2<=4AH. Also C<=sqrt(4H/(3A))+A/4.
This proves finite A,B,C,D bounds without a ternary classification source.
Let M=AC-B^2>0. The projected pairing inequality is

    |2(AF-BD)|<=M.

It gives a finite interval for integer F. For each Delta=1,...,H set

    E=(Delta+C D^2+A F^2-2BDF)/M.

Retain exactly integral E with AE-D^2>=M. The leading minors A,M,Delta
are positive, proving positive definiteness. This enumeration includes a
basis for every lattice in the domain. Additional generated bases are
allowed; no assertion of canonical reduction is needed.

For simple independent replay loops, both C and E are <=2H+A: indeed
S<=sqrt(4H/(3A))<=2H, C=S+B^2/A<=2H+A/4, and
E=Delta/(AS)+Q^2/S+D^2/A<=4H/(3A^2)+S/4+A/4
<=11H/6+A/4<=2H+A. These deliberately loose bounds require no radicals.
The builder and replay use the same proved inequalities but different
enumeration loops (determinant-solving versus bounded diagonal scanning).

## Exact coverage and interruption policy

Enumerate all target vectors with the prescribed source basis norms, then
every compatible ordered triple. Check each required pairing, determinant
+/-1 and U^t H U=G. Every possible isometry sends the three source vectors
to such a triple, proving exhaustive isometry testing. Taking H=G proves
the full automorphism group. A comparison interrupted by a node limit
retains the candidate with uniqueness unresolved; an interrupted group
search supplies no group or negative epsilon claim. Partial shells are
discarded on interruption. Group and shell limits are separate per call.

Use LDL sphere enumeration for construction; replay uses inverse-Gram
coordinate boxes and exhaustive basis-image products. The box inequality
is |y_i|^2<=bound*(G^-1)_ii. The independent checker regenerates the finite
candidate domain, verifies every generation basis change, proves retained
representatives pairwise nonisometric, recomputes full groups and all
stabilizers, and directly sums complete signed shells. Optional LDL shell
replay supplies a second report. The finite domain, every characteristic
verdict and all unresolved cases must be saved.

Initial resource choices: norm4<=32 and 1,000,000 nodes per enumeration
or basis-image search. These are not coefficient identity cutoffs.
The modular identity cutoff and its checked hypotheses are given in the separate coefficient-cutoff appendix.

Controls include Z^3, orthogonal binary/rank-one sums, scaled forms, A3,
noncanonical unimodular bases, and the determinant-29 first-shell trap
(a regression only, outside the census domain).
