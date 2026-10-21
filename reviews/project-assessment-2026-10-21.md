# Project assessment, 21 October 2026

Curated copy (release 1.1.0) of the independent assessment PROJECT-REVIEW-2026-10-21.md, written on 2026-10-21 by a separate Claude session before this release update. Its statements about the state of the workspace and of the release are historical; the release update it recommends is version 1.1.0. File paths are relocated to the release layout.

The original project has reached a defensible scientific stopping point. The
October work has also completed the four research tasks in
the workspace to-do list (not distributed). There is worthwhile mathematics left,
but it is a new phase, with a narrower question to choose. My recommendation
is to consolidate the results first, then pursue the even, odd-determinant
parity question if you want another research phase.

This assessment distinguishes the original symmetry-converse problem from a
complete classification of the parity property in rank four. The former is
settled; the latter is not.

**What the new work accomplished.**

| Question | Current result | Scope |
|---|---|---|
| Does every vanishing even characteristic have an automorphism witness? | Yes in ranks one through three; no in rank four. | The original scientific objective is complete. |
| Can rank-four vanishing persist even when the two associated coset lattices are not isometric? | Yes, at determinant 1729, using a re-embedding of the Conway–Sloane pair. | A stronger obstruction than the original determinant-23,716 construction. |
| Can the original symmetry counterexample be made smaller? | Yes: an even lattice of determinant 1156, with full group {I, −I}. | Smallest in the searched Eisenstein family and range; global minimality is not proved. |
| Does indecomposability plus determinant not divisible by four imply the parity property? | No. The smallest odd-lattice failure under these hypotheses has determinant 49. | Minimality follows from the complete restricted lists and the earlier census. |
| Does the conjecture become true for even lattices of odd determinant? | No. The smallest indecomposable failure has determinant 121. | The only failure among 255 indecomposable classes through determinant 201. |

The determinant-1729 example is the most substantial conceptual advance. The
determinant-1156 example improves the size and presentation of an already
known mechanism. The 1729 example shows that a zero need not even come from
an isometry between the two coset lattices. It still has a piecewise
norm-preserving bijection, so the precise claim is failure of an explanation
by one Euclidean isometry, rather than absence of every geometric
description.

The [literature note](../docs/rank-four-isospectrality-literature.md) now gives a useful framework: a
vanishing characteristic is a theta identity between two lattices in a
Klein-four configuration. The strongest result should be presented through
that bridge, crediting Conway–Sloane for their isospectral lattices and
Cervino–Hein for their general non-isometry theorem. The project supplies
the compatible index-two re-embedding and the characteristic-zero
construction. This is a precise contribution statement; it is not a
priority claim. See [Conway–Sloane, remark (x) and Theorem 1](https://neilsloane.com/doc/nipp.pdf)
and the project's [level-four review](rank-four-level-four-review.md).

The parity work is also more than two negative answers. It proves that an
even indecomposable rank-four lattice of odd determinant has no involution
witness. It proves that an automorphism whose square is −I forces six
zeros on its associated totally singular plane. These explain why the
determinant-121 counterexample is different from the odd determinant-49
example. They are the strongest starting point for another theorem.

**Two useful deductions from the existing examples.**

The determinant-961 pair deserves a short explicit corollary:

> The ordinary theta series does not determine the parity property in rank four.

The first member of the Conway–Sloane remark-(x) pair is decomposable and has
36 vanishing even characteristics. The second is indecomposable and has no
vanishing even characteristics. These decisions are in the parity
certificates and were reproduced by the independent parity check during
this review. The ordinary theta series of the two lattices agree, so their
ordinary representation numbers cannot recover this characteristic-vanishing
behavior. This statement concerns all rank-four lattices; the example does
not settle the same question restricted to pairs of indecomposable lattices.

To avoid relying on a long matching prefix or on the suggestive wording of
remark (x), I checked an explicit modular proof for this particular pair.
For the displayed odd Gram matrices P and Q, use even Gram matrices 2P and
2Q. Both have level 124 and square determinant, so their theta series are
weight-two modular forms with trivial character on Gamma_0(124). Its index
is 192, and the inclusive coefficient bound is 32. Exact LDL enumeration
and independent rectangular-box enumeration agree for every coefficient
from 0 through 32. The theta series are therefore identical. The two
independent isometry searches also find no isometry.

The same small computation resolves the status of the determinant-3844 pair
S1, S2 in section 2 of [rank-four-level-four-zero.md](../docs/rank-four-level-four-zero.md).
That report claims only agreement through norm 2000. Both Grams are already
even, have square determinant and level 62. The group index is 96, so the
inclusive weight-two bound is 16. Both enumeration backends agree through
that bound, and both isometry searches return no isometry. Thus this is an
all-level isospectral non-isometric pair. It is not merely a prefix candidate.
This does not establish that the pair is new to the literature.

The exact general theorem used here is: for a positive definite even Gram
matrix A of rank four and level N (the least N for which N A^{-1} is even
integral), the series sum_x q^{x^t A x/2} belongs to
M_2(Gamma_0(N), chi_det(A)). In these cases the determinant is a square,
so the character is trivial. Two such forms are equal if their coefficients
agree through floor([SL_2(Z):Gamma_0(N)]/6), inclusively. Here
q = exp(2 pi i tau). For the coefficient criterion, see section 8.2 of
[Best et al., Computing classical modular forms](https://people.cs.rutgers.edu/~sg1108/math/Mathematics/cmfs-mf-021120.pdf).
The script uses norm4 = 8 times the Fourier exponent, so the normalization
is explicit.

Reproduction and evidence:

```powershell
python code/check_small_pairs.py
```

The [script](../code/check_small_pairs.py) and
[saved output](../reports/project-review-2026-10-21/small-pairs.json) contain both
Gram pairs, levels, coefficients, and the independent checks. These are
additional deductions made during this assessment, rather than results
already claimed by the October 20 reviews.

**The best next research question.**

I would prioritize this precise conjecture:

> Every even, orthogonally indecomposable, positive definite rank-four
> integral lattice with determinant congruent to 5 modulo 8 has the parity property.

The evidence is all 131 distinct minus-type indecomposable classes through
determinant 201, with every relevant characteristic decided. The report's
larger count of 134 repeats three classes from the old census; the review
already identifies that overlap. A finite census supports the conjecture
but does not prove it.

There are two logically separate tasks. First, close the symmetry argument:
involutions are excluded, and an automorphism squaring to −I is excluded by
minus type, but the report leaves other order-four actions and order-eight
actions untreated. A useful initial target is a complete analysis of those
remaining integral representations and their signs on characteristic
stabilizers. A successful theorem here would classify symmetry witnesses,
even if it did not settle every theta identity.

Second, rule out or construct zeros with no symmetry witness. The
determinant-1156 and determinant-1729 examples demonstrate why excluding
automorphisms cannot finish this step by itself. One possible stronger
target is to prove that every zero in the even, indecomposable,
odd-determinant setting forces span(a, G^{-1}b mod 2) to be a totally
singular plane. If true, minus type follows immediately because it has no
such plane. At present this is a proposed route, not a necessary condition
established for arbitrary identity zeros.

I would start with the remaining order-four/order-eight cases and only
extend the census to test an explicit obstruction arising there. The
completion criterion for that first stage is a witness-classification
theorem or a certified counterexample to it. This is a more informative
objective than another larger table with no failures.

The companion plus-type question is whether the square-root-of-−I
mechanism exhausts all zeros under the same evenness, indecomposability
and odd-determinant hypotheses. It is attractive but broader: the data
currently contain only one failing indecomposable class in that universe.
I would keep it as the second stage.

**A separate, more computational option.**

The hierarchy has certified examples at levels 1, 2 and 4, while level 3
remains unexhibited in the examined domains. The precise question is whether
the associated lattices can be isometric while no isometry between them
preserves their common sign kernel. This is a clean existence question.

A targeted search could start with a lattice M, compare its fifteen
index-two sublattices, and look for isometric kernels in different Aut(M)
orbits. Transporting M along a kernel isometry provides candidate
configurations. Each candidate still needs an integral third lattice and
the full kernel-group check that excludes levels 1 and 2. A bounded failure
would not prove level 3 impossible. I rank this below the minus-type problem
because it fills an explanatory distinction rather than classifying a
substantial family.

Global determinant minimization is also open. The existing reports' safe
general bracket is 24 < det <= 1156. Their family-minimum statement should
not be promoted to a global minimum. The newer restricted censuses give
more exclusions, but do not cover all odd lattices with determinant divisible
by four above 24. I would not make a broad census to 1155 the next default
task. Likewise, moving to rank five is unnecessary to establish the
dimension-four threshold already proved.

**What remains before the work is presented as a current release.**

The local research is substantially ahead of its entry points. The workspace
README still describes the earlier low-rank checkpoint. The local Github
checkout is clean, with its latest commit dated September 27; its README
and manuscript still lead with determinant 23,716. The October builders,
certificates, reviews and structural results have not been incorporated
there. I did not check the remote branch or current hosted CI status.

A focused release update should bring the introduction, theorem statements,
literature discussion and reproduction guide into agreement. I would lead
the rank-four discussion with the 1729 obstruction, retain 1156 as the
smaller symmetry counterexample, and give the odd-determinant parity results
their own treatment. The two old parity conjectures should be marked false
where they are still presented as proposals. The new scripts and certificates
should be part of the curated reproduction path.

Carry the existing review corrections into that presentation: distinguish
distinct classes from repeated observations; distinguish the first stored
witness from all possible witness orders (D4 also has involution witnesses);
and distinguish the 26,775 decisions in Task B from its 26,769 nonzero
decisions, and the 120,120 decisions in Task E from its 120,065 nonzeros.
These are reporting corrections, not failures of the headline results.

For mathematical publication, independent human review of the core proofs
and the attribution is the useful next validation step. The existing reviews
are AI-assisted local reviews with computational replication, not journal
peer review. There is no need to settle every remaining rank-four question
before preparing a coherent manuscript. Authorship and public submission
are separate choices.

**Checks performed for this assessment.**

I read the three October reports and their reviews, the literature bridge,
the earlier closure and continuation notes, the relevant tests and
enumeration code, and the local release history. I ran:

```powershell
python scripts/run_tests.py   # the three October test modules are part of this default suite
python code/rank4_parity_review_check.py
python code/check_small_pairs.py
```

All 36 tests passed in about 70 seconds. The independent parity script
reproduced the determinant-121 and determinant-49 results, the pair-member
parity decisions and the E8 checks. The added modular checks passed for
both small isospectral pairs. I did not regenerate the entire large class
lists; the minimality claims above rely on their existing complete
enumerations and independent review records. No existing research files,
release files or published artifacts were changed.

My recommendation is to finish the current presentation and release, and
then choose the minus-type conjecture if further research is desired. The
project has enough completed mathematics to stop responsibly now; continuing
would pursue a sharper theorem, rather than repair an unfinished original goal.
