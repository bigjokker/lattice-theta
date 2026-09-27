# Checked citations for the rank-three review snapshot

Prepared 2026-09-27 for `papers/low-rank-and-root-lattice-results.md`. Earlier checks and access limits
remain in [docs/literature-and-source-scope.md](literature-and-source-scope.md). This addition records the exact
primary-source passages used for rank-three attribution and the cutoff.
It does not turn a limited search into a novelty claim.

## Classical geometry

[Conway and Sloane, Low-dimensional lattices VI, author text](https://neilsloane.com/doc/fedorov.pdf)
was reopened. The formula used for geometric attribution is equation (5)
in section 2.0.5, the proof following Theorem 3; Theorem 3 itself concerns
Voronoi vectors. Theorem 8 states three-dimensional existence of an obtuse
superbase. Section 8.0.2 defines the labeled Delone graph. These are the
specific citations used in the new manuscript.

[Kurlin, author PDF](https://kurlin.org/projects/lattice-geometry/lattices3Dmaths.pdf)
was reopened at Definition 2.5, Theorem 2.8 and Appendix A, Lemma A.1.
The definition permits zero conorms. The lemma gives the adjacent-superbase
move used here. The manuscript independently proves termination for integral
lattices using the sum of four squared lengths, whose decrease is
2 times the positive pairing. This differs from the partial-sum vonorm
quantity in the source. These citations provide geometric background;
the signed-shell proof is included in the manuscript itself.

## Coset modularity and the complex cutoff

[Kane-Kim, arXiv:2211.03987v2](https://arxiv.org/html/2211.03987v2)
was reopened at the level/discriminant definitions, section 2.2 and
Proposition 2.3. The level clears the inverse Gram denominators. The coset
series uses q raised to the squared norm. In odd rank the space in that
proposition has weight k/2, level 4 N_L a^2 with its stated conductor
intersection, and character chi_(4d_L). The project's common-space and
trivial-coset translations remain in docs/coset-modularity.md; they are deductions,
not a spinor-genus consequence.

[Brunault, Sturm bounds for general congruence subgroups](https://perso.ens-lyon.fr/francois.brunault/recherche/Sturm-bound-general.pdf)
was reopened at Theorem 1 on PDF page 1 and the cusp-width definition.
Its complex identity bound uses the index of the group enlarged by -I
and order measured in the cusp parameter. For Gamma0(M), -I is already
present and infinity has width one. Applying the theorem to the holomorphic
weight-six fourth power gives the inequality used in Theorem 7.
The relation between its order and the signed norm4 index, divisibility
of the index by eight, and the inclusive cutoff are proved in the manuscript.
No finite-field congruence theorem or optimal-level claim is inferred.

## External search and its limits

[Grok's prior-art report](../reviews/low-rank-literature-review.md) is preserved as an
external research report. It records the portions read and unavailable leads.
It found no predecessor within that search and reported no contradiction.
The new manuscript attributes classical background and presents its own proof;
it makes no claim that the signed theorem is new.

Two wording points from the report are resolved during integration:

- An even number of odd factors gives even total parity. In the decomposable
  rank-three even zeros at issue here, exactly two factors are odd.
- The Selling norm expansion in the Conway-Sloane author text is located
  precisely in equation (5), section 2.0.5, rather than the statement of
  Theorem 3.

Selling's French original, Delone, Igusa, Schiemann's primary article,
the additional Mumford pages and other leads in Grok's report are not needed
as new direct citations in this manuscript. Their descriptions retain the
report's external or unread status; this closing check does not certify
them. No further open-ended search is a requirement for this snapshot.
