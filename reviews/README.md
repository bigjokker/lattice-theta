# Local review evidence

These documents record AI-assisted mathematical audits. They are included for
transparency about proof checking and revisions, and are not presented as
journal peer review or independent published research. The manuscripts, exact
certificates, and reproducible verifier outputs carry the claims and can be
checked without trusting the reviewers.

## Release 1.0.0 reviews (Grok and Astra)

The [focused rank-four review](rank-four-counterexample-review.md) accepted
the all-level coset bijection, normalization, complete 315-point enumeration,
Aut(G)={I,-I}, and epsilon zero on both elements. Its scope is this
counterexample; it does not re-review every earlier result or establish
minimum determinant or literature priority.

The remaining 1.0.0 reviews cover the lower-rank proofs, source scope,
modularity, finite domain, and later-shell regressions. Their chronological
observations are retained as audit history. Handoff prompts and assignments
are excluded.

## Release 1.1.0 reviews (Claude Fable 5.1 sessions)

Each October 2026 task was specified in a brief by one Claude Fable 5.1
session, executed by a second session, and then reviewed by the brief-writing
session with its own independent replication code (included under `code/` as
`rank4_*_review_check.py` and `rank4_*_review_enum.py`, outputs under
`reports/`). The review notes record what was checked, with which tools, the
verdict, and precision points:

- [rank-four-level-four-review.md](rank-four-level-four-review.md): the
  determinant-1729 level-four zero and the eight configurations; accepted.
- [rank-four-eisenstein-review.md](rank-four-eisenstein-review.md): the
  determinant-1156 counterexample and the family enumeration; accepted, with
  an independent re-enumeration of the family.
- [rank-four-parity-review.md](rank-four-parity-review.md): the refutation of
  both parity conjectures; accepted, with an independent re-enumeration of all
  three class lists from a different reduction box and a line-by-line check of
  the four theorems.

[project-assessment-2026-10-21.md](project-assessment-2026-10-21.md) is an
independent assessment of the whole project by a separate Claude session. It
reran the test suites and the independent parity check, added the two
modular-form deductions of `code/check_small_pairs.py`, and recommended the
release update that became version 1.1.0 and the research directions now
listed in [deferred research](../docs/deferred-research.md). Its remarks on
the state of the workspace and of the release describe the state before this
release.

## Release 1.2.0

The determinant census through 279 ([report](../docs/rank-four-minimum-determinant.md))
has no separate review note. Its evidence is computational: every one of its
241 per-determinant chunks is a decided class list in the schema reviewed for
release 1.1.0, and each is replayed by the unchanged independent verifier
`code/verify_rank4_parity.py` (`reports/rank4-gap-census/*-replay.json`). The
proof of the minimum also uses the accepted determinant-24 census and the
orthogonal-sum corollary of `docs/rank-four-decomposition-and-fibres.md`.

The briefs themselves are handoff prompts and are not distributed.
