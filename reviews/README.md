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

## Release 1.3.0

[rank-four-stage-two-review.md](rank-four-stage-two-review.md) is a Claude
Opus 5.5 review of the stage-two reductions written by a Codex session
([report](../docs/rank-four-stage-two-reductions.md)). It checks sections 2
to 6 of that note line by line, and Lemma A and Theorem 5 of the
[minus-type report](../docs/rank-four-minus-type.md) (Theorem E of the
follow-up note). It adds an independent computation with its own
automorphism code (`code/rank4_stage2_review_check.py`, 31,500 configurations,
no exception), and runs and records the plus-type class list through 497. No
errors were found. Independent human review of these proofs is still
outstanding.

## Release 1.4.0

[rank-four-theorem-s-final-audit.md](rank-four-theorem-s-final-audit.md) is
the final audit of Theorem S by a Codex session (Astra), which accepts the
theorem and its minus-type corollary as a computer-assisted proof. Before it,
the same reviewer's quick review set three tasks, all completed and preserved
in `reports/rank4-theorem-s/package/`: exact bounds, input checks and
explicit incomplete runs (Claude Opus 5.5; `documents/`); an independent
standard-library replay of all 502,740 cases with certificates for every
exclusion (Grok Bot; `independent_replay_grokbot/`); and a check of the fan
against the primary literature and of the proof chain (Grok;
`literature_audit_grok/`). The final audit's own checker recomputes every
disposition without the classifier or its linear-programming backend; its
release copies run in CI. Independent human review is still outstanding.

The briefs themselves are handoff prompts and are not distributed.
