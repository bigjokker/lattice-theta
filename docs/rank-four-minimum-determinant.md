# The rank-four symmetry converse: the minimum counterexample determinant is 280

Curated copy (release 1.2.0, 2026-10-21) of the research report
RANK4-GAP-CENSUS.md of 2026-10-21. File paths are relocated to the release
layout. Executed by a Claude Fable 5.1 session after the determinant-280
correction of release 1.1.1 had lowered the bracket for the minimum
counterexample determinant to 24 < det <= 280. Python 3.11 standard library,
exact arithmetic throughout, 16 worker processes. Every count below names its
universe and its deduplication criterion.

Notation as in [rank-four-parity-property.md](rank-four-parity-property.md).
The symmetry converse for a lattice L: every even characteristic (a, b) with
Theta_L[a,b] identically zero has a witness g in Aut(L) stabilizing (a, b)
with epsilon(g) = 1 (level 1). A counterexample is a zero of level 2, 3 or 4.

## 1. Verdict

**Theorem.** Every positive definite integral rank-four lattice with
det <= 279 satisfies the symmetry converse. Hence the minimum determinant of
a rank-four counterexample is **280**, attained by the odd lattice
G = [[3,1,1,0],[1,5,-1,2],[1,-1,4,-1],[0,2,-1,7]], a = (1,0,0,1),
b = (1,1,1,1) of release 1.1.1 (`data/rank4-level4/det280-level4.json`;
level 4, Aut = {I,-I}). The bracket 24 < det <= 280 is sharp.

The proof is a finite census (sections 2-4) plus two uniform facts: the
accepted determinant-24 census, and the orthogonal-sum corollary for
decomposable lattices. Section 5 lists what the det-280 chunks themselves
show.

## 2. Coverage of 25 <= det <= 279

A lattice is even or odd; an even rank-four lattice has det = 0 or 1 mod 4
(checked here: the complete generator over even lattices with 25 <= det <= 279
and no determinant filter produces 41,974 candidate bases, 22,982 with
det = 0 mod 4 and 18,992 with det = 1 mod 4, none with det = 2 or 3 mod 4).
So the universe splits into the rows below. "Complete list" means every
isometry class in the row's universe, deduplicated by exact isometry
(complete basis-image search, `rank4_feasibility.isometries`, after the
theta-invariant prefilter), with every indecomposable class decided on all
105 even characteristics.

| Universe | Source | Candidates | Classes | Indecomposable | Zeros | Levels |
|---|---|---:|---:|---:|---:|---|
| all lattices, det <= 24 | accepted census `data/census-rank4/det24-bound32.json` (169 classes, 1,510 candidate bases; `code/verify_rank4_census.py` rejects any zero without an epsilon-one witness, and the summary has no `proved_zero_without_symmetry` verdict; accepted in [rank-four-census-review.md](../reviews/rank-four-census-review.md)) | 1,510 | 169 | 27 | 6,376 even zeros | all 1 |
| odd lattices, 4 !\| det, 25..168 | Task E, `reports/rank4-parity/task-e-det168.json` | 46,494 | 3,487 | 1,144 | 55 | all 1 |
| odd lattices, 4 !\| det, 169..279 | **new**, region `odd-notdiv4-169-279` | 106,794 | 6,658 | 2,900 | 157 | all 1 |
| odd lattices, 4 \| det, 25..279 | **new**, region `odd-div4-25-279` | 98,629 | 7,043 | 2,577 | 4,042 | all 1 |
| even lattices, odd det, <= 201 | Task B, `reports/rank4-parity/task-b-det201.json` | 10,188 | 300 | 255 | 6 | all 1 |
| even lattices, odd det, 202..279 | **new**, region `even-odd-202-279` | 8,972 | 269 | 238 | 0 | - |
| even lattices, 4 \| det, <= 64 | Task C, `reports/rank4-parity/task-c-det64.json` | 1,101 | 69 | 27 | 122 | all 1 |
| even lattices, 4 \| det, 65..279 | **new**, region `even-div4-65-279` | 22,080 | 1,017 | 524 | 1,131 | all 1 |

An earlier workspace census of even lattices with det = 5 mod 8, 202..501
(703 indecomposable classes, no zeros; not distributed) overlaps the new even
odd-determinant region on 202..279, and the two agree. Decomposable classes
are not decided in any list: they satisfy the converse by the uniform
orthogonal-sum corollary
([rank-four-decomposition-and-fibres.md](rank-four-decomposition-and-fibres.md),
which rests on the rank one to three converse). Imprimitive lattices appear
among the candidates and are decided like any other class (they would in any
case reduce to det <= 17, covered by the census).

## 3. The four new regions

Universe of each region: all isometry classes of positive definite integral
rank-four lattices with the stated lattice parity, determinant range and
determinant filter. Generator: `rank4_parity_search.class_list`, the complete
reduced-basis builder of the parity report (the coverage argument is embedded
in every class-list document as `coverage_argument`; see also
[rank-four-census-domain-and-cutoff.md](rank-four-census-domain-and-cutoff.md)),
run once per determinant value with D0 = dmin = det, so each chunk is a
complete class list for one determinant in the standard schema
`work7-rank4-parity-v1-class-list`. Dedup within a chunk is exact isometry as
above; distinct determinants are distinct classes. Decision:
`rank4_parity_search.decide_class_list` on every indecomposable class (full
automorphism group by complete isometry search; for each of the 105 even
characteristics a witness in the stabilizer, else a nonzero complete shell
coefficient through the ladder 100/400/1600, else the modular cutoff and the
kernel level).

| Region | Determinants | Candidates | Classes | Decomposable | Indecomposable | P holds | P fails | Even zeros | Largest first nonzero norm4 | Largest Aut | Build wall / CPU |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| odd, 4 !\| det, 169..279 | 84 | 106,794 | 6,658 | 3,758 | 2,900 | 2,862 | 38 | 157 | 238 | 96 | 105 s / 1,370 s |
| odd, 4 \| det, 25..279 | 63 | 98,629 | 7,043 | 4,466 | 2,577 | 1,665 | 912 | 4,042 | 136 | 96 | 97 s / 1,185 s |
| even, odd det, 202..279 | 39 | 8,972 | 269 | 31 | 238 | 238 | 0 | 0 | 180 | 144 | 12 s / 117 s |
| even, 4 \| det, 65..279 | 53 | 22,080 | 1,017 | 493 | 524 | 289 | 235 | 1,131 | 140 | 384 | 24 s / 276 s |
| **total** | 239 | 236,475 | 14,987 | 8,748 | 6,239 | 5,054 | 1,185 | 5,330 | 238 | 384 | 238 s / 2,948 s |

Decisions: 6,239 x 105 = 655,095 characteristics. 5,330 are zeros, **every
one with a symmetry witness (level 1)**; 649,765 are nonzero by a complete
shell coefficient (the first nonzero shell has norm4 <= 238, inside the
ladder, so the modular cutoff was never needed); 0 unresolved. "P fails"
counts indecomposable classes without the parity property; those failures
are all explained by witnesses and are not counterexamples to the converse.

## 4. Independent replay

`verify_rank4_parity.replay_decided` (FAST_BOX) on every chunk: regenerates
the candidate domain with the loose-box generator (must equal the stored
domain), recomputes every full automorphism group and decomposability by box
search, re-decides every stored decision by box enumeration, and checks
pairwise non-isometry of all representatives of equal determinant. All 239
chunks: `verified`.

| Region | Pairwise non-isometry searches | Replay wall |
|---|---:|---:|
| odd, 4 !\| det, 169..279 | 284,230 | 189 s |
| odd, 4 \| det, 25..279 | 511,590 | 138 s |
| even, odd det, 202..279 | 1,864 | 8 s |
| even, 4 \| det, 65..279 | 11,626 | 13 s |

`reports/rank4-gap-census/closure-25-279.json` ties the chunk SHA-256s of
the build summaries to those of the replays, records the prior lists with
their hashes, and states `symmetry_converse_holds_for_25_to_279: true`.

## 5. Determinant 280 itself

Run after the closure, same method, regions `odd-div4-280` and
`even-div4-280` (one chunk each, both replayed `verified`). Universe: all
isometry classes of positive definite integral rank-four lattices with
det = 280, by lattice parity; dedup by exact isometry as above.

| Universe | Candidates | Classes | Indecomposable | P fails | Even zeros | Levels |
|---|---:|---:|---:|---:|---:|---|
| odd lattices, det 280 | 2,911 | 243 | 93 | 39 | 163 | 162 at level 1, **1 at level 4** |
| even lattices, det 280 | 486 | 27 | 9 | 6 | 24 | all 1 |

The single zero without a witness is in class 202 of the odd chunk,
G' = [[3,-1,-1,0],[-1,4,-1,-1],[-1,-1,5,2],[0,-1,2,7]], a = (1,0,0,1),
b = (1,1,1,1), Aut = {I,-I}, proved zero by the modular cutoff with kernel
level 4. G' is the certified det-280 lattice in another basis: with
T = [[-1,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,1]], T^t G' T = G (exact
isometry search, two isometries, both found). So **at determinant 280 there
is exactly one isometry class that fails the symmetry converse, with exactly
one failing characteristic**, the one of release 1.1.1. Outputs:
`reports/rank4-gap-census/{odd,even}-div4-280/det280.json.gz`,
`*-summary.json`, `*-replay.json`.

## 6. Files and regeneration

Script: `code/rank4_gap_census.py` (tasks `run`, `replay`, `summary`;
`--regions`, `--workers`, `--chunks`, `--out`). Preserved outputs in
`reports/rank4-gap-census/`: one chunk `<region>/detNNN.json.gz` per
determinant (241 chunks, 42 MB gzipped, read compressed; not listed in
`reports/compressed-reports.json`), `<region>-summary.json`,
`<region>-replay.json`, `closure-25-279.json`, `run.log`, `replay.log`. The
original workspace paths inside these JSON files are `results/...`. Fresh
outputs go to `--out` (default `scratch/rank4-gap-census/`); `replay` and
`summary` read the preserved files.

    python code/rank4_gap_census.py replay --regions odd-div4-280 even-div4-280 --workers 4 --out scratch/rank4-gap-census
    python code/rank4_gap_census.py replay --workers 16 --out scratch/rank4-gap-census
    python code/rank4_gap_census.py summary --out scratch/rank4-gap-census
    python code/rank4_gap_census.py run --workers 16 --out scratch/rank4-gap-census

The first command replays the determinant-280 chunks (about a minute, also
run in CI); the second replays all 239 closure chunks (about six minutes on
16 cores); `summary` re-derives the closure from the preserved summaries and
replays, re-hashes the prior lists and checks the even-determinant mod-4
fact; `run` regenerates every chunk from scratch (about four minutes on 16
cores; the chunk bytes differ from the preserved ones only in timing fields).

## 7. Consequences for existing statements

- The bracket 24 < det <= 280 of release 1.1.1 is an equality; the
  follow-up note states it as Theorem D, and the determinant-280 lattice as
  Theorem A''.
- The dated correction notes of 2026-10-21 in the level-four, Eisenstein and
  literature reports and in the reviews are completed with this result.
- The global-minimum question of the deferred-research table is answered.
