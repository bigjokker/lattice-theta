# Rank-one/two census, determinant 1 through 24

Completed initial run and separate closure run 2026-09-27 following
the reviewed blueprint (historical workspace reference). Lattice generation, exact
deduplication, full automorphism groups, minima and characteristic verdicts are
complete for ranks one/two and determinant 1 through 24.

## Final classification and closure run

Subsequent uniform result: [docs/rank-two-classification.md](rank-two-classification.md) now proves the binary
classification and symmetry converse in ranks one/two at every determinant.
The finite census below is preserved as independent computational evidence;
its bound and historical verdicts are unchanged.

| Rank | Isometry classes | Even pairs | Proved zero | Proved nonzero | Unresolved | Odd zero by parity |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 24 | 72 | 0 | 72 | 0 | 24 |
| 2 | 70 | 700 | 44 | 656 | 0 | 420 |
| Total | 94 | 772 | 44 | 728 | 0 | 444 |

The separate closure run uses norm4 bound 25, with the same domain and node limit
1,000,000. Its [pack](../data/census/rank12-det24-closed.json.gz),
[box report](../reports/census-rank12-det24-closed-box.json) and
[LDL report](../reports/census-rank12-det24-closed-ldl.json) resolve every pair.
Both backends agree; no computation reaches a resource limit. The original run
and its reports are preserved without changes. A regression checks that its 739
resolved verdicts and evidence, all groups and all stabilizers remain the same.
The original standalone verifier also rechecks all 33 new nonzero claims.

**Finite-domain theorem.** Every rank-one class satisfies the parity converse.
In rank two, the property holds exactly for the 26 non-diagonal reduced
representatives and fails for the 44 diagonal representatives. Equivalently,
within this finite domain it holds exactly for orthogonally indecomposable
rank-two lattices: an orthogonal decomposition into integral rank-one lattices
would give a diagonal representative already present in the generated domain,
and exact pairwise nonisometry rules this out for the non-diagonal classes.
There are 50 positive and 44 negative lattice classes in total.

The symmetry converse holds for every characteristic in this finite domain:
each even pair is proved nonzero or has a verified epsilon-one symmetry witness;
odd pairs have the witness -I by the parity theorem. This is a finite-domain
conclusion, not the unrestricted symmetry converse or an all-determinant binary
classification.

The formerly unresolved pairs have these complete first-shell coefficients:

| Gram matrix | a | b | norm4 | Signed coefficient | Number |
|---|---|---|---:|---:|---:|
| [d], d=17,...,24 | (1) | (0) | d | 2 | 8 |
| diag(1,c), c=16,...,24 | (1,1) | (0,0) | 1+c | 4 | 9 |
| diag(1,c), c=17,...,24 | (0,1) | (0,0) or (1,0) | c | 2 | 16 |

Twenty-five have b=0 and are already nonvanishing by docs/foundations-and-exact-enumeration.md section 8.
For the other eight, y=(0,+/-1) is the complete minimal shell. Its norm4 is c;
x_1=0, so b=(1,0) gives sign +1 for both vectors. The closure run supplies shell
certificates for all 33 uniformly, without adding a new identity-certificate type.

## Preserved initial run, norm4=16

The prescribed initial run had complete lattice, group and characteristic
coverage but partial theta resolution. No rank, determinant or bound expansion
was used in that run. Its unchanged results are:

| Rank | Generated matrices / isometry classes | Even pairs | Proved zero | Proved nonzero | Unresolved pairs | Odd zero by parity |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 24 / 24 | 72 | 0 | 64 | 8 | 24 |
| 2 | 70 / 70 | 700 | 44 | 631 | 25 | 420 |
| Total | 94 / 94 | 772 | 44 | 695 | 33 | 444 |

All generated matrices are pairwise nonisometric; the exact deduplication search
discards none in this domain. Every even zero has an exact symmetry witness.
Full automorphism orders are 2 (31 lattices), 4 (57), 8 (4), and 12 (2).
Each even pair records its complete stabilizer as indices into the full group,
and epsilon for every stabilizing element.

All 70 rank-two lattice-level parity-converse verdicts are settled: 26 positive,
44 negative. The negative classes are exactly the generated diagonal forms,
each with one even zero. The 25 unresolved rank-two pairs occur in nine already
negative lattices and do not obstruct those lattice-level verdicts. In rank one,
the bounded shell run proves the property for d=1,...,16 and leaves d=17,...,24
unresolved. The established scaling theorem and A1 proof already imply positivity
for all rank-one lattices. These eight unresolved entries describe the run's
bounded evidence, not new mathematical uncertainty.

The saved report therefore counts 42 positive lattices by complete shell evidence,
44 negative by zero witnesses and eight unresolved by bounded evidence. Applying
the existing rank-one theorem gives the finite-domain parity classification:
50 positive and 44 negative classes.

No symmetry-converse counterexample was certified. The unresolved pairs prevented
the initial run alone from proving that converse throughout the domain. Their complete
stabilizers all have epsilon zero, but their empty enumerated spheres do not prove
identical vanishing.

## Initial unresolved pairs, now closed above

Binary a is primal and b dual, as in [docs/foundations-and-exact-enumeration.md](foundations-and-exact-enumeration.md). All these searches
completed through norm4=16 without reaching a resource limit. Every sphere is
empty because the bound lies below the first vector.

| Gram matrix | a | b | Number |
|---|---|---|---:|
| [d], d=17,...,24 | (1) | (0) | 8 |
| diag(1,16) | (1,1) | (0,0) | 1 |
| diag(1,c), c=17,...,24 | (0,1) | (0,0) or (1,0) | 16 |
| diag(1,c), c=17,...,24 | (1,1) | (0,0) | 8 |

The pack and both reports retain every pair separately. No unresolved entry is
treated as a numerical or proved zero.

## Coverage and completeness proofs

Rank one consists of G=[d] for 1<=d<=24. In rank two, subtract multiples of the
first basis vector from the second until |2b|<=a and interchange the vectors when
c<a. Each interchange strictly decreases the positive integer a, so reduction
terminates. A sign change gives 1<=a<=c and 0<=2b<=a. All operations are unimodular.
Since ac-b^2>=3a^2/4, the determinant bound gives a<=5. Enumerating these a,b and
a<=c<=floor((24+b^2)/a), keeping determinants 1 through 24, covers every class.
Odd, even, nonprimitive and decomposable forms are included at their supplied
scale, with no minimum restriction. Exact rational LDL certifies admission.

Every isometry U with U^t H U=G sends each source basis vector to an integer vector
in H with squared norm <=max(G_ii). Enumerate that complete sphere, form every
ordered tuple of images with the required norms, and test the Gram identity and
determinant +/-1. This proves both positive and negative isometry results.
Taking H=G enumerates the full automorphism group. The sphere of bound G_00
contains a nonzero shortest vector, since it contains the first basis vector,
and hence also determines the exact minimum.

The builder uses the project's rational LDL sphere enumeration. Default replay
independently enumerates rectangular boxes with the exact inverse-Gram bound
|y_i|^2 <= (y^t G y)(G^-1)_ii, integer square-root bounds and direct norm tests.
It does not call LDL enumeration or the builder's isometry/classification routines.
Replay independently regenerates the finite domain, checks all saved basis changes,
proves retained classes pairwise nonisometric, recomputes full groups and minima,
and checks every stabilizer and epsilon using inverse transpose. Symmetry witnesses
also pass the original exact verifier. It independently sums each complete shell
histogram. Optional LDL shell replay gives a second report; all per-pair results agree.

Rank n has 2^(n-1)(2^n+1) even pairs. Replay requires every one exactly once, in
strict binary-coordinate order. The other 4^n pairs vanish by the established
odd-parity theorem. Gram hashes identify matrices and bases; exact isometry
searches establish class uniqueness separately.

Each sphere and basis-image search has its own default node limit 1,000,000.
Limits discard partial enumerations. An incomplete isometry comparison retains
the candidate with uniqueness unresolved. An incomplete group search publishes
no partial group or negative witness claim. An incomplete shell search publishes
unresolved evidence without partial coefficients. The saved run hit no limits.

## Artifacts and rerun commands

- [code/lattice_census.py](../code/lattice_census.py): generation, exact deduplication, full
  groups, stabilizers and three-way theta classification.
- [code/verify_lattice_census.py](../code/verify_lattice_census.py): independent box replay
  and optional LDL shell replay.
- [Initial certificate pack](../data/census/rank12-det24.json): generation map, all
  matrices, group lists, minima and characteristic evidence.
- [Box report](../reports/census-rank12-det24-box.json) and
  [LDL report](../reports/census-rank12-det24-ldl.json): separate coverage summaries,
  all replay results, unresolved pairs and pack hashes.
- [Closure pack](../data/census/rank12-det24-closed.json.gz) and its
  [box](../reports/census-rank12-det24-closed-box.json) /
  [LDL](../reports/census-rank12-det24-closed-ldl.json) reports: complete theta resolution.
- [tests/test_lattice_census.py](../tests/test_lattice_census.py): ten regressions covering
  noncanonical bases, insufficient invariant matching, missed stabilizing products,
  malformed data, resource limits, the ternary first-shell trap and CLI exit codes.

```sh
python scripts/prepare_data.py
python code/verify_lattice_census.py data/census/rank12-det24.json
python code/verify_lattice_census.py data/census/rank12-det24.json --backend ldl
python code/verify_lattice_census.py data/census/rank12-det24-closed.json
python code/verify_lattice_census.py data/census/rank12-det24-closed.json --backend ldl
python scripts/run_tests.py --all
```

The dedicated census regressions are included in the complete suite. Initial-run replay commands
intentionally exit 2 for unresolved theta pairs; all closure-run commands exit 0.
`certificate_replay_complete=true` with no `replay_limits`
means all published claims were successfully validated; the unresolved verdict
preserves their scientific status. Exit 1 denotes invalid replay input; exit 0
denotes fully resolved valid evidence. Any later bound increase must use separate
output paths to preserve this initial run.
