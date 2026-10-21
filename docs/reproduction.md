# Reproducing the completed results

Run commands from the repository root with Python 3.11 or later. No third-party
packages, network access, or computer-algebra installation is needed.

## Prepare and test

```sh
python scripts/verify_integrity.py
python scripts/prepare_data.py
python scripts/run_tests.py
python scripts/run_tests.py --all
```

The default suite contains 67 rank-four tests (31 from release 1.0.0 and 36
for the level-four, Eisenstein and parity work). The complete suite contains
151 tests covering the original exact verifier, low ranks, root lattices,
modular cutoffs, decomposition, the rank-four counterexample and the 1.1.0
follow-up.

Large JSON files are stored as lossless `.json.gz` files. Preparation restores
their original bytes and verifies the uncompressed SHA-256 hashes. It does not
modify any preserved report. Existing differing files are rejected. The large
preserved outputs of release 1.1.0 (the complete class lists and family
enumerations, about 860 MB expanded) are listed separately in
`reports/compressed-reports.json` and are expanded only on request:

```sh
python scripts/prepare_data.py --reports
```

## Rank-four counterexample

The self-contained [certificate](../data/rank4-counterexample/det23716.json)
has SHA-256:

```text
a1effafe9871bbb833a39a187d98f718879986d24fd23444e167b8635b55e36c
```

```sh
python code/rank4_counterexample.py --replay data/rank4-counterexample/det23716.json --output scratch/counterexample-replay.json
```

This checks the metric identity, complete coset supports, unimodular bijection,
full original automorphism group, normalization, and both stabilizer signs.
The five dedicated adversarial tests include a replay without positive-norm
coefficient corroboration. The [proof](../papers/rank-four-counterexample.md)
establishes identical vanishing at every norm independently of truncation.
Preserved box replays are in [reports/rank4-counterexample](../reports/rank4-counterexample/).

## Rank-four follow-up (release 1.1.0)

The three results of [the follow-up note](../papers/rank-four-follow-up-results.md)
each have a builder, an independent verifier with its own enumeration and
isometry backends, adversarial tests, certificates under `data/` and preserved
outputs under `reports/`. The verifiers never import the builders' enumeration
code.

Certificate hashes (SHA-256):

```text
76175459bf0b317843fae80ac14da5186ca07aa351d3a40a838c522e69b23668  data/rank4-level4/det1729-level4.json
a279081ac204bd75ab1e3694348dc9c6d07465cd0e265fef4c258b3ba577cd38  data/rank4-level4/cs-bijection.json
d11ad820da6ccb19f814f47f755d529dd0e8eb91f322100ff95121e9bbea8a04  data/rank4-eisenstein/det1156-level2.json
56773197994e9be621db875e3aed79f8a92f89147f9aad4b1bfdc5ec9beda282  data/rank4-parity/c1-counterexamples.json
```

Certificate replays (seconds, not minutes; these also run in CI):

```sh
python code/verify_rank4_level4.py --bijection data/rank4-level4/cs-bijection.json --certificates data/rank4-level4/det1729-level4.json --output scratch/level4-det1729-replay.json
python code/verify_rank4_level4.py --bijection data/rank4-level4/cs-bijection.json --certificates data/rank4-level4/other-level4-certificates.json --output scratch/level4-others-replay.json
python code/verify_rank4_eisenstein.py --certificates data/rank4-eisenstein/det1156-level2.json --output scratch/eisenstein-det1156-replay.json
python code/verify_rank4_parity.py certificates data/rank4-parity/c1-counterexamples.json --output scratch/parity-certificates-replay.json
python code/verify_rank4_parity.py pair-members data/rank4-parity/pair-members-parity.json --output scratch/parity-pair-members-replay.json
python code/check_small_pairs.py
```

The level-four replay checks the Conway-Sloane bijection on the nine cosets,
the modular cutoff with complete signed shells by a box enumeration, the
isometries of the coset lattices to the Conway-Sloane pair, and their
non-isometry. The Eisenstein replay checks the order-three kernel isometry,
the cutoff, and the manuscript-schema document through the unchanged
`code/rank4_counterexample.py`. The parity replay checks each certified zero
(cutoff and witness) and each nonzero (complete shell) with a box backend.
`code/check_small_pairs.py` proves the two isospectral pairs of determinants
961 and 3844 equal at every level by the weight-two Sturm bound and confirms
non-isometry with two independent searches.

Complete-list replays (after `python scripts/prepare_data.py --reports`;
times on a 24-core desktop, single-threaded Python):

```sh
python code/verify_rank4_level4.py --task-a reports/rank4-level4/task-a.json --bijection data/rank4-level4/cs-bijection.json --certificates data/rank4-level4/det1729-level4.json --output scratch/level4-task-a-replay.json
python code/verify_rank4_level4.py --task-d reports/rank4-level4/task-d-det24.json --census data/census-rank4/det24-bound32.json --output scratch/level4-task-d-replay.json
python code/verify_rank4_eisenstein.py --enumeration reports/rank4-eisenstein/task-b.json --certificates data/rank4-eisenstein/det1156-level2.json --task-e reports/rank4-eisenstein/task-e.json --output scratch/eisenstein-b-c-e-replay.json
python code/verify_rank4_eisenstein.py --enumeration reports/rank4-eisenstein/task-d.json --output scratch/eisenstein-d-replay.json
python code/verify_rank4_parity.py task-a reports/rank4-parity/task-a.json --output scratch/parity-task-a-replay.json
python code/verify_rank4_parity.py classes reports/rank4-parity/task-b-det201.json --output scratch/parity-b201-replay.json
python code/verify_rank4_parity.py classes reports/rank4-parity/task-c-det64.json --output scratch/parity-c64-replay.json
python code/verify_rank4_parity.py classes reports/rank4-parity/task-e-det168.json --fast-box --output scratch/parity-e168-replay.json
```

The family replay to D' = 154 takes about four minutes; the odd-lattice class
list (3,487 classes) about seven minutes; the even odd-determinant list (300
classes) about two minutes. Each replay regenerates the candidate domain with
its own reduction box, recomputes every isometry and every decision, and
reports `verified` only if every count matches. The reviewers' independent
re-enumerations with a different (Minkowski) box are
`code/rank4_parity_review_enum.py` and `code/rank4_eisenstein_review_enum.py`;
their outputs are preserved under `reports/`.

Independent review checks (replicate the headline numbers from the accepted
modules only; each prints and exits):

```sh
python code/rank4_level4_review_check.py
python code/rank4_eisenstein_review_check.py
python code/rank4_eisenstein_review_enum.py
python code/rank4_parity_review_check.py
python code/rank4_parity_facts_check.py
python code/rank4_klein_triple.py
python code/rank4_eisenstein_formula_check.py
```

The parity builder and verifier resolve the hash-bound source paths stored in
the documents (`examples/...`, `results/...`) to `data/` and `reports/`; this
is the only change from the reviewed research copies (see
[export provenance](export-provenance.json)). The other builders keep their
research-workspace default paths and take explicit paths on the command line;
inspect `--help` before regenerating, and use fresh output paths under
`scratch/`. The regeneration commands of each report are in its section 8,
with paths relocated to this layout.

## Complete finite censuses

```sh
python code/verify_lattice_census.py data/census/rank12-det24-closed.json
python code/verify_lattice_census.py data/census/rank12-det24-closed.json --backend ldl
python code/verify_rank3_census.py data/census-rank3/det24-bound32.json
python code/verify_rank3_census.py data/census-rank3/det24-bound32.json --backend ldl
python code/verify_rank3_structure.py data/census-rank3/structure-det24.json
python code/verify_rank4_census.py data/census-rank4/det24-bound32.json
python code/verify_rank4_census.py data/census-rank4/det24-bound32.json --backend ldl
python code/verify_rank4_feasibility.py data/rank4-feasibility
```

The censuses cover determinant at most 24. The rank-four census has 169
isometry classes and 43,264 binary characteristics, all resolved; all zeros
in that finite domain have a symmetry witness. This does not imply the
unrestricted converse. The determinant-1156 and determinant-23,716
counterexamples are outside it.

The initial rank-one/two census is retained because tests explicitly check
its unresolved records against the completed closure; the rank-four pilot
is a coverage regression fixture. They are not additional final results.

## Root-lattice checks

```sh
python code/verify_even_pack.py data/e6-even.json
python code/verify_e6_ambient.py data/e6-even.json
python code/verify_exceptional.py data/e7-classification.json
python code/verify_exceptional.py data/e7-classification.json --backend ambient
python code/verify_exceptional.py data/e8-classification.json
python code/verify_exceptional.py data/e8-classification.json --backend ambient
python code/verify_dn.py data/dn/d10.jsonl.gz --backend ambient
```

The complete tests replay the D2–D10 streams and the F4/G2 certificates.
Uniform statements rely on the manuscript and proof appendices, rather than
only the saved finite instances.

## Reports, regeneration, and cleanup

Saved reports in `reports/` retain their original bytes, input hashes, and
historical source-path strings. Map an original `examples/` path to `data/`
and a `results/` path to `reports/`. A historical report's review-status field
describes its original run; subsequent review acceptance is recorded separately.
See [export provenance](export-provenance.json) for every source mapping.

Use a fresh path under `scratch/` for new verifier reports. Builders write
certificates or reports; inspect their `--help` before regeneration. There is
no need to regenerate inputs to independently verify them.

To remove only the verified expanded copies before a browser upload:

```sh
python scripts/prepare_data.py --clean
python scripts/prepare_data.py --clean --reports
```

Compressed inputs remain intact. `scratch/`, expanded inputs, and Python caches
are excluded from Git by `.gitignore`. Upload the contents of this repository
directory; do not include the surrounding historical research workspace.
