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

The default suite contains 31 rank-four tests. The complete suite contains
115 tests covering the original exact verifier, low ranks, root lattices,
modular cutoffs, decomposition, and the rank-four counterexample.

Large JSON files are stored as lossless `.json.gz` files. Preparation restores
their original bytes and verifies the uncompressed SHA-256 hashes. It does not
modify any preserved report. Existing differing files are rejected.

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
unrestricted converse. The determinant-23,716 counterexample is outside it.

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
```

Compressed inputs remain intact. `scratch/`, expanded inputs, and Python caches
are excluded from Git by `.gitignore`. Upload the contents of this repository
directory; do not include the surrounding historical research workspace.
