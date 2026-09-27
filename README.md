# Lattice Theta Functions with Half-Integral Characteristics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23002464.svg)](https://doi.org/10.5281/zenodo.23002464)

Exact proofs, certificates, and reproducible computations for identical
vanishing of lattice theta functions with half-integral characteristics.

The symmetry converse holds in ranks one through three and fails in rank four.
The rank-four counterexample is a primitive, positive definite integral lattice
of determinant **23,716**. An even characteristic has identically zero theta
series, while the complete automorphism group is **{I, −I}**, with trivial sign
on both elements. The proof uses an all-level coset bijection and a complete
315-point enumeration of possible automorphism columns.

The completed work also includes the low-rank parity classification, the
classification of irreducible crystallographic root lattices in the stated
normalization, and exact determinant-at-most-24 censuses through rank four.
No smallest-counterexample-determinant or literature-priority claim is made.

## Read the results

- [Rank-four counterexample manuscript](papers/rank-four-counterexample.md)
- [Low-rank and root-lattice manuscript](papers/low-rank-and-root-lattice-results.md)
- [Reproduction instructions](docs/reproduction.md)
- [Mathematical appendices](docs/README.md)
- [Review evidence and its scope](reviews/README.md)

## Reproduce

Use **Python 3.11 or later**. All computation uses the Python standard library;
there are no third-party dependencies. Run these commands from the repository root.

```sh
python scripts/verify_integrity.py
python scripts/run_tests.py
python code/rank4_counterexample.py --replay data/rank4-counterexample/det23716.json --output scratch/counterexample-replay.json
```

The default suite runs the 31 rank-four tests. The output path for a replay
must be fresh. The test runner automatically prepares compressed datasets.
To run the complete 115-test suite:

```sh
python scripts/run_tests.py --all
```

The written coset bijection proves identical vanishing at every level; a long
list of zero coefficients is not used as the identity proof. Independent box
and LDL enumerations supply separate checks of the finite computational evidence.

## Repository layout

| Directory | Contents |
|---|---|
| `papers/` | Completed manuscripts |
| `docs/` | Proof appendices, source scope, reproduction, and export provenance |
| `code/` | Exact certificate builders and verification programs |
| `tests/` | Coverage, corruption, normalization, and resource-limit regressions |
| `data/` | Exact certificates and compressed census inputs |
| `reports/` | Preserved independent verification reports |
| `reviews/` | AI-assisted local mathematical reviews |
| `scripts/` | Dataset preparation, tests, and file-integrity verification |

This is a curated research release. Exploratory scripts, duplicate snapshots,
handoff prompts, and downloaded third-party papers are excluded. Certificate
schemas retain their original identifiers for compatibility. Scientific work
through the rank-four counterexample is complete; [further questions](docs/deferred-research.md)
are recorded separately and are not active assignments.

## Citation

Archived on Zenodo. Cite the concept DOI to refer to the work in general, or
the version DOI to pin a specific release.

* All versions: [10.5281/zenodo.23002464](https://doi.org/10.5281/zenodo.23002464)
* v1.0.0: [10.5281/zenodo.23002465](https://doi.org/10.5281/zenodo.23002465)

```bibtex
@misc{lattice_theta,
  title  = {Lattice Theta Functions with Half-Integral Characteristics},
  author = {bigjokker},
  year   = {2026},
  doi    = {10.5281/zenodo.23002464},
  url    = {https://github.com/bigjokker/lattice-theta}
}
```

This repository is released under the MIT License. See [LICENSE](LICENSE).
The local reviews are not a claim of journal publication or formal peer review.
