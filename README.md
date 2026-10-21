# Lattice Theta Functions with Half-Integral Characteristics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23002464.svg)](https://zenodo.org/doi/10.5281/zenodo.23002464)

Exact proofs, certificates, and reproducible computations for identical
vanishing of lattice theta functions with half-integral characteristics.

The symmetry converse holds in ranks one through three and fails in rank four.
Release 1.1.0 (October 2026) adds three rank-four results, each with exact
certificates, independent replays and a local review:

- **A zero whose coset lattices are not isometric.** An odd lattice of
  determinant **1729** has an identically vanishing even characteristic whose
  two coset lattices are the Conway–Sloane isospectral pair L⁻(1,7,13,19),
  L⁺(1,7,13,19), which are not isometric. No single isometry, inside the
  lattice's automorphism group or outside it, explains the vanishing. Eight
  such configurations are certified.
- **A smaller symmetry counterexample.** An even lattice of determinant
  **1156** with automorphism group **{I, −I}** has an identically vanishing
  even characteristic, by the same Eisenstein-integer construction as the
  original determinant-23,716 example. The minimum counterexample determinant
  in rank four satisfies 24 < det ≤ 1156; no minimality claim is made in
  between.
- **The parity property.** The conjecture that an indecomposable rank-four
  lattice whose determinant is not divisible by four has no vanishing even
  characteristic is false, by complete class lists: odd lattices fail first at
  determinant **49**, and even lattices of odd determinant fail first at
  determinant **121**, a Gaussian-integer lattice whose witnesses square to
  −I. An even indecomposable rank-four lattice of odd determinant has no
  involution witness. Whether every such lattice of determinant 5 mod 8 has
  the parity property is open.

The original results remain: the determinant-23,716 counterexample with its
all-level coset bijection and complete automorphism group {I, −I}, the
low-rank parity classification, the classification of irreducible
crystallographic root lattices in the stated normalization, and exact
determinant-at-most-24 censuses through rank four. No literature-priority
claim is made; the isospectral pairs are Schiemann's, Conway–Sloane's and
Earnest–Nipp's, and their non-isometry is Cervino–Hein's theorem.

## Read the results

- [Rank-four follow-up results, release 1.1.0](papers/rank-four-follow-up-results.md)
- [Rank-four counterexample manuscript](papers/rank-four-counterexample.md)
- [Low-rank and root-lattice manuscript](papers/low-rank-and-root-lattice-results.md)
- [Reproduction instructions](docs/reproduction.md)
- [Mathematical appendices and full reports](docs/README.md)
- [Review evidence and its scope](reviews/README.md)
- [Changelog](CHANGELOG.md)

## Reproduce

Use **Python 3.11 or later**. All computation uses the Python standard library;
there are no third-party dependencies. Run these commands from the repository root.

```sh
python scripts/verify_integrity.py
python scripts/run_tests.py
python code/rank4_counterexample.py --replay data/rank4-counterexample/det23716.json --output scratch/counterexample-replay.json
python code/verify_rank4_level4.py --bijection data/rank4-level4/cs-bijection.json --certificates data/rank4-level4/det1729-level4.json --output scratch/level4-replay.json
python code/verify_rank4_eisenstein.py --certificates data/rank4-eisenstein/det1156-level2.json --output scratch/eisenstein-replay.json
python code/verify_rank4_parity.py certificates data/rank4-parity/c1-counterexamples.json --output scratch/parity-replay.json
```

The default suite runs the 67 rank-four tests. The output path for a replay
must be fresh. The test runner automatically prepares compressed datasets.
To run the complete 151-test suite:

```sh
python scripts/run_tests.py --all
```

The written coset bijections prove identical vanishing at every level; a long
list of zero coefficients is not used as the identity proof. Independent box
and LDL enumerations supply separate checks of the finite computational
evidence. The complete class lists behind the parity and Eisenstein results
are preserved in compressed form and can be replayed with the verifiers; see
the [reproduction guide](docs/reproduction.md).

## Repository layout

| Directory | Contents |
|---|---|
| `papers/` | Completed manuscripts and the 1.1.0 results note |
| `docs/` | Proof appendices, full research reports, source scope, reproduction, and export provenance |
| `code/` | Exact certificate builders, verification programs, and independent review checks |
| `tests/` | Coverage, corruption, normalization, and resource-limit regressions |
| `data/` | Exact certificates and compressed census inputs |
| `reports/` | Preserved outputs, verifier replays, and compressed class lists |
| `reviews/` | AI-assisted local mathematical reviews and the project assessment |
| `scripts/` | Dataset preparation, tests, and file-integrity verification |

This is a curated research release. Exploratory scripts, duplicate snapshots,
handoff prompts, and downloaded third-party papers are excluded. Certificate
schemas retain their original identifiers for compatibility. The original
scientific goal is complete; [further questions](docs/deferred-research.md),
including the directions proposed by the October 2026 assessment, are
recorded separately and are not active assignments.

## Citation

Archived on Zenodo. Cite the concept DOI to refer to the work in general, or
the version DOI to pin a specific release.

* All versions: [10.5281/zenodo.23002464](https://zenodo.org/doi/10.5281/zenodo.23002464)
* v1.0.0: [10.5281/zenodo.23002465](https://zenodo.org/doi/10.5281/zenodo.23002465)
* v1.1.0: [10.5281/zenodo.23540724](https://zenodo.org/doi/10.5281/zenodo.23540724)

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
