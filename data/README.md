# Exact input certificates

All certificates preserve the original research input bytes. JSON inputs above
one million bytes are stored as `.json.gz`; run `python scripts/prepare_data.py`
from the repository root to expand them and verify their uncompressed hashes.
The Dn certificate streams were already distributed in gzip form.

- `rank4-counterexample/`: determinant-23,716 all-level counterexample.
- `rank4-level4/`: the determinant-1729 level-four certificate, the
  Conway-Sloane bijection certificate, and the seven other level-four
  configurations (release 1.1.0).
- `rank4-eisenstein/`: the determinant-1156 level-two counterexample, with
  the manuscript-schema document embedded (release 1.1.0).
- `rank4-parity/`: the six determinant-121 certificates of the refuted even
  parity conjecture and the parity certificates of the eighteen isospectral
  pair members (release 1.1.0).
- `census/`, `census-rank3/`, `census-rank4/`: finite determinant-24 domains.
- `rank4-feasibility/`: domain-generation and modular-cutoff controls.
- `rank4-remainder/`, `rank4-fibre-obstruction/`: exact cancellation regressions.
- `dn/`, `root-lattices/`, and the E6/E7/E8 packs: root-lattice evidence.
- Small standalone files: normalization, symmetry-witness, and nonzero-shell controls.

Original schema names are retained to keep certificates compatible with the
verifiers. [Compressed file metadata](compressed-files.json) and
[export provenance](../docs/export-provenance.json) record sizes and hashes.
The complete class lists that the parity and Eisenstein results rest on are
preserved outputs, kept under [reports](../reports/README.md).
