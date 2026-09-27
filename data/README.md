# Exact input certificates

All certificates preserve the original research input bytes. JSON inputs above
one million bytes are stored as `.json.gz`; run `python scripts/prepare_data.py`
from the repository root to expand them and verify their uncompressed hashes.
The Dn certificate streams were already distributed in gzip form.

- `rank4-counterexample/`: determinant-23,716 all-level counterexample.
- `census/`, `census-rank3/`, `census-rank4/`: finite determinant-24 domains.
- `rank4-feasibility/`: domain-generation and modular-cutoff controls.
- `rank4-remainder/`, `rank4-fibre-obstruction/`: exact cancellation regressions.
- `dn/`, `root-lattices/`, and the E6/E7/E8 packs: root-lattice evidence.
- Small standalone files: normalization, symmetry-witness, and nonzero-shell controls.

Original schema names are retained to keep certificates compatible with the
verifiers. [Compressed file metadata](compressed-files.json) and
[export provenance](../docs/export-provenance.json) record sizes and hashes.
