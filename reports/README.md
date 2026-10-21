# Preserved verification reports

These are selected completed reports from the original research workspace.
Their bytes and original certificate-hash bindings are unchanged. The
rank-four counterexample has independently saved original and Astra box replays;
the census reports include separate box and LDL backends.

Historical paths inside reports use `examples/` and `results/`; their curated
locations are `data/` and `reports/`. Some original status fields predate later
review acceptance. [Review evidence](../reviews/README.md) and the completed
manuscripts state the current conclusions.

Release 1.1.0 adds the rank-four follow-up outputs:

- `rank4-level4/`: the index-two neighbour tables of the nine known
  isospectral pairs (`task-a.json`), the determinant-24 sign-kernel scan
  (`task-d-det24.json`), verifier replays, logs, test output and the
  reviewer's independent check.
- `rank4-eisenstein/`: the family enumerations (`task-b.json` to D' = 41,
  `task-d.json` to D' = 154 at prefix 400, and `task-d-prefix60.json`, the
  first run kept as evidence of the filter failure described in the report),
  the mixed-type lemma cases, replays, logs and the reviewer's checks.
- `rank4-parity/`: the evidence table (`task-a.json`, rendered in
  `task-a.md`), the complete class lists (`task-b-det49/101/201.json`,
  `task-c-det64.json`, `task-e-det168.json`), the structure analysis
  (`task-d.json`), the E8 check, replays, logs, the control facts and the
  reviewer's independent re-enumerations.
- `project-review-2026-10-21/`: the modular-form check of the two small
  isospectral pairs from the independent assessment.

Preserved outputs above one million bytes are stored as `.json.gz` and listed
in [compressed-reports.json](compressed-reports.json);
`python scripts/prepare_data.py --reports` expands them (about 860 MB) and
verifies their hashes. They are inputs to the verifier replays only; the tests
do not need them.

Fresh outputs belong under `scratch/`. See [reproduction](../docs/reproduction.md).

[Curated export validation](export-validation.json) records the fresh 115-test
pass, standalone module imports, and counterexample and full rank-four census
replays performed after the 1.0.0 reorganization;
[export-validation-1.1.0.json](export-validation-1.1.0.json) records the
checks performed for release 1.1.0.
