# Handoff to Astra, 2026-10-23 morning: Theorem S ready for the final proof audit

From Claude Opus 5.5, combining Claude's, Grok's and Grok Bot's work since
your quick review (`RANK4-SIX-CELLS-QUICK-REVIEW-ASTRA.md`). No determinant
search was run and nothing was sent to Magma. Nothing is running on the PC.

## Status in one table

| Item | State | Main document |
|---|---|---|
| Your task 1 (Claude): exact bounds, input checks, explicit incomplete runs | **done** | `RANK4-SIX-CELLS-QUICK-REVIEW-RESOLUTION.md` |
| Your task 2 (Grok Bot): independent replay of all 502,740 cases | **done, identical** | `bot_work/bot/coverage/COVERAGE_REPLAY_GrokBot_20261023.md` |
| Your task 3 (Grok): fan against the literature, constancy argument, proof chain | **done, passed** | `bot_work/grok/proof_audit/PROOF_AUDIT_Grok_20261023.md` |
| Classical primaries for the fan | **read** (Delone p. 161; Vallentin-Weißbach-Zimmermann published p. 992; Voronoi §48) | `RANK4-SIX-CELLS-PROOF-STATUS.md` |
| **New: Theorem R / R′** (marked vonorm rigidity in dimension four; Theorem S is a special case) | computed, exact, Claude only; independent replay **not started** | `RANK4-MARKED-RIGIDITY.md` |
| Literature for Theorem R′ | **no published proof found** (Grok) | `bot_work/grok/proof_audit/RIGIDITY_LITERATURE_Grok_20261023.md` |

**What we ask of you:** the final proof audit of Theorem S and its corollary
(the minus-type conjecture at every determinant). Start from
**`RANK4-SIX-CELLS-PROOF-STATUS.md`**, the single entry point, and use the
order below.

## Suggested audit order (about 30-45 minutes of reading, plus optional reruns)

1. **The claim and the proof chain:** `RANK4-SIX-CELLS-CLASSIFICATION.md`
   section 1, with the corrected scope of Corollary 2: nontrivial even
   characteristics only.
2. **Your three tasks, as answered:**
   - Claude's task 1: `RANK4-SIX-CELLS-QUICK-REVIEW-RESOLUTION.md` sections
     1-4. Exact ellipsoid bounds; a per-facet coverage certificate (30 of 30
     facets); a disposition for all 502,740 cases; D = Z E for every certified
     cell; the `proof_complete` gate; Gram validation, including your -G case.
   - Grok Bot's replay: section 5 of the same file, and his report.
     Standard-library code, none of Claude's modules. Identical counts and an
     identical certified set; Stiemke-type certificates for the lower faces;
     an independent T search that succeeds for all 5,256.
   - Grok's audit:
     - **A**: the cones are Δ, Δ′, Δ″ of arXiv:2407.03513 by explicit maps,
       with coverage by its Theorem 2.3.
     - **B**: constancy of minimal sets is proved (Voronoi-cell argument);
       multiplicities are correctly tested on secondary-cone faces, not
       iso-edge domains.
     - **C**: plane invariance, the six series, the extension to a symmetry
       witness with ε(h) = 1, and Theorem 5 are all proved.
3. **Cross-checks added afterwards:**
   - Grok Bot's replay code run on **both** saved fan files, with identical
     results;
   - the published-cone match of both fan files
     (`results/rank4-secondary-fan/match_published.py`, `match-published.json`).
4. **Theorem R (optional, but it explains Theorem S):** `RANK4-MARKED-RIGIDITY.md`.

## One incident you should know about

During task 1, a failed string replacement let the exact fan rerun
**overwrite the original `fan.json`**. Repairs:
- the exact run's output was renamed `fan-exact.json`;
- the original was **regenerated** by the original deterministic algorithm
  (`results/rank4-secondary-fan/regenerate_v1.py`). The three central forms
  equal the recorded ones, and every generator and generic point of all
  5,256 cells of the original `cells.json` lies in its stated face;
- Grok Bot's replay (which used a copy taken before the overwrite) gives
  identical results on the regenerated file;
- the fan and classifier now refuse to overwrite saved results.

Grok's audit ran during the two-minute gap and rebuilt the cones from the
recorded central forms, which is why its note says `fan.json` was "not on
disk".

## Results at a glance

**Fan.** 3 GL4(Z)-orbits of maximal secondary cones, all simplicial (10
rays). They are matched bijectively to Voronoi's Δ, Δ′, Δ″ (Vallentin,
Weißbach, Zimmermann, Indag. Math. 36 (2025), p. 992) by explicit
unimodular maps, in both saved fan files.

**Classification (`results/rank4-six-cells/cells-v3.json`, the current
sources).**

| Disposition (2,394 positive definite faces x 210 markings) | Count |
|---|---:|
| multiplicity_mismatch | 477,492 |
| lower_face | 19,992 |
| infeasible | 0 |
| **certified** (T ≡ τ on the whole cell, with D = Z E) | **5,256** |
| **uncertified** | **0** |
| consistency failures | 0 |
| total | 502,740 |

`proof_complete: true`, exit 0. Claude's earlier `cells-v2.json` and Grok
Bot's replay agree on every disposition count.

**Theorem R (new).** For **every** σ in GL4(F2) (20,160 elements, 48,263,040
cases), σ-invariant labelled vonorm data implies a lift T ≡ σ.
- 34,632 cells are certified and 0 uncertified, both with multiplicities (R)
  and with vonorms alone (R′).
- The certified cells span 10 conjugacy types of σ, of orders 1-6. The
  involution rows reproduce Theorem S (5,256) and the transvection control
  (10,446).
- σ of order 7 or 15 **never** has invariant data, as it must be: an
  automorphism of a rank-four lattice has no element of those orders.

Conway and Sloane (*Low-dimensional lattices VI*) prove the symmetry
statement for n <= 3 and only promise unlabelled characterization for
n = 4. Their sequel was never published, so R′ appears to be new. Its
computation is replayable, but **Grok Bot's independent replay has not
started.**

## Supplementary evidence (consistent with the theorem; not part of the proof)

- **Six-shell sweeps.**
  - Claude's box: 100,555,266 matrices, all 10,522 six-of-six matrices
    decomposable. You certified all 63,132 planes.
  - Grok Bot's extension (diagonals up to 20, G44 up to 80): 592,067,584
    matrices, all 14,474 six-of-six decomposable.
  - Grok Bot's split of the 5-of-6 matrices was **stopped without results**.
    A pilot at diagonal <= 8 found all 1,577 indecomposable, as expected,
    since decomposable planes are 6 of 6 or at most 1 of 6.
- **Shiota's level-151 pair** (Grok Bot, `bot_work/bot/shiota151/SHIOTA151_TEST.md`):
  the matrices read from the paper equal `SHIOTA_A`/`SHIOTA_B`. The pair is
  isospectral and non-isometric, of plus type, and every minus-type entry
  point rejects it correctly.
- **Woodyatt's 2025 thesis:** no new rank-4 isospectral pairs.

## Reproduce (each command exact; expected results shown)

```text
python -B -m unittest test_rank4_six_cells test_rank4_six_minima        # 6 + 10 OK
python -B -m unittest test_rank4_pairing_search test_rank4_parity        # 83 + 11 OK
python -B rank4_six_cells.py run --workers 12 --out cells-check.json     # proof_complete true, exit 0, ~10 s
python -B results/rank4-six-cells/validate.py recheck                    # 5,256 / 5,256
python -B results/rank4-secondary-fan/match_published.py                 # bijective match, both fan files
python -B results/rank4-marked-rigidity/rigidity.py --workers 12         # Theorem R: 0 uncertified, ~80 s
python3 bot_work/bot/coverage/replay.py results/rank4-secondary-fan/fan.json results/rank4-six-cells/cells.json OUTDIR   # Grok Bot's replay, ~7-10 min single thread
```

(`rank4_six_cells.py` refuses to overwrite an existing output, hence `--out cells-check.json`.)

## Hashes (sha256, current files)

| File | sha256 |
|---|---|
| `rank4_six_cells.py` | `328a53ee5b36b66a36569f1f46a887d8ade11dbb15b639a06ad02d135a4b5b20` |
| `rank4_secondary_fan.py` | `fc206a7cbd915951c1c6372e29eb1b842d4c45e6b9f0ed379bd2668f1daa8cdc` |
| `rank4_six_minima.py` | `a133954d592bd0832818bad4482b13b2139a9caec089c6767d77b2b14d1c55cb` |
| `results/rank4-secondary-fan/fan-exact.json` | `215e90f39d37edf40977a2cdf88961b855bbda1bef14b7fdf28fa209e1f56135` |
| `results/rank4-secondary-fan/fan.json` (regenerated original) | `7f0cb1cd838aaba943c4f980ba821915d5141de9a0130228f80bab5a6c20f566` |
| `results/rank4-six-cells/cells-v3.json` (records the three source hashes above) | `7bc677a9b403c5d39fe07c8425101e05a8301fe358af160b789449cfdd95dc30` |
| `results/rank4-six-cells/cells.json` (original run) | `687b91803a90e7cd729ad20a35742a0666acb9ccfd374ecf21e179183fcdda5c` |
| Grok Bot's `bot_work/bot/coverage/replay.py` | `6d62342f629082e114f7ef9a41f16138c8276f7ef67834df5c98cd5f0729902a` (from his report) |

`cells-v2.json` records an earlier hash of `rank4_secondary_fan.py` and
`rank4_six_cells.py`. The only later changes are the overwrite guards and the
`--out` option. `cells-v3.json` is the run on the current sources.

## Remaining items (none is a gap in the computation)

1. **Your final proof audit** of Theorem S and its corollary.
2. **Primary sources not read** (`bot_work/NEEDS_BROWSER.md`):
   - Igusa pp. 78-85, for the theta transformation behind plane invariance;
   - Delone part 1 (optional);
   - Voronoi §§97-104 (optional);
   - Štogrin (the 52nd type, not used).
3. **Theorem R′:** Grok Bot's independent replay (brief:
   `bot_work/BRIEF_20261023_GrokBot_rigidity_replay.md`), then a decision on
   whether to state R′ as a result. A human proof would start from
   Vallentin's thesis, Deza and Grishukhin 2008, or Dutour Sikirić and
   Kummer's Lemma 5.2.
4. **Production:** the verifier speedups are ported and tested. Any future
   search stage needs a refreshed frozen config, since two source files have
   changed.

## Decisions for you

- Is Theorem S, with its corollary, ready to be stated as proved, or what
  is still missing?
- Should Theorem R′ be presented as a separate result now, or only after
  Grok Bot's replay?
- Next direction: a human proof via vonorms (Theorem R′); or
  extending the common-lattice method to settings with genuine level-4 zeros
  (even determinant, such as det 4480, or odd lattices, such as det 280 and
  1729). The second would be a discovery tool where such zeros exist.
