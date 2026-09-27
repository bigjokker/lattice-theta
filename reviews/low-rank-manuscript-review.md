# Final reader and reproducibility review

Prepared 2026-09-27 from `MANUSCRIPT.md` and `MANUSCRIPT-SUPPLEMENT.md`, compared with `PROOF-STATUS.md`, `docs/foundations-and-exact-enumeration.md` §2 and §15, `docs/root-lattice-classification.md`, `docs/rank-one-and-two-census.md`, `docs/literature-and-source-scope.md`, `docs/theta-transformations.md`, `docs/coset-modularity.md`, and `README.md`. `docs/foundations-and-exact-enumeration.md` §15 cites `docs/rank-two-classification.md` for the expanded binary proof; this pass used the theorem and witness statement in §15. The manuscript, proof notes, code, manifest, certificate inputs, saved reports, and census packs were left unchanged. The completed binary, `GL(n,ℤ)`, and rescaling audits are not reopened. This pass found their incorporated statements still in agreement with the current manuscript.

**Checked** means a file read, a byte hash, or a saved field. **Deduction** means a comparison of those texts. **Unrun** means a classification or test command that this pass did not execute.

## 1. Issue table

| Location | Statement | Severity | Reasoning | Correction |
|---|---|---|---|---|
| `MANUSCRIPT-SUPPLEMENT.md`, finite-classification table, historical census row; replay section, the sentence “Replaying the preserved bound-16 census intentionally returns 2.” | The table names `data/census/rank12-det24.json` and the box and LDL reports. The displayed command block has no command for that pack. | reproducibility gap | A reader can name the input and both reports. The print-only block gives E6, E7, E8, the closed census, the 64-test command, Dn ranks 2 through 10, and the F4/G2 one-liner. The historical replay is the finite assertion whose command is missing from that block. `code/verify_lattice_census.py` maps the verdict `unresolved` to exit code 2, and both saved bound-16 reports record that verdict. This pass did not execute the replay. | Add these lines to the print-only block, with no `--output`: `python code/verify_lattice_census.py data/census/rank12-det24.json --backend box` and the same command with `--backend ldl`. Keep them beside the exit-code sentence. |
| `MANUSCRIPT.md`, “Supplement and reproducibility,” the sentence that `README.md` supplies the individual replay commands for the preserved initial and closed census runs | The manuscript sends the reader to the README command block for those preserved runs. | reproducibility gap | `README.md` lines that rebuild the census call `code/lattice_census.py --output` on `data/census/rank12-det24.json` and `data/census/rank12-det24-closed.json`. The following four lines call `code/verify_lattice_census.py --output` on the four saved census reports. Those commands replace the preserved files. The supplement’s command block omits `--output` and prints a new summary. The same README block writes the exceptional reports through `--output` as well. | Point preservation-safe replay to `MANUSCRIPT-SUPPLEMENT.md`. In the same sentence, say that the README block also contains the builders and the `--output` paths above, and that those paths replace the preserved census inputs and reports. |
| `MANUSCRIPT.md`, Proposition 5 proof, “in odd rank the symbol changes by `(8u/v)/(u/v)=(8/v)`” | The sketch records the character factor and then points to `docs/coset-modularity.md` for exact branch conventions. | optional editorial | The proposition statement already gives `H = Γ_0(16N_{L_0}) ∩ Γ^0(8)`, character `χ` for even rank and `χ χ_8` for odd rank, cusp width eight, and the width-eight parameter `exp(2π i τ/8)`. The proof also requires `d_c > 0` and the factor `(d_c/8)^{-n/2}`. Those statements match `docs/coset-modularity.md` (3), (4), and (5). The principal branch, negative `v`, and the `-I` check are written in that note and are omitted from the sketch. | Optional. Copy the note’s paragraph on `Arg ∈ (−π, π]`, negative `v`, and `-I` into the proposition proof. The stated law can stay as it is. |
| `reports/root-lattices-audit.json`, the F4 and G2 objects; `reports/manuscript-artifacts.json`, `report_bindings` | The supplement names the F4 and G2 packs and this audit file as their saved report. | optional editorial | Those two objects record `coverage_complete`, 136 even pairs with 9 zeros and 127 nonzeros, and 10 even pairs with 10 nonzeros. They contain no certificate hash. The 28 report bindings cover E6, E7, E8, both censuses, and Dn ranks 2 through 10. They include no F4 or G2 row. The inventory still hashes `data/root-lattices/f4-classes.json`, `g2-classes.json`, and the audit file as separate files. The supplement’s `verify_root_pack` command is present. | Optional. On a later intentional inventory regeneration, store each pack hash on those two audit objects and add the two bindings. Until then, one sentence in the supplement can say that those rows are count summaries and that pack identity is the inventory hash together with `verify_root_pack`. |

No row is an incorrect theorem, count, character, or citation. The first two rows are the reader-index edits that make every finite assertion in the supplement copyable without replacing a preserved file.

## 2. Agreement of the mathematical statements

Checked against `PROOF-STATUS.md`, `docs/foundations-and-exact-enumeration.md`, `docs/root-lattice-classification.md`, and `docs/rank-one-and-two-census.md`.

Proposition 1 states that epsilon is a homomorphism to `F_2` and that a nontrivial value forces identical vanishing. The coordinate formula is `a^t ((T^{-t}-I)b)/2` modulo 2. The manuscript assigns `-I` the value `-p` modulo 2, which is the same residue as `p` modulo 2, and uses it for odd-p vanishing. `docs/foundations-and-exact-enumeration.md` §2 states the same homomorphism, the same coordinate formula, and the same value of epsilon at `-I`.

Proposition 2 states that every nontrivial orthogonal sum of positive-rank integral lattices fails `P(L)`. The proof builds one odd characteristic on each of two summands. That is the argument indexed in `PROOF-STATUS.md` for `docs/foundations-and-exact-enumeration.md` §7.

Theorem 3 states that a positive definite integral rank-two lattice has `P(L)` exactly when it is orthogonally indecomposable, that a decomposable rank-two lattice has exactly one even zero, and that the symmetry converse holds in ranks one and two at every determinant. The reduced Gram is `[[A,B],[B,C]]` with `1 ≤ A ≤ C` and `0 ≤ 2B ≤ A`. The three nontrivial shells are the coefficients `+2`, `+2`, and `-2` at norm4 `A`, `C`, and `A+C-2B`. The even zero on a diagonal Gram is witnessed by `T = diag(-1,1)`, with dual shift `(-1,0)`. `docs/foundations-and-exact-enumeration.md` §15 records that same witness. Rank one is the constant coefficient 1 for `a = 0` and the coefficient 2 for `a = 1`, `b = 0`, with `-I` on the odd pairs.

Theorem 4 uses short roots of squared length 2. The positive types are `A_n`, `E_6`, `C_3`, and `G_2`, with the rank-one labels `B_1` and `C_1` included as `A_1`. The proof identifies `B_n` with an orthogonal sum of `n` copies of `A_1`, `C_n` with `D_n` for `n ≥ 2`, `F_4` with `D_4`, and `G_2` with `A_2`. `C_3` is positive and `C_n` for `n ≥ 4` is negative. The F4 and G2 packs are the complete sets of 136 and 10 even pairs. `docs/root-lattice-classification.md` states the same normalization, the same positive list, and the same 9 and 0 even-zero counts for `F_4` and `G_2`. Positive rescaling is stated to preserve the parity verdicts in both documents.

The exceptional table in the manuscript is:

| Lattice | Even pairs | Symmetry zeros | Nonzero shells | Squared-norm bound |
|---|---:|---:|---:|---:|
| E6 | 2080 | 0 | 2080 | 1 |
| E7 | 8256 | 1260 | 6996 | 3/2 |
| E8 | 32896 | 9450 | 23446 | 1 |

These are the counts in `PROOF-STATUS.md` and in the supplement. The even-pair totals match `2^{n-1}(2^n+1)`. The manuscript limits the exceptional certificates to the parity verdicts and the recorded witnesses. The supplement says the root packs certify witnesses and that full exceptional automorphism groups are outside those packs. `PROOF-STATUS.md` still lists the full automorphism group of `D_4`, of order 1152, as uncertified. The manuscript does not claim that group.

Section 5 states a domain of 94 classes, 24 of rank one and 70 of rank two, with the closed pack resolving 772 even pairs as 728 nonzero and 44 symmetry zeros, and with 50 positive and 44 negative classes. It states that each diagonal rank-two class contributes one even zero, that the census is independent of the proof of Theorem 3, and that the domain is a rank-one/two determinant bound of 24. `docs/rank-one-and-two-census.md` states the same 26 positive and 44 diagonal rank-two classes, hence 50 positive classes after the 24 rank-one lattices, and one even zero on each diagonal class. The manuscript retains the initial norm4 `16` run and its 33 unresolved records, and it assigns the completing coefficients to a separate norm4 `25` run. Orders 4, 8, and 16 are supplied witnesses. The manuscript says they do not certify minimum order, unbounded minimum order, the glue criterion, or CM provenance. Those exclusions match the still-open rows of `PROOF-STATUS.md`.

The supplement’s Dn row states that all 700,070 even pairs in `D_2` through `D_10` are resolved and that the all-rank count is symbolic. `PROOF-STATUS.md` states the same split. This pass did not rerun those streams.

## 3. Supplement, manifest, and census bytes

Checked by reading the supplement, the verifier, the manifest bindings, and the four census reports, and by hashing the two census packs.

For E6, E7, E8, F4, G2, `D_2` through `D_10`, and the closed census, the supplement names the module, the certificate input, the saved report pair, and a command. The commands call `code/verify_even_pack.py`, `code/verify_e6_ambient.py`, `code/verify_exceptional.py`, `code/verify_lattice_census.py`, `code/verify_dn.py`, and `nonsimply_laced.verify_root_pack`. The Dn reports named by the pattern `reports/dn-dN-{ldl,ambient}.json` are the eighteen files `reports/dn-d2-ldl.json` through `reports/dn-d10-ambient.json` present in the manifest. The seven unittest modules named in the manuscript and the supplement contain 64 methods whose names begin with `test_`: 10 in `tests/test_exact_theta.py`, 9 in `tests/test_e6.py`, 11 in `tests/test_exceptional.py`, 12 in `tests/test_dn.py`, 6 in `tests/test_root_lattices.py`, 10 in `tests/test_lattice_census.py`, and 6 in `tests/test_rank_two.py`. This pass counted those definitions and did not run the suite. The sentence “All 64 tests pass” remains the statement in `PROOF-STATUS.md` and `README.md`.

The historical census is the row whose command is missing, as in the issue table. Its inputs and report pair are named.

`python scripts/verify_integrity.py

```
{"verdict": "file_bytes_verified", "files_checked": 114, "total_bytes": 10868976, "report_bindings_checked": 28, "mathematical_replay_performed": false}
```

and exited 0. The script checks the schema `work7-manuscript-artifacts-v1`, relative forward-slash paths, byte length, SHA-256, and each report binding’s `hash_field`. The manifest includes the current `Grok.md`, the manuscript, the supplement, both census packs, and the inventory verifier. It excludes itself.

An independent SHA-256 of the census bytes, separate from that script, gave:

| File | Bytes | SHA-256 |
|---|---:|---|
| `data/census/rank12-det24-closed.json` | 1023995 | `9b359c8736fe8687ad904a448382f6475dd955795ea015a3e4ca180d8df1e571` |
| `data/census/rank12-det24.json` | 947569 | `c8d121c7736791a86990800ddae198547689c1c617f27c7959b82cf3f8c9aa03` |

Both closed reports, box and LDL, record the first digest as `pack_sha256`. Both bound-16 reports record the second. Those four bindings are in the manifest, and the verifier’s exit 0 includes them. The saved summaries, read and not recomputed, are:

| Saved report | Verdict | Even pairs | Proved zero | Proved nonzero | Unresolved | Parity proved | Parity disproved | Parity unresolved |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| closed box and closed LDL | resolved | 772 | 44 | 728 | 0 | 50 | 44 | 0 |
| bound-16 box and bound-16 LDL | unresolved | 772 | 44 | 695 | 33 | 42 | 44 | 8 |

Both closed summaries record 94 retained representatives, with rank counts 24 and 70. The manuscript’s closed-pack sentence matches the 772, 728, 44, 50, and 44 fields. The manuscript’s 33 unresolved records match the bound-16 field. `docs/rank-one-and-two-census.md` explains the bound-16 lattice-level split 42, 44, and 8: the eight unresolved lattice verdicts are the rank-one cases `d = 17, …, 24`, and the rank-one theorem supplies the passage from 42 shell-positive lattices to the 50-class domain count. The manuscript attaches 50 and 44 to the closed pack. These are stored fields. They are not a new classification.

## 4. The three references

Checked by comparing `MANUSCRIPT.md` “References” with `docs/literature-and-source-scope.md` and with the formula numbers cited in `docs/theta-transformations.md`. This pass did not reopen the author PDF, the DLMF pages, or the Kane–Kim HTML.

The manuscript’s reference list has three items.

Chihara and Stanton are cited for equations (2.2) and (2.3) and for the remark after Theorem 3.4, with the author PDF identified as the pagination used in `docs/literature-and-source-scope.md`. The manuscript’s An section says the shell sum is a separate argument from `k_k(s,2,2k)`. That is the boundary recorded for the polynomial and its odd roots.

The DLMF citation lists (21.2.5), (21.2.6), (21.3.4), (21.3.6), (21.5.5), and (21.5.9). The manuscript body uses (21.2.5) and (21.3.6) for the series and for half-characteristic parity along `Ω = τ G`, and (21.5.5) and (21.3.4) for the unreduced basis change and the stabilizer sign. `docs/theta-transformations.md` cites (21.2.6) and (21.5.9) in that comparison. The manuscript describes the sign as a specialization of those laws and keeps the lattice classification in the lattice arguments. `docs/literature-and-source-scope.md` records the same six formulas as the ones read, and it keeps the Igusa and Mumford books off the list of checked book pages. The manuscript does not cite those books.

Kane and Kim, arXiv:2211.03987v2, are cited for §§2.1–2.2 and Proposition 2.3. The manuscript says the ternary spinor-genus theorem is unused. Proposition 5 uses §2.2 to mark the translation-by-one convention as inapplicable to `H` until cusp width is allowed to be arbitrary. That is the boundary between their modularity conventions and their ternary theorem.

## 5. Proposition 5

Checked by reading the proposition against `docs/coset-modularity.md` §§1–4. The earlier valuation gap and the unsigned cusp-frame gap are present in the note in repaired form: §2 concludes `g = 1` from the displayed valuation `0 = n v_ℓ(N) = v_ℓ(d) + v_ℓ(det H) ≥ n`, and §4 replaces `σ'` by `-σ'` when needed so that `d_c > 0` and `a_c > 0`. The manuscript’s factor `(d_c/8)^{-n/2}` is `(d_c/8)^{-k}` with `k = n/2`.

The common `z`-space in the manuscript is weight `n/2` on `Γ_0(16 N_{L_0})`, with `χ_{4 det L_0}` for odd rank and `χ_{(-1)^{n/2} 4 det L_0}` for even rank. A trivial doubled coset uses base `2L_0` and conductor one; a nontrivial coset uses base `L_0` and conductor two; the series identity is `z = τ/8`. The note’s §3 records the same space and records that `Γ_1(2)` is redundant because the diagonal entries are odd. The manuscript states the resulting group and omits the redundant intersection. The level of `2L_0` is the note’s theorem that this level equals `4N`; the manuscript states the group `Γ_0(16N_{L_0})` that this level produces.

Cusp width is stated in the proposition: width eight at infinity, arbitrary-width modular forms, and the source convention that requires translation by one. The odd-rank character is `χ_8(v) = (8/v)`. The note also treats the principal branch, negative lower-right entries, and `-I`. The manuscript proof defers those branch conventions in one sentence. That deferral is the optional editorial row. It leaves the stated transformation law in agreement with (3), (4), and (5).

The proposition claims holomorphy on this sufficient group. It claims neither cuspidality, nor a minimum level, nor a coefficient cutoff. The note’s truncated checks remain corroboration. This pass did not recompute them.

## 6. What this pass ran and read

Ran: `python scripts/verify_integrity.py

Read: `Grok.md`; `MANUSCRIPT.md`; `MANUSCRIPT-SUPPLEMENT.md`; `PROOF-STATUS.md`; `docs/coset-modularity.md`; the citation sections of `docs/literature-and-source-scope.md`; the `-I` and binary-witness lines of `docs/foundations-and-exact-enumeration.md`; the classification statement of `docs/root-lattice-classification.md`; the census verdict paragraphs of `docs/rank-one-and-two-census.md`; the README verifier block; the manifest path list and `report_bindings`; the F4 and G2 objects in `reports/root-lattices-audit.json`; and `verify_root_pack` in `code/nonsimply_laced.py`. The formula numbers (21.2.6) and (21.5.9) were confirmed as citations inside `docs/theta-transformations.md`.

Unrun: every exceptional, Dn, root-pack, and census replay, and `python -m unittest`. The saved reports were read as records. Their classifications were not reproduced by a new computation. The relative errors in `reports/coset-rescaling-checks.json` were not recomputed.

## 7. Remaining steps for the review draft

The reader-facing draft becomes complete as a review draft when the two reproducibility rows are edited in place:

1. Add the two print-only bound-16 commands to `MANUSCRIPT-SUPPLEMENT.md`.
2. Revise the README sentence in `MANUSCRIPT.md` so preservation-safe replay is the supplement block, and so the README `--output` and `code/lattice_census.py --output` lines are identified as commands that replace preserved census files.

The Proposition 5 branch paragraph and the F4/G2 hash bindings are optional exposition and inventory strengthenings. A later inventory regeneration belongs with an intentional snapshot, because the current manifest matches the current bytes, including this brief’s `Grok.md` as it stood when the verifier ran.

A fresh execution of the named 64-test command, and any fresh census or exceptional replay, would be new runs. This pass counted the tests and checked bytes. `PROOF-STATUS.md` remains the place that records the earlier passing run.

These edits finish the reader index of the review draft. The manuscript already places the rank-at-least-three classification, the symmetry converse in those ranks, minimum witness orders, the glue criterion, and CM provenance among the open questions, and it already states that the draft makes no priority claim. Completing the reader index leaves that scope and that claim as they stand. It is a completed review draft in the sense of this pass, and the manuscript’s own header remains the right description of its status.
