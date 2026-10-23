"""Freeze the audited Theorem S package (run once; refuses to overwrite).

Copies, unchanged, every file that the audited proof rests on into this folder, and writes MANIFEST.json with
the sha256 of each copy and of its source.  Source modules are found by importing the three audited modules and
collecting every project module they load (transitive dependencies), so the snapshot can run on its own.
"""
import hashlib
import json
from pathlib import Path
import shutil
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))

AUDITED = {   # hashes recorded in cells-v3.json and verified in RANK4-THEOREM-S-FINAL-AUDIT-ASTRA.md
    "rank4_six_cells.py": "328a53ee5b36b66a36569f1f46a887d8ade11dbb15b639a06ad02d135a4b5b20",
    "rank4_secondary_fan.py": "fc206a7cbd915951c1c6372e29eb1b842d4c45e6b9f0ed379bd2668f1daa8cdc",
    "rank4_six_minima.py": "a133954d592bd0832818bad4482b13b2139a9caec089c6767d77b2b14d1c55cb",
    "results/rank4-secondary-fan/fan-exact.json": "215e90f39d37edf40977a2cdf88961b855bbda1bef14b7fdf28fa209e1f56135",
}

FILES = {
    "data": ["results/rank4-secondary-fan/fan-exact.json", "results/rank4-secondary-fan/fan.json",
             "results/rank4-six-cells/cells.json", "results/rank4-six-cells/cells-v2.json",
             "results/rank4-six-cells/cells-v3.json", "results/rank4-six-cells/cells-astra-final-20261023.json"],
    "independent_replay_grokbot": ["bot_work/bot/coverage/replay.py", "bot_work/bot/coverage/ownT.py",
                                   "bot_work/bot/coverage/dispositions.jsonl.gz", "bot_work/bot/coverage/summary.json",
                                   "bot_work/bot/coverage/ownT.json",
                                   "bot_work/bot/coverage/COVERAGE_REPLAY_GrokBot_20261023.md"],
    "literature_audit_grok": ["bot_work/grok/proof_audit/PROOF_AUDIT_Grok_20261023.md",
                              "bot_work/grok/proof_audit/check_proof_audit.py",
                              "bot_work/grok/proof_audit/check_proof_audit.json",
                              "bot_work/grok/proof_audit/match_primitive_rays.py",
                              "bot_work/grok/proof_audit/primitive_match.json",
                              "bot_work/grok/proof_audit/adjacent_cones.py"],
    "final_audit_astra": ["RANK4-THEOREM-S-FINAL-AUDIT-ASTRA.md",
                          "results/rank4-final-audit-20261023/check_bundle.py",
                          "results/rank4-final-audit-20261023/bundle-audit.json",
                          "results/rank4-final-audit-20261023/check_published.py",
                          "results/rank4-final-audit-20261023/published-cones-and-hashes.json",
                          "results/rank4-final-audit-20261023/manifest.json"],
    "checks_claude": ["results/rank4-secondary-fan/match_published.py", "results/rank4-secondary-fan/match-published.json",
                      "results/rank4-secondary-fan/regenerate_v1.py", "results/rank4-secondary-fan/regenerate_v1.log",
                      "results/rank4-secondary-fan/fan-exact.log", "results/rank4-six-cells/validate.py",
                      "results/rank4-six-cells/validate-recheck.json", "results/rank4-six-cells/validate-controls.json",
                      "results/rank4-six-cells/validate-fan2.json",
                      "results/rank4-six-cells/grokbot-replay-crosscheck/summary-on-regenerated-fan.json",
                      "results/rank4-six-cells/grokbot-replay-crosscheck/summary-on-fan-exact.json",
                      "test_rank4_six_cells.py", "test_rank4_six_minima.py"],
    "documents": ["RANK4-SIX-CELLS-CLASSIFICATION.md", "RANK4-SIX-CELLS-PROOF-STATUS.md",
                  "RANK4-SIX-CELLS-QUICK-REVIEW-ASTRA.md", "RANK4-SIX-CELLS-QUICK-REVIEW-RESOLUTION.md",
                  "RANK4-SIX-MINIMA-COMMON-LATTICE.md", "RANK4-MINUS-TYPE.md", "RANK4-STAGE2-STRATEGY.md",
                  "RANK4-PARITY-CONJECTURE.md", "Astra-handoff-20261023-morning.md"],
}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def project_modules():
    before = set(sys.modules)
    import rank4_six_cells, rank4_secondary_fan, rank4_six_minima    # noqa: F401,E401
    import verify_lattice_census, verify_rank4_feasibility, verify_rank4_parity    # noqa: F401,E401
    out = []
    for name, mod in list(sys.modules.items()):
        f = getattr(mod, "__file__", None)
        if f and Path(f).resolve().parent == ROOT and f.endswith(".py"):
            out.append(Path(f).resolve().name)
    return sorted(set(out))


def main():
    manifest_path = HERE / "MANIFEST.json"
    assert not manifest_path.exists(), "package already built: refusing to overwrite"
    for rel, h in AUDITED.items():
        assert sha(ROOT / rel) == h, f"{rel} differs from the audited hash: do not package"
    FILES["sources"] = project_modules()
    entries = []
    for group, rels in FILES.items():
        for rel in rels:
            src = ROOT / rel
            assert src.exists(), rel
            dst = HERE / group / Path(rel).name
            dst.parent.mkdir(parents=True, exist_ok=True)
            assert not dst.exists(), f"name clash in {group}: {dst.name}"
            shutil.copy2(src, dst)
            h = sha(dst)
            assert h == sha(src)
            entries.append({"group": group, "file": f"{group}/{dst.name}", "source": rel.replace(chr(92), "/"),
                            "sha256": h, "bytes": dst.stat().st_size, "audited_hash": AUDITED.get(rel)})
    doc = {"package": "Theorem S and the minus-type corollary (computer-assisted proof), audited 2026-10-23",
           "audit": "RANK4-THEOREM-S-FINAL-AUDIT-ASTRA.md (accepted)", "audited_hashes": AUDITED,
           "files": entries, "count": len(entries)}
    manifest_path.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({"files": len(entries), "groups": {g: sum(e["group"] == g for e in entries) for g in FILES},
                      "sources": FILES["sources"]}))


if __name__ == "__main__":
    main()
