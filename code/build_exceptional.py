"""Build complete E7/E8 even-characteristic classifications with exact evidence."""

import argparse
from collections import Counter
import json
from pathlib import Path
from time import perf_counter

from exact_theta import verify
from exceptional_batch import (export_certificate, ldl_histograms, root_gram,
                                symmetry_forest, validate_pack, walsh)

ROOT = Path(__file__).resolve().parent.parent


def build(rank, bound=8):
    G = root_gram(rank)
    start = perf_counter()
    histograms, stats = ldl_histograms(G, bound)
    coefficients = {key: walsh(histogram, rank) for key, histogram in histograms.items()}
    shell_seconds = perf_counter() - start
    entries, candidates, counts = [], [], Counter()
    for am in range(1 << rank):
        norms = sorted(N for aa, N in coefficients if aa == am)
        for bm in range(1 << rank):
            if (am & bm).bit_count() % 2:
                continue
            N = next((N for N in norms if coefficients[(am, N)][bm]), None)
            if N is None:
                candidates.append((am, bm))
                evidence = None
            else:
                evidence = {"kind": "shell", "norm4": N, "coefficient": coefficients[(am, N)][bm]}
                counts[N] += 1
            entries.append({"a_mask": am, "b_mask": bm, "evidence": evidence})
    print(json.dumps({"stage": "shells", "rank": rank, **stats, "seconds": round(shell_seconds, 3),
                      "nonzero": sum(counts.values()), "candidates": len(candidates)}), flush=True)
    start = perf_counter()
    generators, seeds, nodes, indices = symmetry_forest(G, candidates)
    for entry in entries:
        if entry["evidence"] is None:
            entry["evidence"] = {"kind": "symmetry", "node": indices[(entry["a_mask"], entry["b_mask"])]}
    pack = {"schema": "work7-theta-classification-v1", "name": f"E{rank} even-characteristic classification",
            "G": G, "representatives": "binary-primal-and-dual-masks", "norm4_bound": bound,
            "provenance": "simple roots ordered Bourbaki (1,3,4,...,n,2); see E7-E8-MILESTONE.md",
            "generators": generators, "seeds": seeds, "zero_nodes": nodes, "classes": entries}
    validate_pack(pack)
    print(json.dumps({"stage": "symmetry", "rank": rank, "zero_certificates": len(nodes),
                      "seed_witnesses": len(seeds), "seconds": round(perf_counter() - start, 3)}), flush=True)
    return pack


def write_pack(pack, path):
    header = {k: v for k, v in pack.items() if k not in ("zero_nodes", "classes")}
    text = json.dumps(header, indent=2)[:-2]
    for name in ("zero_nodes", "classes"):
        text += f',\n  "{name}": [\n'
        text += ",\n".join("    " + json.dumps(row, separators=(",", ":")) for row in pack[name])
        text += "\n  ]"
    text += "\n}\n"
    if json.loads(text) != pack:
        raise ValueError("serialization mismatch")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ranks", type=int, nargs="+", choices=(7, 8))
    parser.add_argument("--norm4-bound", type=int, default=8)
    args = parser.parse_args()
    for rank in args.ranks:
        pack = build(rank, args.norm4_bound)
        write_pack(pack, ROOT / "data" / f"e{rank}-classification.json")
        zero = next(e for e in pack["classes"] if e["evidence"]["kind"] == "symmetry")
        nonzero = max((e for e in pack["classes"] if e["evidence"]["kind"] == "shell"),
                      key=lambda e: e["evidence"]["norm4"])
        for label, entry in (("zero", zero), ("nonzero", nonzero)):
            doc = export_certificate(pack, entry)
            verify(doc)
            (ROOT / "data" / f"e{rank}-{label}.json").write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
