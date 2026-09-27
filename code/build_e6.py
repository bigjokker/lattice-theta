"""Construct one complete nonzero-shell certificate for each even E6 class."""

import argparse
from collections import Counter
import json
from pathlib import Path

from exact_theta import require, verify

ROOT = Path(__file__).resolve().parent.parent
E6_GRAM = [
    [2, -1, 0, 0, 0, 0],
    [-1, 2, -1, 0, 0, 0],
    [0, -1, 2, -1, 0, -1],
    [0, 0, -1, 2, -1, 0],
    [0, 0, 0, -1, 2, 0],
    [0, 0, -1, 0, 0, 2],
]


def bits(mask):
    return [(mask >> i) & 1 for i in range(6)]


def build(norm4_bound=4):
    entries, histogram, unresolved = [], Counter(), []
    for am in range(64):
        a = bits(am)
        for bm in range(64):
            if (am & bm).bit_count() % 2:
                continue
            b = bits(bm)
            document = {"schema": "work7-theta-v1", "G": E6_GRAM, "a": a, "b": b,
                        "evidence": {"kind": "search", "norm4_bound": norm4_bound}}
            result = verify(document)
            if result["verdict"] != "does_not_vanish":
                unresolved.append([am, bm])
                continue
            N, c = result["nonzero_norm4"], result["signed_coefficient"]
            entries.append({"a": a, "b": b, "evidence": {"kind": "shell", "norm4": N,
                                                              "coefficient": c}})
            histogram[N] += 1
    print(json.dumps({"resolved": len(entries), "witness_norm4_counts": dict(sorted(histogram.items())),
                      "unresolved": unresolved}), flush=True)
    require(not unresolved, "E6 builder left unresolved classes; no complete pack written")
    return {"schema": "work7-even-theta-pack-v1", "name": "E6 parity converse",
            "G": E6_GRAM, "representatives": "binary-primal-and-dual",
            "provenance": "simple roots ordered along a five-node chain, then its central branch; see E6-MILESTONE.md",
            "certificates": entries}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--norm4-bound", type=int, default=4)
    parser.add_argument("--output", type=Path, default=ROOT / "data" / "e6-even.json")
    args = parser.parse_args()
    pack = build(args.norm4_bound)
    # One compact entry per line makes the complete finite table easy to inspect.
    header = {k: v for k, v in pack.items() if k != "certificates"}
    text = json.dumps(header, indent=2)[:-2] + ',\n  "certificates": [\n'
    text += ",\n".join("    " + json.dumps(entry, separators=(",", ":"))
                       for entry in pack["certificates"])
    text += "\n  ]\n}\n"
    # Check the serialization before writing the artifact.
    require(json.loads(text) == pack, "pack serialization mismatch")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
