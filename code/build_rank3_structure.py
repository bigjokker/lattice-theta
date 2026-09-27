"""Certify the existing 120-class domain using exact structural proofs."""

import argparse
from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path

from exact_theta import (EnumerationLimit, dot, inverse, matvec, multiply, require,
                         integer, symmetry_data, transpose)
from rank3_structure import decomposition, graph_type, obtuse_superbase, shell_claim

SCHEMA = "work7-rank3-structure-v1"


def factor_pair(G, a, b, split):
    B = split["basis"]
    ap = [int(x) % 2 for x in matvec(inverse(B), a)]
    bp = [int(x) % 2 for x in matvec(transpose(B), b)]
    blocks = [[0], [1], [2]] if split["type"] == "three_lines" else [[0], [1, 2]]
    odd = next((block for block in blocks if sum(ap[i]*bp[i] for i in block) % 2), None)
    if odd is None:
        return {"a": a, "b": b, "verdict": "proved_nonzero", "proof": "rank12-product",
                "factor_a": ap, "factor_b": bp}
    D = [[(-1 if i in odd else 1) * int(i == j) for j in range(3)] for i in range(3)]
    T = multiply(multiply(B, D), inverse(B))
    require(all(x.denominator == 1 for row in T for x in row), "nonintegral factor witness")
    T = [[int(x) for x in row] for row in T]
    require(symmetry_data(G, a, b, T)["epsilon"] == 1, "factor witness has zero sign")
    return {"a": a, "b": b, "verdict": "proved_zero", "proof": "odd-factor", "T": T,
            "factor_a": ap, "factor_b": bp}


def build(census, census_sha256, max_nodes=1_000_000):
    integer(max_nodes, "max_nodes")
    require(max_nodes > 0, "max_nodes must be positive")
    require(census["schema"] == "work7-rank3-census-v1" and census["summary"]["symmetry_converse_proved_in_domain"],
            "resolved rank-three census required")
    rows, limits = [], []
    for i, original in enumerate(census["lattices"]):
        G = original["G"]
        row = {"lattice": i, "G": G, "decomposition": None, "superbase": None, "pairs": []}
        rows.append(row)
        try:
            row["decomposition"] = decomposition(G, max_nodes)
            if row["decomposition"]["verdict"] == "indecomposable":
                row["superbase"] = obtuse_superbase(G, max_nodes)
                require(graph_type(row["superbase"]["conorms"]) != "articulation", "indecomposable graph articulates")
            for a in product((0, 1), repeat=3):
                for b in product((0, 1), repeat=3):
                    a, b = list(a), list(b)
                    if dot(a, b) % 2:
                        record = {"a": a, "b": b, "verdict": "proved_zero", "proof": "odd-parity",
                                  "T": [[-int(i == j) for j in range(3)] for i in range(3)]}
                    elif row["decomposition"]["verdict"] == "decomposable":
                        record = factor_pair(G, a, b, row["decomposition"]["split"])
                    else:
                        record = {"a": a, "b": b, "verdict": "proved_nonzero", "proof": "obtuse-graph-shell",
                                  "shell": shell_claim(G, a, b, row["superbase"])}
                    if dot(a, b) % 2 == 0:
                        previous = next(p for p in original["pairs"] if p["a"] == a and p["b"] == b)
                        require(record["verdict"] == previous["verdict"], "structure contradicts saved census")
                    row["pairs"].append(record)
        except EnumerationLimit as error:
            row["pairs"] = []
            limits.append({"lattice": i, "reason": str(error)})
    types = Counter(r["decomposition"]["split"]["type"] if r["decomposition"] and
                    r["decomposition"]["verdict"] == "decomposable" else
                    r["decomposition"]["verdict"] if r["decomposition"] else "unresolved" for r in rows)
    counts = Counter(p["verdict"] for r in rows for p in r["pairs"])
    summary = {"classes": len(rows), "structure_types": dict(types), "even_pairs": sum(
        dot(p["a"], p["b"]) % 2 == 0 for r in rows for p in r["pairs"]),
        "all_pairs": sum(len(r["pairs"]) for r in rows), "proved_zero": counts["proved_zero"],
        "proved_nonzero": counts["proved_nonzero"], "complete": not limits,
        "positive_graphs": dict(Counter(graph_type(r["superbase"]["conorms"]) for r in rows if r["superbase"])),
        "symbolic_shell_cases": dict(Counter(p["shell"]["case"] for r in rows for p in r["pairs"] if "shell" in p))}
    return {"schema": SCHEMA, "census_sha256": census_sha256,
            "domain": {"ranks": [3], "determinant_bound": 24}, "max_nodes": max_nodes,
            "proof_note": "RANK3-STRUCTURE.md", "rows": rows, "limits": limits, "summary": summary}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--census", type=Path, default=Path("data/census-rank3/det24-bound32.json"))
    parser.add_argument("--output", type=Path, default=Path("data/census-rank3/structure-det24.json"))
    parser.add_argument("--max-nodes", type=int, default=1_000_000)
    args = parser.parse_args()
    raw = args.census.read_bytes()
    document = build(json.loads(raw), hashlib.sha256(raw).hexdigest(), args.max_nodes)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(document["summary"]))
    return 0 if document["summary"]["complete"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
