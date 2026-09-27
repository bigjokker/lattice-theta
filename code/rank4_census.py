"""Exact bounded rank-four census; incomplete searches never certify a census."""
import argparse
from collections import Counter
import json
from pathlib import Path
import time

from exact_theta import dot, identity, require
from rank4_feasibility import (BITS, action_tables, admit, candidates, characteristic,
                               cutoff_build, isometries, modular_data)


def build(H=4, bound=32, max_nodes=1_000_000):
    require(type(bound) is int and bound >= 0, "invalid shell bound")
    start = time.perf_counter()
    forms = candidates(H, max_nodes)
    rows, generation = [], []
    for G in forms:
        delta = admit(G)
        match = None
        for i, row in enumerate(rows):
            if row["determinant"] == delta:
                maps = isometries(G, row["G"], max_nodes)
                if maps:
                    match = (i, maps[0])
                    break
        if match is None:
            match = (len(rows), identity(4))
            rows.append({"G": G, "determinant": delta})
        generation.append({"G": G, "representative": match[0], "T": match[1]})
    # Any interrupted generation, comparison or group search aborts the build.
    # No partially generated or partially deduplicated pack is written.
    for row in rows:
        G = row["G"]
        group = isometries(G, G, max_nodes)
        tables = action_tables(group)
        pairs = []
        for a in BITS:
            for b in BITS:
                p = characteristic(G, a, b, group, tables, bound, max_nodes)
                if p["verdict"] == "unresolved" and dot(a, b) % 2 == 0 and any(b):
                    cert = cutoff_build(G, a, b, max_nodes=max_nodes)
                    p["cutoff_certificate"] = cert
                    if cert["evidence"]["verdict"] == "proved_zero_by_modular_cutoff":
                        p["verdict"] = "proved_zero_without_symmetry"
                    elif cert["evidence"]["verdict"] == "proved_nonzero":
                        p["verdict"] = "proved_nonzero_by_modular_replay"
                pairs.append(p)
        row.update(automorphisms=group, pairs=pairs,
                   modular_kernels=[{"b": b, **modular_data(G, [0]*4, b)} for b in BITS if any(b)])
    counts = Counter(p["verdict"] for r in rows for p in r["pairs"] if dot(p["a"], p["b"]) % 2 == 0)
    return {"schema": "work7-rank4-census-v1", "domain": {"rank": 4, "determinant_bound": H},
            "run": {"norm4_bound": bound, "max_nodes": max_nodes},
            "generation": generation, "lattices": rows,
            "summary": {"candidate_bases": len(forms), "lattice_classes": len(rows),
                        "even_pairs": 136*len(rows), "odd_pairs": 120*len(rows),
                        "even_verdicts": dict(counts)},
            "seconds": time.perf_counter()-start}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--det", type=int, default=4)
    parser.add_argument("--bound", type=int, default=32)
    parser.add_argument("--max-nodes", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), "choose a fresh output path")
    doc = build(args.det, args.bound, args.max_nodes)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(doc, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"summary": doc["summary"], "seconds": doc["seconds"]}))


if __name__ == "__main__":
    main()
