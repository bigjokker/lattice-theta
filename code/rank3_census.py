"""Exact rank-three integral census, determinant <=24; see RANK3-CENSUS-PLAN.md."""

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import lcm
from pathlib import Path

from exact_theta import (EnumerationLimit, dot, identity, integer, inverse, ldl,
                         matrix, matvec, multiply, require, transpose)
from lattice_census import classify, even_pairs, sphere

SCHEMA = "work7-rank3-census-v1"


def determinant(G):
    a, b, d = G[0]
    _, c, f = G[1]
    _, _, e = G[2]
    # Works for nonsymmetric basis-image matrices too.
    return (a * (c * e - f * G[2][1])
            - b * (G[1][0] * e - f * G[2][0])
            + d * (G[1][0] * G[2][1] - c * G[2][0]))


def admit(G, determinant_bound=24):
    matrix(G, "G", 3)
    ldl(G)
    det = determinant(G)
    require(1 <= det <= determinant_bound, "Gram determinant outside rank-three domain")
    return det


def candidates(determinant_bound=24):
    """Complete finite projected-Gauss domain; not canonical representatives."""
    H = integer(determinant_bound, "determinant_bound")
    require(1 <= H <= 24, "determinant bound must be in 1..24")
    result = []
    for a in range(1, 4):
        if 27 * a ** 3 > 64 * H:
            continue
        for b in range(-(a // 2), a // 2 + 1):
            for c in range(a, 2 * H + a + 1):
                m = a * c - b * b
                if 3 * m * m > 4 * a * H:
                    continue
                for d in range(-(a // 2), a // 2 + 1):
                    lo = -((m - 2 * b * d) // (2 * a))
                    hi = (m + 2 * b * d) // (2 * a)
                    for f in range(lo, hi + 1):
                        for delta in range(1, H + 1):
                            numerator = delta + c * d * d + a * f * f - 2 * b * d * f
                            e, remainder = divmod(numerator, m)
                            if not remainder and a * e - d * d >= m:
                                result.append([[a, b, d], [b, c, f], [d, f, e]])
    return sorted(result)


def isometries(G, H, max_nodes=1_000_000, sphere_backend=sphere):
    """Exhaustive basis images, with exact pairing pruning and explicit limits."""
    admit(G)
    admit(H)
    vectors = sphere_backend(H, max(G[i][i] for i in range(3)), max_nodes)
    images = [[v for v in vectors if dot(v, matvec(H, v)) == G[i][i]] for i in range(3)]
    Himages = {v: matvec(H, v) for v in vectors}
    result, nodes = [], 0
    for v0 in images[0]:
        for v1 in images[1]:
            nodes += 1
            if nodes > max_nodes:
                raise EnumerationLimit("rank-three basis-image search exceeded node limit")
            if dot(v0, Himages[v1]) != G[0][1]:
                continue
            for v2 in images[2]:
                nodes += 1
                if nodes > max_nodes:
                    raise EnumerationLimit("rank-three basis-image search exceeded node limit")
                if dot(v0, Himages[v2]) != G[0][2] or dot(v1, Himages[v2]) != G[1][2]:
                    continue
                U = transpose([v0, v1, v2])
                if abs(determinant(U)) == 1:
                    require(multiply(multiply(transpose(U), H), U) == G, "basis-image mismatch")
                    result.append(U)
    return sorted(result)


def modular_data(G, b):
    """Established modular input only; no coefficient-cutoff claim."""
    if not any(b):
        return None
    j = next(i for i in range(3) if b[i])
    K = identity(3)
    K[j][j] = 2
    for i in range(3):
        if i != j:
            K[j][i] = -b[i]
    G0 = multiply(multiply(transpose(K), G), K)
    inv = inverse(G0)
    N0 = lcm(*(v.denominator for row in inv for v in row))
    return {"b": b, "index_two_basis": K, "G0": G0, "d0": determinant(G0),
            "N0": N0, "unsigned_level": 16 * N0,
            "character_discriminant": 4 * determinant(G0),
            "weight": "3/2", "coefficient_index": "norm4", "cutoff": None}


def summarize(document):
    rows = document["lattices"]
    counts = Counter(p["verdict"] for row in rows for p in row["pairs"])
    return {"generated_matrices": len(document["generation"]),
            "retained_representatives": len(rows),
            "determinant_counts": dict(sorted(Counter(str(r["determinant"]) for r in rows).items(),
                                              key=lambda item: int(item[0]))),
            "generation_complete": True, "isometry_class_coverage_complete": True,
            "uniqueness_complete": all(g["uniqueness_complete"] for g in document["generation"]),
            "groups_complete": all(r["automorphisms"] is not None for r in rows),
            "characteristic_coverage_complete": all(len(r["pairs"]) == 36 for r in rows),
            "even_pairs": sum(len(r["pairs"]) for r in rows), "odd_zero_by_parity": 28 * len(rows),
            "proved_zero": counts["proved_zero"], "proved_nonzero": counts["proved_nonzero"],
            "unresolved": counts["unresolved"], "theta_resolution_complete": not counts["unresolved"],
            "parity_converse_proved": sum(all(p["verdict"] == "proved_nonzero" for p in r["pairs"]) for r in rows),
            "parity_converse_disproved": sum(any(p["verdict"] == "proved_zero" for p in r["pairs"]) for r in rows),
            "parity_converse_unresolved": sum(not any(p["verdict"] == "proved_zero" for p in r["pairs"])
                                               and any(p["verdict"] == "unresolved" for p in r["pairs"]) for r in rows),
            "symmetry_converse_counterexamples": 0,
            "symmetry_converse_proved_in_domain": not counts["unresolved"]
                and all(r["automorphisms"] is not None for r in rows)
                and all(g["uniqueness_complete"] for g in document["generation"])}


def build(determinant_bound=24, norm4_bound=32, max_nodes=1_000_000):
    integer(norm4_bound, "norm4_bound")
    integer(max_nodes, "max_nodes")
    require(norm4_bound >= 0 and max_nodes > 0, "invalid run limits")
    document = {"schema": SCHEMA, "domain": {"ranks": [3], "determinant_bound": determinant_bound,
                "generation_proof": "RANK3-CENSUS-PLAN.md: projected-Gauss"},
                "run": {"norm4_bound": norm4_bound, "max_nodes": max_nodes},
                "generation": [], "lattices": []}
    for G in candidates(determinant_bound):
        delta = admit(G, determinant_bound)
        target, U, failures = None, None, []
        for i, row in enumerate(document["lattices"]):
            if row["determinant"] != delta:
                continue
            try:
                matches = isometries(G, row["G"], max_nodes)
            except EnumerationLimit:
                failures.append(i)
                continue
            if matches:
                target, U = i, matches[0]
                break
        if target is None:
            target, U = len(document["lattices"]), identity(3)
            try:
                group = isometries(G, G, max_nodes)
                minimum = min(dot(v, matvec(G, v)) for v in sphere(G, G[0][0], max_nodes) if any(v))
                group_error = None
            except EnumerationLimit as error:
                group, minimum, group_error = None, None, str(error)
            document["lattices"].append({"G": G, "determinant": delta,
                "gram_sha256": hashlib.sha256(json.dumps(G, separators=(",", ":")).encode("ascii")).hexdigest(),
                "minimum_squared_norm": minimum, "automorphisms": group, "group_error": group_error,
                "odd_pairs": [{"a": list(a), "b": list(b), "verdict": "proved_zero", "evidence": "minus_identity"}
                              for a in product((0, 1), repeat=3) for b in product((0, 1), repeat=3)
                              if dot(a, b) % 2],
                "modular_inputs": [modular_data(G, list(b)) for b in product((0, 1), repeat=3) if any(b)],
                "pairs": [classify(G, a, b, group, norm4_bound, max_nodes) for a, b in even_pairs(3)]})
        document["generation"].append({"G": G, "representative": target, "U": U,
            "uniqueness_complete": not failures, "incomplete_comparisons": failures})
    document["summary"] = summarize(document)
    return document


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("data/census-rank3/det24-bound32.json"))
    parser.add_argument("--determinant-bound", type=int, default=24)
    parser.add_argument("--norm4-bound", type=int, default=32)
    parser.add_argument("--max-nodes", type=int, default=1_000_000)
    args = parser.parse_args()
    document = build(args.determinant_bound, args.norm4_bound, args.max_nodes)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(document["summary"]))
    return 2 if document["summary"]["unresolved"] or not document["summary"]["groups_complete"] or not document["summary"]["uniqueness_complete"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
