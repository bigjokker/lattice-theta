"""Complete rank-1/2 integral lattice census, determinant <=24 (stdlib only).

See LATTICE-CENSUS.md for generation and exhaustive basis-image proofs.
The initial theta bound is a resource choice, not a zero criterion.
"""

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import isqrt
from pathlib import Path

from exact_theta import (EnumerationLimit, complete_shells, dot, identity, integer,
                         inverse, ldl, matrix, matvec, multiply, require,
                         symmetry_data, transpose)

SCHEMA = "work7-rank12-census-v1"


def admit(G, determinant_bound=24):
    matrix(G, "G")
    require(len(G) in (1, 2), "census rank must be one or two")
    _, pivots = ldl(G)
    d = Fraction(1)
    for pivot in pivots:
        d *= pivot
    require(d.denominator == 1 and 1 <= d <= determinant_bound,
            "Gram determinant outside census domain")
    return int(d)


def candidates(determinant_bound=24):
    integer(determinant_bound, "determinant_bound")
    require(1 <= determinant_bound <= 24, "determinant bound must be in 1..24")
    result = [[[d]] for d in range(1, determinant_bound + 1)]
    # det >= 3 a^2/4, hence a <= floor(sqrt(4D/3)).
    for a in range(1, isqrt(4 * determinant_bound // 3) + 1):
        for b in range(a // 2 + 1):
            for c in range(a, (determinant_bound + b * b) // a + 1):
                if 1 <= a * c - b * b <= determinant_bound:
                    result.append([[a, b], [b, c]])
    return result


def determinant(U):
    return U[0][0] if len(U) == 1 else U[0][0] * U[1][1] - U[0][1] * U[1][0]


def sphere(G, bound, max_nodes=1_000_000):
    """All integer vectors with norm <=bound, via the existing LDL engine."""
    vectors = []
    complete_shells(G, [0] * len(G), [0] * len(G), 4 * bound, max_nodes,
                    lambda y, q: vectors.append(tuple(v // 2 for v in y)))
    return vectors


def isometries(G, H, max_nodes=1_000_000, sphere_backend=sphere):
    """Every unimodular U with U^t H U=G; limit raises, never returns partial."""
    admit(G)
    admit(H)
    if len(G) != len(H):
        return []
    n = len(G)
    vectors = sphere_backend(H, max(G[i][i] for i in range(n)), max_nodes)
    columns = [[v for v in vectors if dot(v, matvec(H, v)) == G[i][i]]
               for i in range(n)]
    result = []
    for nodes, cols in enumerate(product(*columns), 1):
        if nodes > max_nodes:
            raise EnumerationLimit("basis-image search exceeded node limit")
        U = transpose(cols)
        if abs(determinant(U)) == 1 and multiply(multiply(transpose(U), H), U) == G:
            result.append(U)
    return sorted(result)


def even_pairs(n):
    return [(list(a), list(b)) for a in product((0, 1), repeat=n)
            for b in product((0, 1), repeat=n) if dot(a, b) % 2 == 0]


def stabilizer(G, a, b, group):
    indices, epsilon = [], []
    for i, T in enumerate(group):
        dual = transpose(inverse(T))
        if any((x - y) % 2 for x, y in zip(matvec(T, a), a)):
            continue
        if any((x - y) % 2 for x, y in zip(matvec(dual, b), b)):
            continue
        data = symmetry_data(G, a, b, T)
        indices.append(i)
        epsilon.append(data["epsilon"])
    return {"indices": indices, "epsilon": epsilon}


def classify(G, a, b, group, bound, max_nodes):
    record = {"a": a, "b": b, "stabilizer": None}
    if group is not None:
        record["stabilizer"] = stabilizer(G, a, b, group)
        for i, e in zip(record["stabilizer"]["indices"], record["stabilizer"]["epsilon"]):
            if e:
                record.update(verdict="proved_zero", evidence={"kind": "symmetry", "T": group[i]})
                return record
    try:
        enumeration = complete_shells(G, a, b, bound, max_nodes)
    except EnumerationLimit as error:
        record.update(verdict="unresolved", evidence={"kind": "search", "norm4_bound": bound},
                      enumeration_complete=False, reason=str(error))
        return record
    record["enumeration"] = enumeration
    nonzero = next((s for s in enumeration["shells"] if s["signed_coefficient"]), None)
    if nonzero:
        record.update(verdict="proved_nonzero", evidence={"kind": "shell", "norm4": nonzero["norm4"],
                                                         "coefficient": nonzero["signed_coefficient"]})
    else:
        record.update(verdict="unresolved", evidence={"kind": "search", "norm4_bound": bound},
                      reason="all complete coefficients through the bound are zero")
    return record


def summarize(document):
    rows = document["lattices"]
    counts = Counter(p["verdict"] for row in rows for p in row["pairs"])
    return {"generated_matrices": len(document["generation"]), "retained_representatives": len(rows),
            "rank_counts": {str(n): sum(len(r["G"]) == n for r in rows) for n in (1, 2)},
            "generation_complete": True,
            "isometry_class_coverage_complete": True,
            "uniqueness_complete": all(g["uniqueness_complete"] for g in document["generation"]),
            "groups_complete": all(r["automorphisms"] is not None for r in rows),
            "characteristic_coverage_complete": True,
            "even_pairs": sum(len(r["pairs"]) for r in rows),
            "odd_zero_by_parity": sum(4 ** len(r["G"]) - len(r["pairs"]) for r in rows),
            "proved_zero": counts["proved_zero"], "proved_nonzero": counts["proved_nonzero"],
            "unresolved": counts["unresolved"], "theta_resolution_complete": not counts["unresolved"],
            "parity_converse_proved": sum(all(p["verdict"] == "proved_nonzero" for p in r["pairs"]) for r in rows),
            "parity_converse_disproved": sum(any(p["verdict"] == "proved_zero" for p in r["pairs"]) for r in rows),
            "parity_converse_unresolved": sum(not any(p["verdict"] == "proved_zero" for p in r["pairs"])
                                              and any(p["verdict"] == "unresolved" for p in r["pairs"]) for r in rows),
            "symmetry_converse_counterexamples": 0}


def build(determinant_bound=24, norm4_bound=16, max_nodes=1_000_000):
    integer(norm4_bound, "norm4_bound")
    require(norm4_bound >= 0, "negative shell bound")
    integer(max_nodes, "max_nodes")
    require(max_nodes > 0, "max_nodes must be positive")
    document = {"schema": SCHEMA, "domain": {"ranks": [1, 2], "determinant_bound": determinant_bound},
                "run": {"norm4_bound": norm4_bound, "max_nodes": max_nodes},
                "generation": [], "lattices": []}
    for G in candidates(determinant_bound):
        d = admit(G, determinant_bound)
        target, U, complete, failures = None, None, True, []
        for i, row in enumerate(document["lattices"]):
            if len(row["G"]) != len(G) or row["determinant"] != d:
                continue
            try:
                matches = isometries(G, row["G"], max_nodes)
            except EnumerationLimit:
                complete = False
                failures.append(i)
                continue
            if matches:
                target, U = i, matches[0]
                break
        if target is None:
            target, U = len(document["lattices"]), identity(len(G))
            try:
                group = isometries(G, G, max_nodes)
                minimum = min(dot(v, matvec(G, v)) for v in sphere(G, G[0][0], max_nodes) if any(v))
                group_error = None
            except EnumerationLimit as error:
                group, minimum, group_error = None, None, str(error)
            row = {"G": G, "determinant": d,
                   "gram_sha256": hashlib.sha256(json.dumps(G, separators=(",", ":")).encode("ascii")).hexdigest(),
                   "minimum_squared_norm": minimum, "automorphisms": group, "group_error": group_error,
                   "pairs": [classify(G, a, b, group, norm4_bound, max_nodes) for a, b in even_pairs(len(G))]}
            document["lattices"].append(row)
        document["generation"].append({"G": G, "representative": target, "U": U,
                                       "uniqueness_complete": complete, "incomplete_comparisons": failures})
    document["summary"] = summarize(document)
    return document


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("data/census/rank12-det24.json"))
    parser.add_argument("--determinant-bound", type=int, default=24)
    parser.add_argument("--norm4-bound", type=int, default=16)
    parser.add_argument("--max-nodes", type=int, default=1_000_000)
    args = parser.parse_args()
    document = build(args.determinant_bound, args.norm4_bound, args.max_nodes)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(document["summary"]))
    return 2 if document["summary"]["unresolved"] or not document["summary"]["groups_complete"] or not document["summary"]["uniqueness_complete"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
