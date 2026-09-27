"""Replay census coverage, full groups and theta evidence using exact boxes.

The box proof is |y_i|^2 <= bound*(G^-1)_ii. No LDL enumeration
or census classification routine is used by the default replay backend.
"""

import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import isqrt
from pathlib import Path

from exact_theta import (EnumerationLimit, InvalidCertificate, complete_shells,
                         dot, identity, integer, inverse, matrix, matvec,
                         multiply, require, transpose, vector, verify)
from lattice_census import SCHEMA, admit, summarize


def box_shells(G, a, b, bound, max_nodes=1_000_000):
    """Independent complete rectangular enumeration with rational bounds."""
    inv = inverse(G)
    radii = [isqrt((bound * inv[i][i]).numerator // (bound * inv[i][i]).denominator)
             for i in range(len(G))]
    ranges = [range(-r + (ai + r) % 2, r + 1, 2) for r, ai in zip(radii, a)]
    counts, coeff = defaultdict(int), defaultdict(int)
    for nodes, y in enumerate(product(*ranges), 1):
        if nodes > max_nodes:
            raise EnumerationLimit("independent box exceeded node limit")
        q = dot(y, matvec(G, y))
        if q > bound:
            continue
        counts[q] += 1
        x = [(yi - ai) // 2 for yi, ai in zip(y, a)]
        coeff[q] += (-1) ** (dot(x, b) % 2)
    return [{"norm4": q, "norm": str(Fraction(q, 4)), "vector_count": counts[q],
             "signed_coefficient": coeff[q]} for q in sorted(counts)]


def box_vectors(G, bound, max_nodes):
    inv = inverse(G)
    radii = [isqrt((bound * inv[i][i]).numerator // (bound * inv[i][i]).denominator)
             for i in range(len(G))]
    vectors = []
    for nodes, v in enumerate(product(*(range(-r, r + 1) for r in radii)), 1):
        if nodes > max_nodes:
            raise EnumerationLimit("independent group box exceeded node limit")
        if dot(v, matvec(G, v)) <= bound:
            vectors.append(v)
    return vectors


def box_isometries(G, H, max_nodes):
    n = len(G)
    if n != len(H):
        return []
    vectors = box_vectors(H, max(G[i][i] for i in range(n)), max_nodes)
    columns = [[v for v in vectors if dot(v, matvec(H, v)) == G[i][i]] for i in range(n)]
    matches = []
    for nodes, cols in enumerate(product(*columns), 1):
        if nodes > max_nodes:
            raise EnumerationLimit("independent basis-image search exceeded node limit")
        T = transpose(cols)
        det = T[0][0] if n == 1 else T[0][0] * T[1][1] - T[1][0] * T[0][1]
        if abs(det) == 1 and multiply(multiply(transpose(T), H), T) == G:
            matches.append(T)
    return sorted(matches)


def _replay(document, max_nodes, backend):
    require(isinstance(document, dict) and document.get("schema") == SCHEMA, "unknown census schema")
    require(backend in ("box", "ldl"), "unknown backend")
    integer(max_nodes, "max_nodes")
    require(max_nodes > 0, "max_nodes must be positive")
    domain = document["domain"]
    require(domain["ranks"] == [1, 2], "wrong census ranks")
    D = integer(domain["determinant_bound"], "determinant_bound")
    require(1 <= D <= 24, "wrong determinant domain")
    bound = integer(document["run"]["norm4_bound"], "norm4_bound")
    integer(document["run"]["max_nodes"], "builder max_nodes")
    require(bound >= 0 and document["run"]["max_nodes"] > 0, "invalid run limits")
    # Independently reconstruct the finite reduction domain.
    expected = [[[d]] for d in range(1, D + 1)]
    for a in range(1, D + 1):
        if 3 * a * a > 4 * D:
            continue
        for b in range(a + 1):
            if 2 * b > a:
                continue
            for c in range(a, D + b * b + 1):
                if 1 <= a * c - b * b <= D:
                    expected.append([[a, b], [b, c]])
    generation, rows = document["generation"], document["lattices"]
    require(isinstance(generation, list) and isinstance(rows, list), "census lists required")
    require([g["G"] for g in generation] == expected, "generation coverage mismatch")
    seen = set()
    for g in generation:
        G = g["G"]
        admit(G, D)
        target = integer(g["representative"], "representative")
        require(0 <= target < len(rows), "bad representative reference")
        U = matrix(g["U"], "U", len(G))
        require(multiply(multiply(transpose(U), rows[target]["G"]), U) == G,
                "invalid generation basis change")
        require(all(x.denominator == 1 for row in inverse(U) for x in row), "nonunimodular basis change")
        if target not in seen:
            require(target == len(seen) and G == rows[target]["G"] and U == identity(len(G)),
                    "representatives must be retained in generation order")
            seen.add(target)
        require(type(g["uniqueness_complete"]) is bool, "uniqueness flag must be boolean")
        failures = g["incomplete_comparisons"]
        require(isinstance(failures, list) and len(set(failures)) == len(failures), "bad incomplete comparisons")
        for index in failures:
            integer(index, "comparison index")
            require(0 <= index < target and len(rows[index]["G"]) == len(G)
                    and rows[index]["determinant"] == admit(G, D), "bad incomplete comparison reference")
        require(g["uniqueness_complete"] == (not failures), "inconsistent uniqueness flag")
    require(seen == set(range(len(rows))), "unused or missing representative")
    uniqueness_replayed = True
    for i, row in enumerate(rows):
        G = row["G"]
        d = admit(G, D)
        require(type(row["determinant"]) is int and row["determinant"] == d, "determinant mismatch")
        digest = hashlib.sha256(json.dumps(G, separators=(",", ":")).encode("ascii")).hexdigest()
        require(row["gram_sha256"] == digest, "Gram hash mismatch")
        origin = next(g for g in generation if g["representative"] == i)
        for j in range(i):
            H = rows[j]["G"]
            if len(H) != len(G) or rows[j]["determinant"] != d:
                continue
            try:
                matches = box_isometries(G, H, max_nodes)
            except EnumerationLimit:
                uniqueness_replayed = False
                continue
            require(not matches or not origin["uniqueness_complete"], "isometric duplicate with uniqueness claim")
    replay_rows, replay_unresolved = [], []
    for i, row in enumerate(rows):
        G, n = row["G"], len(row["G"])
        group = row["automorphisms"]
        group_checked = False
        try:
            actual_group = box_isometries(G, G, max_nodes)
            if group is not None:
                require(isinstance(group, list), "group must be a list")
                for T in group:
                    matrix(T, "automorphism", n)
                require(group == actual_group, "full automorphism list mismatch")
                group_checked = True
            minimum = min(dot(v, matvec(G, v)) for v in box_vectors(G, G[0][0], max_nodes) if any(v))
            if row["minimum_squared_norm"] is not None:
                require(type(row["minimum_squared_norm"]) is int and row["minimum_squared_norm"] == minimum,
                        "minimum norm mismatch")
        except EnumerationLimit as error:
            replay_unresolved.append({"lattice": i, "kind": "group/minimum", "reason": str(error)})
        require(group is not None or isinstance(row["group_error"], str), "missing incomplete-group reason")
        pairs = row["pairs"]
        expected_pairs = [(list(a), list(b)) for a in product((0, 1), repeat=n)
                          for b in product((0, 1), repeat=n) if dot(a, b) % 2 == 0]
        require([(p["a"], p["b"]) for p in pairs] == expected_pairs, "even characteristic coverage mismatch")
        pair_results = []
        for p in pairs:
            a, b = vector(p["a"], "a", n), vector(p["b"], "b", n)
            if group is None:
                require(p["stabilizer"] is None, "stabilizer claimed without complete group")
            elif group_checked:
                indices, epsilon = [], []
                for k, T in enumerate(group):
                    dual = transpose(inverse(T))
                    if all((x - y) % 2 == 0 for x, y in zip(matvec(T, a), a)) and all(
                            (x - y) % 2 == 0 for x, y in zip(matvec(dual, b), b)):
                        shift = [(x - y) / 2 for x, y in zip(matvec(dual, b), b)]
                        require(all(x.denominator == 1 for x in shift), "nonintegral dual shift")
                        indices.append(k)
                        epsilon.append(int(dot(a, shift)) % 2)
                require(p["stabilizer"] == {"indices": indices, "epsilon": epsilon}, "stabilizer mismatch")
            evidence, verdict = p["evidence"], p["verdict"]
            if verdict == "proved_zero":
                require(evidence["kind"] == "symmetry", "zero lacks identity proof")
                result = verify({"schema": "work7-theta-v1", "G": G, "a": a, "b": b, "evidence": evidence}, max_nodes)
                require(result["verdict"] == "vanishes_by_symmetry", "bad zero evidence")
                pair_results.append({"a": a, "b": b, "verdict": verdict})
                continue
            require(verdict in ("proved_nonzero", "unresolved"), "unknown verdict")
            require(evidence["kind"] == ("shell" if verdict == "proved_nonzero" else "search"), "wrong evidence kind")
            if verdict == "unresolved":
                require(type(evidence["norm4_bound"]) is int and evidence["norm4_bound"] == bound, "wrong search bound")
            try:
                shells = (box_shells(G, a, b, bound, max_nodes) if backend == "box" else
                          complete_shells(G, a, b, bound, max_nodes)["shells"])
            except EnumerationLimit as error:
                replay_unresolved.append({"lattice": i, "a": a, "b": b, "kind": "shell", "reason": str(error)})
                pair_results.append({"a": a, "b": b, "verdict": "replay_unresolved"})
                continue
            if "enumeration" in p:
                require(p["enumeration"]["complete"] is True and p["enumeration"]["norm4_bound"] == bound,
                        "invalid completeness metadata")
                require(p["enumeration"]["shells"] == shells, "shell histogram mismatch")
            else:
                require(verdict == "unresolved" and p["enumeration_complete"] is False, "missing complete shells")
            if verdict == "proved_nonzero":
                q = integer(evidence["norm4"], "norm4")
                coeff = integer(evidence["coefficient"], "coefficient")
                require(0 <= q <= bound and coeff != 0 and any(s["norm4"] == q and s["signed_coefficient"] == coeff
                        for s in shells), "nonzero coefficient mismatch")
            elif "enumeration" in p:
                require(not any(s["signed_coefficient"] for s in shells), "unresolved claim contradicts complete shells")
            pair_results.append({"a": a, "b": b, "verdict": verdict, "shells": shells})
        replay_rows.append({"lattice": i, "G": G, "group_replayed": group_checked, "pairs": pair_results})
    require(document["summary"] == summarize(document), "summary mismatch")
    complete = uniqueness_replayed and not replay_unresolved
    return {"schema": "work7-rank12-census-result-v1", "backend": backend,
            "certificate_replay_complete": complete, "summary": document["summary"],
            "uniqueness_replayed": uniqueness_replayed, "replay_limits": replay_unresolved,
            "verdict": "unresolved" if not complete or document["summary"]["unresolved"]
            or not document["summary"]["groups_complete"] or not document["summary"]["uniqueness_complete"]
            else "resolved", "lattices": replay_rows}


def replay(document, max_nodes=1_000_000, backend="box"):
    try:
        return _replay(document, max_nodes, backend)
    except (KeyError, TypeError, IndexError, ZeroDivisionError, AttributeError) as error:
        raise InvalidCertificate(f"malformed census: {error}") from error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pack", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-nodes", type=int, default=1_000_000)
    parser.add_argument("--backend", choices=("box", "ldl"), default="box")
    args = parser.parse_args()
    try:
        raw = args.pack.read_bytes()
        result = replay(json.loads(raw), args.max_nodes, args.backend)
        result["pack_sha256"] = hashlib.sha256(raw).hexdigest()
    except (InvalidCertificate, OSError, json.JSONDecodeError) as error:
        result = {"verdict": "invalid_certificate", "reason": str(error)}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "lattices"}))
    return {"resolved": 0, "invalid_certificate": 1, "unresolved": 2}[result["verdict"]]


if __name__ == "__main__":
    raise SystemExit(main())
