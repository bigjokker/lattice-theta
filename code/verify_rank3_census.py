"""Independent finite-domain, full-group and signed-shell replay for rank three."""

import argparse
from collections import Counter
import hashlib
from itertools import product
import json
from pathlib import Path

from exact_theta import (EnumerationLimit, InvalidCertificate, complete_shells, dot,
                         identity, integer, inverse, ldl, matrix, matvec, multiply,
                         require, symmetry_data, transpose, vector, verify)
from rank3_census import SCHEMA, summarize
from verify_lattice_census import box_shells, box_vectors


def det3(A):
    return sum(A[0][j] *
               (A[1][(j + 1) % 3] * A[2][(j + 2) % 3]
                - A[1][(j + 2) % 3] * A[2][(j + 1) % 3])
               for j in range(3))


def admission(G, H):
    matrix(G, "G", 3)
    _, pivots = ldl(G)
    d = pivots[0] * pivots[1] * pivots[2]
    require(d.denominator == 1 and 1 <= d <= H, "wrong Gram domain")
    return int(d)


def independent_candidates(H):
    """Scan bounded diagonals, rather than solve E from a chosen determinant."""
    found = []
    for a in range(1, H + 1):
        if 27 * a ** 3 > 64 * H:
            continue
        for b in range(-a, a + 1):
            if 2 * abs(b) > a:
                continue
            for c in range(a, 2 * H + a + 1):
                m = a * c - b * b
                if 3 * m * m > 4 * a * H:
                    continue
                for d in range(-a, a + 1):
                    if 2 * abs(d) > a:
                        continue
                    for f in range(-H - a, H + a + 1):
                        if 2 * abs(a * f - b * d) > m:
                            continue
                        for e in range(1, 2 * H + a + 1):
                            delta = m * e - c * d * d - a * f * f + 2 * b * d * f
                            if a * e - d * d >= m and 1 <= delta <= H:
                                found.append([[a, b, d], [b, c, f], [d, f, e]])
    return sorted(found)


def box_isometries(G, H, max_nodes):
    vectors = box_vectors(H, max(G[i][i] for i in range(3)), max_nodes)
    images = [[v for v in vectors if dot(v, matvec(H, v)) == G[i][i]] for i in range(3)]
    found = []
    # Unpruned exhaustive tuples deliberately differ from the builder.
    for nodes, cols in enumerate(product(*images), 1):
        if nodes > max_nodes:
            raise EnumerationLimit("independent rank-three basis-image limit")
        U = transpose(cols)
        if abs(det3(U)) == 1 and multiply(multiply(transpose(U), H), U) == G:
            found.append(U)
    return sorted(found)


def check_modular_inputs(G, inputs):
    expected_b = [list(b) for b in product((0, 1), repeat=3) if any(b)]
    require([item["b"] for item in inputs] == expected_b, "modular b coverage mismatch")
    for item in inputs:
        b, K = vector(item["b"], "modular b", 3), matrix(item["index_two_basis"], "index-two basis", 3)
        require(abs(det3(K)) == 2 and all(v % 2 == 0 for v in matvec(transpose(K), b)),
                "basis is not the index-two kernel")
        G0 = multiply(multiply(transpose(K), G), K)
        require(item["G0"] == G0 and item["d0"] == det3(G0) == 4 * det3(G), "wrong G0")
        inv = inverse(G0)
        # Incremental gcd arithmetic, independently of the builder's math.lcm.
        from math import gcd
        N = 1
        for row in inv:
            for v in row:
                N = N * v.denominator // gcd(N, v.denominator)
        require(item["N0"] == N and item["unsigned_level"] == 16 * N
                and item["character_discriminant"] == 4 * det3(G0)
                and item["weight"] == "3/2" and item["coefficient_index"] == "norm4"
                and item["cutoff"] is None, "incorrect or unproved modular metadata")


def _replay(document, max_nodes, backend):
    require(document["schema"] == SCHEMA and document["domain"]["ranks"] == [3], "wrong schema/domain")
    H = integer(document["domain"]["determinant_bound"], "determinant_bound")
    bound = integer(document["run"]["norm4_bound"], "norm4_bound")
    integer(max_nodes, "max_nodes")
    builder_nodes = integer(document["run"]["max_nodes"], "builder max_nodes")
    require(1 <= H <= 24 and bound >= 0 and max_nodes > 0 and builder_nodes > 0, "invalid bounds")
    require(backend in ("box", "ldl"), "invalid backend")
    generation, rows = document["generation"], document["lattices"]
    require([g["G"] for g in generation] == independent_candidates(H), "generation coverage mismatch")
    seen = set()
    for g in generation:
        G, target = g["G"], integer(g["representative"], "representative")
        admission(G, H)
        require(0 <= target < len(rows), "bad representative reference")
        U = matrix(g["U"], "basis change", 3)
        require(abs(det3(U)) == 1 and multiply(multiply(transpose(U), rows[target]["G"]), U) == G,
                "bad generation basis change")
        if target not in seen:
            require(target == len(seen) and G == rows[target]["G"] and U == identity(3), "retention order mismatch")
            seen.add(target)
        failures = g["incomplete_comparisons"]
        require(isinstance(failures, list) and len(set(failures)) == len(failures)
                and type(g["uniqueness_complete"]) is bool and g["uniqueness_complete"] == (not failures),
                "invalid uniqueness metadata")
        for k in failures:
            integer(k, "comparison index")
            require(0 <= k < target and rows[k]["determinant"] == admission(G, H), "bad failed comparison")
    require(seen == set(range(len(rows))), "unused representative")
    limits, replay_rows, uniqueness = [], [], True
    for i, row in enumerate(rows):
        G = row["G"]
        d = admission(G, H)
        require(row["determinant"] == d and type(row["determinant"]) is int, "wrong determinant")
        require(row["gram_sha256"] == hashlib.sha256(json.dumps(G, separators=(",", ":")).encode("ascii")).hexdigest(),
                "wrong Gram hash")
        origin = next(g for g in generation if g["representative"] == i)
        for j in range(i):
            if rows[j]["determinant"] != d:
                continue
            try:
                matches = box_isometries(G, rows[j]["G"], max_nodes)
                require(not matches or not origin["uniqueness_complete"], "duplicate with uniqueness claim")
            except EnumerationLimit as error:
                uniqueness = False
                limits.append({"lattice": i, "comparison": j, "reason": str(error)})
        group_checked, group = False, row["automorphisms"]
        try:
            actual = box_isometries(G, G, max_nodes)
            if group is not None:
                for T in group:
                    matrix(T, "automorphism", 3)
                require(group == actual, "full automorphism mismatch")
                group_checked = True
            minimum = min(dot(v, matvec(G, v)) for v in box_vectors(G, G[0][0], max_nodes) if any(v))
            require(row["minimum_squared_norm"] is None or (type(row["minimum_squared_norm"]) is int
                    and row["minimum_squared_norm"] == minimum), "minimum mismatch")
        except EnumerationLimit as error:
            limits.append({"lattice": i, "kind": "group/minimum", "reason": str(error)})
        require(group is not None or isinstance(row["group_error"], str), "missing group-limit reason")
        all_pairs = [(list(a), list(b)) for a in product((0, 1), repeat=3) for b in product((0, 1), repeat=3)]
        require([(p["a"], p["b"]) for p in row["pairs"]] == [(a, b) for a, b in all_pairs if dot(a, b) % 2 == 0],
                "even coverage mismatch")
        require(row["odd_pairs"] == [{"a": a, "b": b, "verdict": "proved_zero", "evidence": "minus_identity"}
                                      for a, b in all_pairs if dot(a, b) % 2], "odd coverage mismatch")
        check_modular_inputs(G, row["modular_inputs"])
        for p in row["odd_pairs"]:
            vector(p["a"], "odd a", 3)
            vector(p["b"], "odd b", 3)
            require(symmetry_data(G, p["a"], p["b"], [[-int(i == j) for j in range(3)] for i in range(3)])["epsilon"] == 1,
                    "odd parity witness mismatch")
        results = []
        for p in row["pairs"]:
            a, b, verdict, evidence = vector(p["a"], "a", 3), vector(p["b"], "b", 3), p["verdict"], p["evidence"]
            if group is None:
                require(p["stabilizer"] is None, "stabilizer without full group")
            elif group_checked:
                indices, signs = [], []
                for k, T in enumerate(group):
                    dual = transpose(inverse(T))
                    if all((x-y) % 2 == 0 for x, y in zip(matvec(T, a), a)) and all(
                            (x-y) % 2 == 0 for x, y in zip(matvec(dual, b), b)):
                        shift = [(x-y)/2 for x, y in zip(matvec(dual, b), b)]
                        require(all(v.denominator == 1 for v in shift), "nonintegral dual shift")
                        indices.append(k)
                        signs.append(int(dot(a, shift)) % 2)
                require(p["stabilizer"] == {"indices": indices, "epsilon": signs}, "stabilizer mismatch")
            require(verdict in ("proved_zero", "proved_nonzero", "unresolved"), "unknown verdict")
            if verdict == "proved_zero":
                require(evidence["kind"] == "symmetry", "zero without identity proof")
                checked = verify({"schema": "work7-theta-v1", "G": G, "a": a, "b": b, "evidence": evidence}, max_nodes)
                require(checked["verdict"] == "vanishes_by_symmetry", "invalid zero")
            else:
                require(evidence["kind"] == ("shell" if verdict == "proved_nonzero" else "search"), "wrong evidence")
                if verdict == "unresolved":
                    require(evidence["norm4_bound"] == bound, "wrong search bound")
                try:
                    shells = box_shells(G, a, b, bound, max_nodes) if backend == "box" else complete_shells(G, a, b, bound, max_nodes)["shells"]
                    if "enumeration" in p:
                        require(p["enumeration"]["complete"] is True and p["enumeration"]["norm4_bound"] == bound
                                and p["enumeration"]["shells"] == shells, "complete shell mismatch")
                    else:
                        require(verdict == "unresolved" and p["enumeration_complete"] is False, "missing shell evidence")
                    if verdict == "proved_nonzero":
                        N = integer(evidence["norm4"], "norm4")
                        c = integer(evidence["coefficient"], "coefficient")
                        require(0 <= N <= bound and c != 0 and any(s["norm4"] == N and s["signed_coefficient"] == c for s in shells),
                                "nonzero coefficient mismatch")
                    elif "enumeration" in p:
                        require(not any(s["signed_coefficient"] for s in shells), "unresolved contradicts shell")
                except EnumerationLimit as error:
                    limits.append({"lattice": i, "a": a, "b": b, "reason": str(error)})
            results.append({"a": a, "b": b, "verdict": verdict})
        replay_rows.append({"lattice": i, "G": G, "group_replayed": group_checked, "pairs": results})
    require(document["summary"] == summarize(document), "summary mismatch")
    complete = uniqueness and not limits
    resolved = complete and document["summary"]["symmetry_converse_proved_in_domain"]
    return {"schema": "work7-rank3-census-result-v1", "backend": backend,
            "certificate_replay_complete": complete, "uniqueness_replayed": uniqueness,
            "replay_limits": limits, "summary": document["summary"],
            "verdict": "resolved" if resolved else "unresolved", "lattices": replay_rows}


def replay(document, max_nodes=1_000_000, backend="box"):
    try:
        return _replay(document, max_nodes, backend)
    except (KeyError, TypeError, IndexError, ZeroDivisionError, AttributeError, ValueError) as error:
        raise InvalidCertificate(f"malformed rank-three census: {error}") from error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pack", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--backend", choices=("box", "ldl"), default="box")
    parser.add_argument("--max-nodes", type=int, default=1_000_000)
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
    return {"resolved": 0, "unresolved": 2, "invalid_certificate": 1}[result["verdict"]]


if __name__ == "__main__":
    raise SystemExit(main())
