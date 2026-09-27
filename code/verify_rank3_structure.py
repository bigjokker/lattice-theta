"""Independent box and unimodular replay of rank-three structural certificates."""

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import combinations, product
import json
from math import gcd
from pathlib import Path

from exact_theta import (EnumerationLimit, InvalidCertificate, dot, identity, integer,
                         inverse, ldl, matrix, matvec, multiply, require, symmetry_data,
                         transpose, vector)
from verify_lattice_census import box_shells, box_vectors


def det(A):
    return sum(A[0][j] * (A[1][(j+1) % 3]*A[2][(j+2) % 3]
                         - A[1][(j+2) % 3]*A[2][(j+1) % 3]) for j in range(3))


def splitting_vectors(G, max_nodes):
    bound = max(G[i][i] for i in range(3))
    points = box_vectors(G, bound, max_nodes)
    candidates = []
    for v in points:
        if gcd(gcd(v[0], v[1]), v[2]) != 1:
            continue
        d = dot(v, matvec(G, v))
        if d and all(Fraction(x, d).denominator == 1 for x in matvec(G, v)):
            candidates.append(v)
    return candidates


def check_split(G, data):
    v = vector(data["v"], "splitting vector", 3)
    d = integer(data["norm"], "summand norm")
    require(gcd(gcd(v[0], v[1]), v[2]) == 1 and d == dot(v, matvec(G, v)) > 0,
            "nonprimitive or wrong-norm summand")
    lam = vector(data["projection_row"], "projection row", 3)
    require(lam == [Fraction(x, d) for x in matvec(G, v)] and dot(lam, v) == 1, "wrong integral projection")
    B = matrix(data["basis"], "split basis", 3)
    require(abs(det(B)) == 1 and transpose(B)[0] == v, "complement is not saturated/unimodular")
    H = multiply(multiply(transpose(B), G), B)
    require(data["block_gram"] == H and H[0][0] == d
            and H[0][1] == H[0][2] == H[1][0] == H[2][0] == 0, "block identity mismatch")
    A, C, E = H[1][1], H[1][2], H[2][2]
    require(1 <= A <= E and 0 <= 2*C <= A, "binary block is not reduced")
    kind = "three_lines" if C == 0 else "line_and_indecomposable_binary"
    require(data["type"] == kind, "wrong orthogonal decomposition type")
    return kind


def check_superbase(G, data, max_nodes):
    S = [[-1,-1,-1], [1,0,0], [0,1,0], [0,0,1]]
    require(isinstance(data["steps"], list), "Selling steps must be a list")
    if len(data["steps"]) > max_nodes:
        raise EnumerationLimit("Selling trace exceeds replay limit")
    for step in data["steps"]:
        require(isinstance(step, list) and len(step) == 2, "wrong Selling step")
        i, j = (integer(x, "Selling index") for x in step)
        require(0 <= i < j < 4, "bad Selling indices")
        epsilon = dot(S[i], matvec(G, S[j]))
        require(epsilon > 0, "Selling step does not reduce an acute pair")
        energy = sum(dot(v, matvec(G, v)) for v in S)
        old = S[i][:]
        S = [[-x for x in old] if k == i else S[k][:] if k == j
             else [x+y for x, y in zip(S[k], old)] for k in range(4)]
        require(sum(dot(v, matvec(G, v)) for v in S) == energy-2*epsilon, "energy decrease mismatch")
    require(data["vectors"] == S and all(sum(v[j] for v in S) == 0 for j in range(3)), "superbase mismatch")
    for omit in range(4):
        require(abs(det(transpose([S[i] for i in range(4) if i != omit]))) == 1, "superbase has wrong index")
    W = [[0 if i == j else -dot(S[i], matvec(G, S[j])) for j in range(4)] for i in range(4)]
    require(data["conorms"] == W and all(v >= 0 for row in W for v in row), "wrong/negative conorm")
    def reachable(exclude=None):
        vertices = set(range(4)) - ({exclude} if exclude is not None else set())
        seen = {min(vertices)}
        while True:
            larger = seen | {j for j in vertices for i in seen if W[i][j] > 0}
            if larger == seen:
                return seen == vertices
            seen = larger
    require(reachable() and all(reachable(i) for i in range(4)), "positive graph disconnected/articulated")
    edges = sum(W[i][j] > 0 for i, j in combinations(range(4), 2))
    require(edges in (4, 5, 6), "invalid indecomposable graph type")
    return {4: "cycle", 5: "diamond", 6: "complete"}[edges]


def _replay(document, census, census_sha256, max_nodes):
    require(document["schema"] == "work7-rank3-structure-v1", "wrong structural schema")
    require(document["census_sha256"] == census_sha256, "source census hash mismatch")
    require(document["domain"] == {"ranks": [3], "determinant_bound": 24}, "wrong structural domain")
    require(document["proof_note"] == "RANK3-STRUCTURE.md", "wrong proof convention")
    integer(max_nodes, "max_nodes")
    integer(document["max_nodes"], "builder max_nodes")
    require(max_nodes > 0 and document["max_nodes"] > 0, "invalid resource limit")
    rows = document["rows"]
    require([r["G"] for r in rows] == [r["G"] for r in census["lattices"]], "source class coverage mismatch")
    require([r["lattice"] for r in rows] == list(range(len(rows))), "class order mismatch")
    require(isinstance(document["limits"], list), "limits must be a list")
    indexed_limits = {r["lattice"]: r["reason"] for r in document["limits"]}
    require(len(indexed_limits) == len(document["limits"]), "duplicate builder limits")
    for i, reason in indexed_limits.items():
        integer(i, "limited lattice")
        require(0 <= i < len(rows) and isinstance(reason, str), "bad builder limit")
    types, graphs, cases, verdicts = Counter(), Counter(), Counter(), Counter()
    report, replay_limits = [], []
    for i, row in enumerate(rows):
        G = matrix(row["G"], "G", 3)
        _, pivots = ldl(G)
        require(1 <= pivots[0]*pivots[1]*pivots[2] <= 24, "wrong determinant")
        decomposition = row["decomposition"]
        dtype = (decomposition["split"]["type"] if decomposition and decomposition["verdict"] == "decomposable"
                 else decomposition["verdict"] if decomposition else "unresolved")
        types[dtype] += 1
        if row["superbase"]:
            # Recorded graph types are counted even when a later builder limit occurred.
            W = row["superbase"]["conorms"]
            graphs[{4: "cycle", 5: "diamond", 6: "complete"}[sum(W[j][k] > 0 for j, k in combinations(range(4), 2))]] += 1
        for p in row["pairs"]:
            verdicts[p["verdict"]] += 1
            if "shell" in p:
                cases[p["shell"]["case"]] += 1
        if i in indexed_limits:
            require(row["pairs"] == [], "limited row has partial characteristic claims")
            report.append({"lattice": i, "verdict": "builder_unresolved"})
            continue
        require(decomposition and decomposition["sphere_complete"] is True
                and decomposition["search_bound"] == max(G[j][j] for j in range(3)), "incomplete splitting search")
        try:
            candidates = splitting_vectors(G, max_nodes)
            if decomposition["verdict"] == "decomposable":
                kind = check_split(G, decomposition["split"])
                require(tuple(decomposition["split"]["v"]) in candidates, "split absent from complete sphere")
                require(row["superbase"] is None, "unexpected positive superbase proof on split row")
            else:
                require(decomposition["verdict"] == "indecomposable" and not candidates, "false indecomposability")
                kind = "indecomposable"
                check_superbase(G, row["superbase"], max_nodes)
            pairs = [(list(a), list(b)) for a in product((0,1), repeat=3) for b in product((0,1), repeat=3)]
            require([(p["a"], p["b"]) for p in row["pairs"]] == pairs, "characteristic coverage mismatch")
            for p in row["pairs"]:
                a, b = vector(p["a"], "a", 3), vector(p["b"], "b", 3)
                odd = dot(a, b) % 2
                if odd:
                    require(p["proof"] == "odd-parity" and p["verdict"] == "proved_zero", "wrong odd verdict")
                    require(p["T"] == [[-int(j == k) for k in range(3)] for j in range(3)], "odd witness is not -I")
                elif kind != "indecomposable":
                    B = decomposition["split"]["basis"]
                    ap = [int(x) % 2 for x in matvec(inverse(B), a)]
                    bp = [x % 2 for x in matvec(transpose(B), b)]
                    require(p["factor_a"] == ap and p["factor_b"] == bp, "wrong factor transport")
                    zero = (any(ap[j]*bp[j] for j in range(3)) if kind == "three_lines"
                            else (ap[0]*bp[0] % 2 or (ap[1]*bp[1]+ap[2]*bp[2]) % 2))
                    require(p["verdict"] == ("proved_zero" if zero else "proved_nonzero")
                            and p["proof"] == ("odd-factor" if zero else "rank12-product"), "wrong factor verdict")
                else:
                    require(p["proof"] == "obtuse-graph-shell" and p["verdict"] == "proved_nonzero", "wrong graph verdict")
                    shell = p["shell"]
                    N = integer(shell["norm4"], "shell norm4")
                    coeff = integer(shell["coefficient"], "coefficient")
                    require(N >= 0 and coeff != 0 and shell["complete_by"] == "rank3-obtuse-graph-shell-v1", "invalid shell claim")
                    supplied = shell["selected_vectors"]
                    require(isinstance(supplied, list) and len({tuple(v) for v in supplied}) == len(supplied), "repeated selected vector")
                    selected_coeff = 0
                    for y in supplied:
                        vector(y, "selected vector", 3)
                        require(dot(y, matvec(G, y)) == N and all((yi-ai) % 2 == 0 for yi, ai in zip(y, a)), "invalid selected vector")
                        selected_coeff += (-1)**(dot([(yi-ai)//2 for yi, ai in zip(y, a)], b) % 2)
                    require(selected_coeff == coeff, "selected signs mismatch")
                    complete_shells = box_shells(G, a, b, N, max_nodes)
                    actual = next((s["signed_coefficient"] for s in complete_shells if s["norm4"] == N), 0)
                    require(actual == coeff, "independent complete coefficient mismatch")
                if p["verdict"] == "proved_zero":
                    require(symmetry_data(G, a, b, p["T"])["epsilon"] == 1, "invalid structural witness")
                if not odd:
                    source = next(q for q in census["lattices"][i]["pairs"] if q["a"] == a and q["b"] == b)
                    require(source["verdict"] == p["verdict"], "saved census comparison mismatch")
            report.append({"lattice": i, "structure": kind, "splitting_sphere_complete": True,
                           "pairs_replayed": len(row["pairs"]), "verdict": "verified"})
        except EnumerationLimit as error:
            replay_limits.append({"lattice": i, "reason": str(error)})
            report.append({"lattice": i, "verdict": "replay_unresolved"})
    expected = {"classes": len(rows), "structure_types": dict(types),
                "even_pairs": sum(dot(p["a"], p["b"]) % 2 == 0 for r in rows for p in r["pairs"]),
                "all_pairs": sum(len(r["pairs"]) for r in rows), "proved_zero": verdicts["proved_zero"],
                "proved_nonzero": verdicts["proved_nonzero"], "complete": not document["limits"],
                "positive_graphs": dict(graphs), "symbolic_shell_cases": dict(cases)}
    require(document["summary"] == expected, "structural summary mismatch")
    return {"schema": "work7-rank3-structure-result-v1", "summary": expected, "backend": "box",
            "certificate_replay_complete": not replay_limits, "replay_limits": replay_limits,
            "verdict": "resolved" if not replay_limits and not document["limits"] else "unresolved",
            "census_sha256": census_sha256, "rows": report,
            "symbolic_theorem_externally_reviewed": False}


def replay(document, census, census_sha256, max_nodes=1_000_000):
    try:
        return _replay(document, census, census_sha256, max_nodes)
    except (KeyError, TypeError, IndexError, AttributeError, ValueError, ZeroDivisionError) as error:
        raise InvalidCertificate(f"malformed structural certificate: {error}") from error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pack", type=Path)
    parser.add_argument("--census", type=Path, default=Path("data/census-rank3/det24-bound32.json"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-nodes", type=int, default=1_000_000)
    args = parser.parse_args()
    try:
        raw, source = args.pack.read_bytes(), args.census.read_bytes()
        result = replay(json.loads(raw), json.loads(source), hashlib.sha256(source).hexdigest(), args.max_nodes)
        result["pack_sha256"] = hashlib.sha256(raw).hexdigest()
    except (InvalidCertificate, OSError, json.JSONDecodeError) as error:
        result = {"verdict": "invalid_certificate", "reason": str(error)}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "rows"}))
    return {"resolved": 0, "unresolved": 2, "invalid_certificate": 1}[result["verdict"]]


if __name__ == "__main__":
    raise SystemExit(main())
