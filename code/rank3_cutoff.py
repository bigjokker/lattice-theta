"""Build/replay exact rank-three modular identity certificates; RANK3-CUTOFF.md."""

import argparse
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path

from exact_theta import (EnumerationLimit, InvalidCertificate, complete_shells, dot,
                         integer, inverse, ldl, matrix, multiply, require, transpose, vector)
from rank3_census import admit, modular_data
from verify_lattice_census import box_shells

SCHEMA = "work7-rank3-cutoff-v1"
THEOREM = "rank3-fourth-power-complex-identity-v1"


def factor(n):
    integer(n, "factor input")
    require(n > 0, "positive factor input required")
    result, p = [], 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            result.append([p, e])
        p += 1 if p == 2 else 2
    if n > 1:
        result.append([n, 1])
    return result


def index_gamma0(M):
    result = 1
    for p, e in factor(M):
        result *= p ** (e - 1) * (p + 1)
    return result


def cutoff(N0):
    integer(N0, "N0")
    require(N0 > 0, "N0 must be positive")
    mu = index_gamma0(16 * N0)
    require(mu % 8 == 0, "index divisibility failure")
    return mu // 8


def inputs(G, a, b):
    admit(G)
    vector(a, "a", 3)
    vector(b, "b", 3)
    require(all(v in (0, 1) for v in a + b), "binary representatives required")
    require(any(b) and dot(a, b) % 2 == 0, "nonzero b and even pairing required")


def metadata(G, b):
    item = modular_data(G, b)
    item["cutoff"] = cutoff(item["N0"])
    item["level_factors"] = factor(item["unsigned_level"])
    item["index"] = index_gamma0(item["unsigned_level"])
    item["cleared_inverse"] = [[int(v * item["N0"]) for v in row] for row in inverse(item["G0"])]
    item["theorem"] = THEOREM
    return item


def verdict(shells, bound, B):
    nonzero = next((s for s in shells if s["signed_coefficient"]), None)
    if nonzero:
        return {"verdict": "proved_nonzero", "nonzero_norm4": nonzero["norm4"],
                "coefficient": nonzero["signed_coefficient"]}
    if bound >= B:
        return {"verdict": "proved_zero_by_modular_cutoff"}
    return {"verdict": "unresolved", "reason": "complete zero search below the inclusive cutoff"}


def build(G, a, b, norm4_bound=None, max_nodes=1_000_000):
    inputs(G, a, b)
    item = metadata(G, b)
    bound = item["cutoff"] if norm4_bound is None else integer(norm4_bound, "norm4_bound")
    integer(max_nodes, "max_nodes")
    require(bound >= 0 and max_nodes > 0, "invalid run limits")
    result = {"schema": SCHEMA, "G": G, "a": a, "b": b, "modular": item,
              "run": {"norm4_bound": bound, "max_nodes": max_nodes}}
    try:
        enumeration = complete_shells(G, a, b, bound, max_nodes)
        result["evidence"] = {"enumeration_complete": True, "enumeration": enumeration,
                              **verdict(enumeration["shells"], bound, item["cutoff"])}
    except EnumerationLimit as error:
        result["evidence"] = {"enumeration_complete": False, "verdict": "unresolved", "reason": str(error)}
    return result


def independent_modular_check(G, b, item):
    """Check any supplied kernel basis and recompute denominator/index data."""
    K = matrix(item["index_two_basis"], "K", 3)
    # A separate determinant expansion, without the rank-three builder.
    def determinant(A):
        return sum(A[0][j] * (A[1][(j+1) % 3] * A[2][(j+2) % 3]
                             - A[1][(j+2) % 3] * A[2][(j+1) % 3]) for j in range(3))
    require(abs(determinant(K)) == 2 and all(dot(col, b) % 2 == 0 for col in transpose(K)),
            "K does not span the index-two kernel")
    G0 = multiply(multiply(transpose(K), G), K)
    require(item["G0"] == G0, "G0 mismatch")
    _, pivots = ldl(G0)
    d0 = int(pivots[0] * pivots[1] * pivots[2])
    require(integer(item["d0"], "d0") == d0 == 4 * admit(G), "d0 mismatch")
    inv = inverse(G0)
    N0 = 1
    for row in inv:
        for v in row:
            N0 = N0 * v.denominator // gcd(N0, v.denominator)
    require(integer(item["N0"], "N0") == N0, "level mismatch")
    M = 16 * N0
    # Distinct-prime rational product, rather than the builder's prime-power formula.
    n, primes, p = M, [], 2
    while p * p <= n:
        if n % p == 0:
            primes.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        primes.append(n)
    index = Fraction(M)
    for p in primes:
        index *= Fraction(p + 1, p)
    require(index.denominator == 1 and int(index) % 8 == 0, "index arithmetic failure")
    expected_factors = []
    for p in primes:
        n, e = M, 0
        while n % p == 0:
            n //= p
            e += 1
        expected_factors.append([p, e])
    require(item["level_factors"] == expected_factors, "level factorization mismatch")
    require(integer(item["unsigned_level"], "unsigned_level") == M
            and integer(item["index"], "index") == int(index)
            and integer(item["cutoff"], "cutoff") == int(index) // 8, "cutoff mismatch")
    cleared = matrix(item["cleared_inverse"], "cleared inverse", 3)
    require(cleared == [[N0 * v for v in row] for row in inv], "cleared inverse mismatch")
    require(integer(item["character_discriminant"], "character discriminant") == 4 * d0
            and item["b"] == b and item["weight"] == "3/2"
            and item["coefficient_index"] == "norm4" and item["theorem"] == THEOREM,
            "modular convention mismatch")
    return int(index) // 8


def _replay(document, max_nodes, backend):
    require(document["schema"] == SCHEMA, "unknown cutoff schema")
    G, a, b = document["G"], document["a"], document["b"]
    inputs(G, a, b)
    B = independent_modular_check(G, b, document["modular"])
    bound = integer(document["run"]["norm4_bound"], "norm4_bound")
    nodes = integer(document["run"]["max_nodes"], "builder max_nodes")
    integer(max_nodes, "max_nodes")
    require(bound >= 0 and nodes > 0 and max_nodes > 0 and backend in ("box", "ldl"), "bad replay bounds/backend")
    evidence = document["evidence"]
    require(type(evidence["enumeration_complete"]) is bool, "invalid completeness flag")
    if not evidence["enumeration_complete"]:
        require(evidence["verdict"] == "unresolved" and "enumeration" not in evidence
                and isinstance(evidence["reason"], str), "interrupted evidence has a resolved/partial claim")
        return {"verdict": "unresolved", "certificate_replay_complete": True,
                "coefficient_replay_complete": False, "reason": "builder enumeration interrupted",
                "cutoff": B, "backend": backend}
    enumeration = evidence["enumeration"]
    require(enumeration["complete"] is True and integer(enumeration["norm4_bound"], "enumeration bound") == bound,
            "incomplete coefficient evidence")
    try:
        shells = box_shells(G, a, b, bound, max_nodes) if backend == "box" else complete_shells(G, a, b, bound, max_nodes)["shells"]
    except EnumerationLimit as error:
        return {"verdict": "unresolved", "certificate_replay_complete": False,
                "coefficient_replay_complete": False, "reason": str(error), "cutoff": B, "backend": backend}
    require(enumeration["shells"] == shells, "complete coefficient histogram mismatch")
    actual = verdict(shells, bound, B)
    for key, value in actual.items():
        require(evidence.get(key) == value, "coefficient verdict mismatch")
    return {**actual, "certificate_replay_complete": True, "coefficient_replay_complete": True,
            "cutoff": B, "norm4_bound": bound, "backend": backend, "shells": shells,
            "full_stabilizer_checked": False}


def replay(document, max_nodes=1_000_000, backend="box"):
    try:
        return _replay(document, max_nodes, backend)
    except (KeyError, TypeError, IndexError, AttributeError, ValueError, ZeroDivisionError) as error:
        raise InvalidCertificate(f"malformed cutoff certificate: {error}") from error


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("build", "verify"))
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--norm4-bound", type=int)
    parser.add_argument("--max-nodes", type=int, default=1_000_000)
    parser.add_argument("--backend", choices=("box", "ldl"), default="box")
    args = parser.parse_args()
    try:
        raw = args.input.read_bytes()
        document = json.loads(raw)
        if args.mode == "build":
            result = build(document["G"], document["a"], document["b"], args.norm4_bound, args.max_nodes)
            status = result["evidence"]["verdict"]
        else:
            require(args.norm4_bound is None, "verify uses the certificate's recorded bound")
            result = replay(document, args.max_nodes, args.backend)
            result["certificate_sha256"] = hashlib.sha256(raw).hexdigest()
            status = result["verdict"]
    except (InvalidCertificate, OSError, json.JSONDecodeError, KeyError) as error:
        result, status = {"verdict": "invalid_certificate", "reason": str(error)}, "invalid_certificate"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result))
    return 1 if status == "invalid_certificate" else 2 if status == "unresolved" else 0


if __name__ == "__main__":
    raise SystemExit(main())
