"""Exact certificates for half-integral lattice theta functions (stdlib only).

Columns a and b are coordinates of 2xi in a lattice basis and 2delta
in its dual basis. Coefficients are integers after removing the phase i**(a.b).
See PROOFS.md for the mathematical conventions and enumeration proof.
"""

import argparse
from collections import defaultdict
from fractions import Fraction
import hashlib
import json
from math import isqrt
from pathlib import Path


class InvalidCertificate(ValueError):
    pass


class EnumerationLimit(RuntimeError):
    pass


def require(condition, message):
    if not condition:
        raise InvalidCertificate(message)


def integer(value, label):
    require(type(value) is int, f"{label} must be an integer")
    return value


def matrix(value, label, n=None):
    require(isinstance(value, list) and bool(value), f"{label} must be nonempty")
    n = len(value) if n is None else n
    require(len(value) == n, f"{label} has wrong dimension")
    for row in value:
        require(isinstance(row, list) and len(row) == n, f"{label} must be square")
        for entry in row:
            integer(entry, label)
    return value


def vector(value, label, n):
    require(isinstance(value, list) and len(value) == n, f"{label} has wrong dimension")
    for entry in value:
        integer(entry, label)
    return value


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def transpose(A):
    return [list(row) for row in zip(*A)]


def multiply(A, B):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*B)] for row in A]


def matvec(A, v):
    return [sum(x * y for x, y in zip(row, v)) for row in A]


def dot(u, v):
    return sum(x * y for x, y in zip(u, v))


def inverse(A):
    n = len(A)
    aug = [[Fraction(x) for x in row + unit] for row, unit in zip(A, identity(n))]
    for j in range(n):
        pivot = next((i for i in range(j, n) if aug[i][j]), None)
        require(pivot is not None, "matrix is singular")
        aug[j], aug[pivot] = aug[pivot], aug[j]
        scale = aug[j][j]
        aug[j] = [x / scale for x in aug[j]]
        for i in range(n):
            if i != j:
                scale = aug[i][j]
                aug[i] = [x - scale * y for x, y in zip(aug[i], aug[j])]
    return [row[n:] for row in aug]


def ldl(G):
    """G = L D L^t over Q; positive pivots certify positive definiteness."""
    require(G == transpose(G), "G must be symmetric")
    n = len(G)
    L = [[Fraction(x) for x in row] for row in identity(n)]
    D = []
    for i in range(n):
        pivot = Fraction(G[i][i]) - sum(L[i][k] ** 2 * D[k] for k in range(i))
        require(pivot > 0, "G must be positive definite")
        D.append(pivot)
        for j in range(i + 1, n):
            L[j][i] = (G[j][i] - sum(L[j][k] * L[i][k] * D[k] for k in range(i))) / pivot
    return L, D


def matrix_power(A, exponent):
    require(type(exponent) is int and exponent >= 0, "exponent must be nonnegative")
    result = identity(len(A))
    while exponent:
        if exponent & 1:
            result = multiply(result, A)
        A = multiply(A, A)
        exponent //= 2
    return result


def verify_order(T, order):
    integer(order, "order")
    require(order > 0, "order must be positive")
    I = identity(len(T))
    require(matrix_power(T, order) == I, "claimed order does not annihilate T")
    remaining, divisor = order, 2
    while divisor * divisor <= remaining:
        if remaining % divisor == 0:
            require(matrix_power(T, order // divisor) != I, "claimed order is not minimal")
            while remaining % divisor == 0:
                remaining //= divisor
        divisor += 1
    if remaining > 1:
        require(matrix_power(T, order // remaining) != I, "claimed order is not minimal")


def symmetry_data(G, a, b, T):
    n = len(G)
    matrix(T, "T", n)
    require(multiply(multiply(transpose(T), G), T) == G, "T does not preserve G")
    inv = inverse(T)
    require(all(x.denominator == 1 for row in inv for x in row), "T is not unimodular")
    dual = transpose([[int(x) for x in row] for row in inv])
    da = [x - y for x, y in zip(matvec(T, a), a)]
    db = [x - y for x, y in zip(matvec(dual, b), b)]
    require(all(x % 2 == 0 for x in da), "T does not stabilize a mod 2")
    require(all(x % 2 == 0 for x in db), "T does not stabilize b under inverse transpose mod 2")
    mu = [x // 2 for x in db]
    return {"epsilon": dot(a, mu) % 2, "dual_shift": mu,
            "primal_shift": [x // 2 for x in da]}


def ceil_fraction(x):
    return -((-x.numerator) // x.denominator)


def ceil_sqrt(x):
    require(x >= 0, "negative square-root bound")
    r = isqrt(x.numerator // x.denominator)
    return r if r * r == x else r + 1


def complete_shells(G, a, b, norm4_bound, max_nodes=1_000_000, vector_observer=None):
    """Enumerate ALL y=2x+a with y^t G y <= bound, by exact LDL recursion.

    No partial coefficients are returned if the resource limit is reached.
    Optional vector_observer(y_tuple, norm4) supports private batch histograms.
    A caller must discard its observer data if enumeration raises a limit.
    The function accepts validated square integer G and integer a,b.
    """
    integer(norm4_bound, "norm4_bound")
    require(norm4_bound >= 0, "norm4_bound must be nonnegative")
    integer(max_nodes, "max_nodes")
    require(max_nodes > 0, "max_nodes must be positive")
    L, D = ldl(G)
    n = len(G)
    y = [0] * n
    counts, coefficients = defaultdict(int), defaultdict(int)
    nodes = 0

    def visit(i, remaining):
        nonlocal nodes
        nodes += 1
        if nodes > max_nodes:
            raise EnumerationLimit(f"enumeration exceeded {max_nodes} nodes")
        if i < 0:
            norm4 = dot(y, matvec(G, y))
            require(0 <= norm4 <= norm4_bound, "internal enumeration error")
            x = [(yi - ai) // 2 for yi, ai in zip(y, a)]
            counts[norm4] += 1
            coefficients[norm4] += -1 if dot(x, b) % 2 else 1
            if vector_observer is not None:
                vector_observer(tuple(y), norm4)
            return
        center = sum(L[j][i] * y[j] for j in range(i + 1, n))
        radius = ceil_sqrt(remaining / D[i])
        lo = ceil_fraction(-center - radius)
        hi = (-center + radius).__floor__()
        lo += (a[i] - lo) % 2
        for yi in range(lo, hi + 1, 2):
            nodes += 1
            if nodes > max_nodes:
                raise EnumerationLimit(f"enumeration exceeded {max_nodes} nodes")
            used = D[i] * (yi + center) ** 2
            if used <= remaining:
                y[i] = yi
                visit(i - 1, remaining - used)

    visit(n - 1, Fraction(norm4_bound))
    return {"algorithm": "exact-rational-ldl-v1", "complete": True,
            "norm4_bound": norm4_bound, "visited_nodes": nodes,
            "shells": [{"norm4": q, "norm": str(Fraction(q, 4)),
                        "vector_count": counts[q], "signed_coefficient": coefficients[q]}
                       for q in sorted(counts)]}


def verify(document, max_nodes=1_000_000):
    require(isinstance(document, dict), "certificate must be an object")
    require(document.get("schema") == "work7-theta-v1", "unknown certificate schema")
    G = matrix(document.get("G"), "G")
    n = len(G)
    a = vector(document.get("a"), "a", n)
    b = vector(document.get("b"), "b", n)
    _, D = ldl(G)
    determinant = Fraction(1)
    for d in D:
        determinant *= d
    raw = json.dumps(G, separators=(",", ":"), ensure_ascii=True).encode("ascii")
    result = {"schema": "work7-theta-result-v1", "name": document.get("name", "unnamed"),
              "rank": n, "gram_sha256": hashlib.sha256(raw).hexdigest(),
              "determinant": int(determinant), "a": a, "b": b,
              "pairing": dot(a, b), "phase_i_exponent_mod4": dot(a, b) % 4,
              "basis_convention": "a primal; b dual; T acts on primal columns"}
    evidence = document.get("evidence")
    require(isinstance(evidence, dict), "evidence must be an object")
    kind = evidence.get("kind")
    if kind == "symmetry":
        data = symmetry_data(G, a, b, evidence.get("T"))
        require(data["epsilon"] == 1, "epsilon is zero: this is not a vanishing certificate")
        if "order" in evidence:
            verify_order(evidence["T"], evidence["order"])
            data["witness_order"] = evidence["order"]
        result.update(verdict="vanishes_by_symmetry", symmetry=data,
                      full_automorphism_group_verified=False)
    elif kind in ("shell", "search"):
        bound = evidence.get("norm4" if kind == "shell" else "norm4_bound")
        integer(bound, "norm4 bound")
        require(bound >= 0, "norm4 bound must be nonnegative")
        expected = None
        if kind == "shell":
            expected = integer(evidence.get("coefficient"), "coefficient")
            require(expected != 0, "a nonvanishing certificate requires a nonzero coefficient")
        try:
            enumeration = complete_shells(G, a, b, bound, max_nodes)
        except EnumerationLimit as error:
            result.update(verdict="unresolved", reason=str(error), enumeration_complete=False)
            return result
        result["enumeration"] = enumeration
        if kind == "shell":
            actual = next((s["signed_coefficient"] for s in enumeration["shells"]
                           if s["norm4"] == bound), 0)
            require(actual == expected, f"shell coefficient mismatch: expected {expected}, found {actual}")
            result.update(verdict="does_not_vanish", nonzero_norm4=bound,
                          signed_coefficient=actual)
        else:
            witness = next((s for s in enumeration["shells"] if s["signed_coefficient"]), None)
            if witness:
                result.update(verdict="does_not_vanish", nonzero_norm4=witness["norm4"],
                              signed_coefficient=witness["signed_coefficient"])
            else:
                result.update(verdict="unresolved", reason="all complete coefficients through the bound are zero")
    else:
        raise InvalidCertificate("unknown evidence kind")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificates", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, help="save a JSON array of results")
    parser.add_argument("--max-nodes", type=int, default=1_000_000)
    args = parser.parse_args()
    results, invalid, unresolved = [], False, False
    for path in args.certificates:
        try:
            with path.open(encoding="utf-8") as stream:
                result = verify(json.load(stream), args.max_nodes)
        except (InvalidCertificate, OSError, json.JSONDecodeError) as error:
            result = {"verdict": "invalid_certificate", "reason": str(error)}
            invalid = True
        result["certificate"] = str(path)
        unresolved |= result["verdict"] == "unresolved"
        results.append(result)
        print(json.dumps(result, ensure_ascii=True))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    return 1 if invalid else (2 if unresolved else 0)


if __name__ == "__main__":
    raise SystemExit(main())
