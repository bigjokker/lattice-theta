"""Exact supplementary checks for docs/rank-four-stage-two-reductions.md.

These checks exercise the constructions in the proposed proofs. They do not
prove the conjecture or extend a census. Existing source data are read only.

Release layout (run from the repository root): inputs are read from data/ and
reports/ (a preserved .json.gz is read decompressed); write fresh outputs under
scratch/. This is the research-workspace script of 2026-10-21 with relocated
input paths; the preserved outputs under reports/ were produced by the
workspace copy, so paths recorded inside them are the original workspace paths.
"""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import json
from math import isqrt
from pathlib import Path

from exact_theta import (dot, identity, inverse, matvec, multiply,
                         symmetry_data, transpose)
from rank4_feasibility import determinant, isometries
from rank4_level4_search import gram_of, kernel_basis, lll_reduce
from rank4_minus_type import BITS, I, NEG_I, action, order
from rank4_parity_search import quadratic_value, span_basis

ROOT = Path(__file__).resolve().parent.parent


def _source_bytes(path):
    """Bytes of a preserved input in the release layout.

    A research-workspace path (examples/, results/) is mapped to data/ or
    reports/; when only the lossless .json.gz copy is present, its
    decompressed bytes are returned, so recorded hashes are unchanged.
    """
    p = Path(path)
    if not p.exists():
        s = str(path).replace("\\", "/")
        for old, new in (("examples/", "data/"), ("results/", "reports/")):
            if s.startswith(old):
                p = ROOT / (new + s[len(old):])
        if not p.exists() and not p.is_absolute():
            p = ROOT / p
    if p.exists():
        return p.read_bytes()
    compressed = Path(str(p) + ".gz")
    if compressed.exists():
        return gzip.decompress(compressed.read_bytes())
    raise FileNotFoundError(path)


def integral(A):
    assert all(Fraction(x).denominator == 1 for row in A for x in row)
    return [[int(x) for x in row] for row in A]


def eigensublattice(h, sign):
    """Saturated Z-basis of the sign eigenspace of an integral involution."""
    twice_projection = [[int(i == j) + sign * h[i][j]
                         for j in range(4)] for i in range(4)]
    generators = transpose(twice_projection)
    for v in BITS:
        w = matvec(twice_projection, v)
        if all(x % 2 == 0 for x in w):
            generators.append([x // 2 for x in w])
    return span_basis(generators)


def mixed_reflections(G, g):
    h = multiply(g, g)
    N, P = eigensublattice(h, -1), eigensublattice(h, 1)
    assert len(N[0]) == len(P[0]) == 2
    U, reduced = lll_reduce(gram_of(G, N))
    c = reduced[0][0]
    assert reduced == [[c, 0], [0, c]] and c % 2 == 0
    e = transpose(multiply(N, U))[0]
    f = matvec(g, e)
    B = transpose([e, f] + transpose(P))
    index = abs(determinant(B))
    assert index in (1, 2)
    Binv = inverse(B)
    T = integral(multiply(multiply(Binv, g), B))
    assert [row[:2] for row in T[:2]] == [[0, -1], [1, 0]]
    assert all(T[i][j] == T[j][i] == 0 for i in range(2) for j in range(2, 4))
    out = []
    for sign in (1, -1):
        R0 = [[0] * 4 for _ in range(4)]
        R0[0][0], R0[1][1] = sign, -sign
        for i in range(2, 4):
            for j in range(2, 4):
                R0[i][j] = T[i][j]
        R = integral(multiply(multiply(B, R0), Binv))
        assert multiply(R, R) == I and gram_of(G, R) == G
        out.append(R)
    return out, index


def mixed_and_mapping_checks(doc):
    counts = Counter()
    for row in doc["lattices"]:
        G, group = row["G"], row["automorphisms"]
        d0 = determinant(G)
        if d0 % 8 != 4:
            continue
        counts["kernel_classes"] += 1
        for g in group:
            if order(g) != 4 or multiply(g, g) == NEG_I:
                continue
            reflections, index = mixed_reflections(G, g)
            counts["mixed_order_four_elements"] += 1
            counts[f"eigenlattice_index_{index}"] += 1
            for v in BITS[1:]:
                w = action(g, v)
                if w == v or action(g, w) == v:
                    continue
                assert any(action(R, v) == w for R in reflections)
                counts["four_cycle_steps_with_constructed_reflection"] += 1

        d = d0 // 4
        if isqrt(d) ** 2 == d:
            continue
        counts["nonsquare_kernel_classes"] += 1
        even_extensions = {v for v in BITS[1:]
                           if all(x % 2 == 0 for x in matvec(G, v))
                           and dot(v, matvec(G, v)) % 8 == 0}
        assert len(even_extensions) <= 2
        swaps = {(v, action(g, v)) for g in group for v in BITS[1:]
                 if action(g, action(g, v)) == v}
        for g in group:
            for v in BITS[1:]:
                w = action(g, v)
                if w == v:
                    continue
                t = tuple((x + y) % 2 for x, y in zip(v, w))
                if t in even_extensions:
                    counts["maps_with_even_third_lattice"] += 1
                    assert (v, w) in swaps
    return dict(counts)


def plane_transform_checks(doc):
    counts = Counter()
    J8 = [[0] * 8 for _ in range(8)]
    for i in range(4):
        J8[i][i + 4], J8[i + 4][i] = 1, -1
    for row in doc["lattices"]:
        G = row["G"]
        d = determinant(G)
        assert d % 4 == 1
        H = integral([[d * x for x in r] for r in inverse(G)])
        assert all(H[i][i] % 2 == 0 for i in range(4))
        delta = (3 * d + 1) // 2
        assert delta % 2 == 0 and 2 * delta - 3 * d == 1
        Gamma = [[2 * x for x in I[i]] + [3 * x for x in G[i]]
                 for i in range(4)]
        Gamma += [H[i] + [delta * x for x in I[i]] for i in range(4)]
        assert multiply(multiply(Gamma, J8), transpose(Gamma)) == J8
        counts["lattices"] += 1
        for a in BITS[1:]:
            for c in BITS[1:]:
                b = action(G, c)
                if dot(a, b) % 2:
                    continue
                anew = tuple((delta * x - y) % 2
                             for x, y in zip(a, matvec(H, b)))
                bnew = tuple((-3 * x + 2 * y) % 2
                             for x, y in zip(matvec(G, a), b))
                assert anew == c and bnew == action(G, a)
                assert quadratic_value(H, b) == quadratic_value(G, c)
                counts["characteristic_swaps_and_dual_Q_checks"] += 1
    return dict(counts)


def decomposable_control():
    G = [[2, 1, 0, 0], [1, 2, 0, 0],
         [0, 0, 2, 1], [0, 0, 1, 4]]
    a, c = (1, 0, 1, 0), (0, 1, 0, 1)
    b = action(G, c)
    T = [[(-1 if i < 2 else 1) * int(i == j)
          for j in range(4)] for i in range(4)]
    plane_Q = [quadratic_value(G, v) for v in
               (a, c, tuple(x + y for x, y in zip(a, c)))]
    assert determinant(G) == 21 and plane_Q == [0, 1, 1]
    assert dot(a, b) % 2 == 0
    witness = symmetry_data(G, a, b, T)
    assert witness["epsilon"] == 1
    return {"G": G, "a": a, "b": b, "c": c, "determinant": 21,
            "Q_on_nonzero_plane_vectors": plane_Q,
            "witness": T, "witness_epsilon": witness["epsilon"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("scratch/rank4-stage2/strategy-checks.json"))
    args = parser.parse_args()
    sources = {"odd_det": Path("reports/rank4-parity/task-b-det201.json"),
               "kernels": Path("reports/rank4-parity/task-c-det64.json")}
    data = {name: json.loads(_source_bytes(path))
            for name, path in sources.items()}
    control = decomposable_control()
    control_kernel = gram_of(control["G"], kernel_basis(control["b"]))
    kernel_data = {"lattices": data["kernels"]["lattices"] + [
        {"G": control_kernel, "automorphisms": isometries(control_kernel, control_kernel)}]}
    result = {
        "status": "all supplementary checks passed; conjecture remains open",
        "sources": {name: {"path": str(path),
                            "sha256": hashlib.sha256(_source_bytes(path)).hexdigest()}
                    for name, path in sources.items()},
        "plane_transform": plane_transform_checks(data["odd_det"]),
        "mixed_kernel_and_mapping_universe": "32 stored kernels plus the computed det-84 control kernel",
        "mixed_kernel_and_mapping": mixed_and_mapping_checks(kernel_data),
        "decomposable_control": control,
    }
    target = args.output
    if target.exists():
        raise FileExistsError(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes((json.dumps(result, indent=2) + "\n").encode("utf-8"))
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("sources", "decomposable_control")}, indent=2))


if __name__ == "__main__":
    main()
