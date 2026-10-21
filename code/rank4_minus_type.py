"""Checks for the new rank-four odd-determinant symmetry theorem.

docs/rank-four-minus-type.md contains the proof; finite checks are supplementary.
No search below is a proof that identity zeros without witnesses are absent.

Release layout (run from the repository root): inputs are read from data/ and
reports/ (a preserved .json.gz is read decompressed); write fresh outputs under
scratch/. This is the research-workspace script of 2026-10-21 with relocated
input paths; the preserved outputs under reports/ were produced by the
workspace copy, so paths recorded inside them are the original workspace paths.
"""

import argparse
from collections import Counter
import hashlib
import gzip
from itertools import product
import json
from pathlib import Path
import time

from exact_theta import identity, inverse, matvec, multiply, transpose
from rank4_feasibility import determinant

I = identity(4)
NEG_I = [[-x for x in row] for row in I]
BITS = list(product((0, 1), repeat=4))

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


def power(T, n):
    result = I
    while n:
        if n & 1:
            result = multiply(result, T)
        T = multiply(T, T)
        n //= 2
    return result


def order(T, modulo=None):
    current = I
    for n in range(1, 65):
        current = multiply(current, T)
        if modulo:
            current = [[x % modulo for x in row] for row in current]
        if current == I:
            return n
    raise ValueError("matrix order exceeds 64")


def action(T, x):
    return tuple(v % 2 for v in matvec(T, x))


def dot(x, y):
    return sum(a*b for a, b in zip(x, y))


def symplectic_form(x, y):
    return (x[0]*y[1] + x[1]*y[0] + x[2]*y[3] + x[3]*y[2]) % 2


def quadratic(x, minus=False):
    return (x[0]*x[1] + x[2]*x[3] + (x[0]+x[1] if minus else 0)) % 2


def finite_group_check():
    """Exhaust all 15^4 nonzero column choices, testing the defining form."""
    basis = transpose(I)
    target = [[symplectic_form(x, y) for y in basis] for x in basis]
    tables = {name: Counter() for name in ("Sp4_F2", "O4_plus_F2", "O4_minus_F2")}
    examples = {}
    for cols in product(BITS[1:], repeat=4):
        if any(symplectic_form(cols[i], cols[j]) != target[i][j]
               for i in range(4) for j in range(i)):
            continue
        T = transpose(cols)
        images = [action(T, x) for x in BITS]
        assert len(set(images)) == 16
        n = order(T, 2)
        fixed = sum(x == y for x, y in zip(BITS, images))
        key = f"order={n},fixed_vectors={fixed}"
        tables["Sp4_F2"][key] += 1
        if n == 4:
            assert fixed == 2, "order-four symplectic map must fix one line"
        examples.setdefault(key, T)
        for minus, name in ((False, "O4_plus_F2"), (True, "O4_minus_F2")):
            if all(quadratic(x, minus) == quadratic(y, minus) for x, y in zip(BITS, images)):
                tables[name][key] += 1
    assert [sum(t.values()) for t in tables.values()] == [720, 72, 120]
    return {"column_choices": 15**4, "group_orders": {k: sum(v.values()) for k, v in tables.items()},
            "tables": {k: dict(sorted(v.items())) for k, v in tables.items()},
            "symplectic_representatives": examples,
            "order_four_fixed_line_lemma": "verified on all 720 symplectic maps"}


def audit_lattice(G, group):
    """Check all stabilizer signs directly, with a separate dual-sign check."""
    assert all(G[i][i] % 2 == 0 for i in range(4)) and determinant(G) % 2 == 1
    split = any(T not in (I, NEG_I) and multiply(T, T) == I
                and all((T[i][j] + I[i][j]) % 2 == 0 for i in range(4) for j in range(4))
                for T in group)
    Q = lambda x: (dot(x, matvec(G, x)) // 2) % 2
    minus = sum(Q(x) == 0 for x in BITS) == 6
    actual, gaussian = set(), set()
    shapes, witnesses = Counter(), Counter()
    for T in group:
        assert multiply(multiply(transpose(T), G), T) == G
        n = order(T)
        r = order(T, 2)
        fixed = [x for x in BITS[1:] if action(T, x) == x]
        shapes[f"order={n},mod2_order={r},fixed_nonzero={len(fixed)}"] += 1
        if r == 4:
            assert len(fixed) == 1
        if n == 8:
            assert power(T, 4) == NEG_I and r == 4
        root = multiply(T, T) == NEG_I
        if root:
            assert len(fixed) == 3 and all(Q(x) == 0 for x in fixed) and not minus
        dual = transpose(inverse(T))
        for a, c in product(fixed, repeat=2):
            b = tuple(x % 2 for x in matvec(G, c))
            if dot(a, b) % 2:
                continue
            shifted = [(v - x) // 2 for v, x in zip(matvec(T, a), a)]
            eps = dot(shifted, b) % 2
            dual_shift = [(int(v) - x) // 2 for v, x in zip(matvec(dual, b), b)]
            assert eps == dot(a, dual_shift) % 2
            if root and a != c:
                gaussian.add((a, c))
                assert eps == 1
            if eps:
                actual.add((a, c))
                witnesses[str(n)] += 1
                if not split:
                    odd_part = n
                    while odd_part % 2 == 0:
                        odd_part //= 2
                    h = power(T, odd_part)
                    assert multiply(h, h) == NEG_I
    if not split:
        assert actual == gaussian
        if minus:
            assert not actual
    return {"determinant": determinant(G), "minus_type": minus, "decomposable": split,
            "automorphisms": len(group), "symmetry_zero_count": len(actual),
            "gaussian_zero_count": len(gaussian), "witness_orders_with_multiplicity": dict(witnesses),
            "automorphism_shapes": dict(sorted(shapes.items()))}


def audit_census(path):
    raw = _source_bytes(path)
    doc = json.loads(gzip.decompress(raw) if str(path).endswith(".gz") else raw)
    rows = []
    for row in doc["lattices"]:
        result = audit_lattice(row["G"], row["automorphisms"])
        assert result["decomposable"] == row["decomposable"]
        result["source_index"] = row["index"]
        rows.append(result)
    indec = [r for r in rows if not r["decomposable"]]
    return {"source": str(path), "sha256": hashlib.sha256(raw).hexdigest(),
            "classes": len(rows), "indecomposable": len(indec),
            "minus_indecomposable": sum(r["minus_type"] for r in indec),
            "minus_indecomposable_symmetry_zeros": sum(r["symmetry_zero_count"] for r in indec if r["minus_type"]),
            "indecomposable_symmetry_zeros": sum(r["symmetry_zero_count"] for r in indec),
            "symmetry_theorem_matches_every_indecomposable_class": True, "lattices": rows}


def mixed_kernel_check(path="reports/rank4-parity/task-c-det64.json"):
    """Test the mixed order-four kernel lemma on complete saved groups."""
    raw = _source_bytes(path)
    doc = json.loads(raw)
    rows, total_elements, total_orbit_steps = [], 0, 0
    for row in doc["lattices"]:
        if row["determinant"] % 8 != 4:
            continue  # v_2(det M) = 2, exactly the lemma's hypothesis.
        G, group = row["G"], row["automorphisms"]
        involutions = [R for R in group if multiply(R, R) == I]
        elements, steps = 0, 0
        for T in group:
            if order(T) != 4 or multiply(T, T) == NEG_I:
                continue
            elements += 1
            for v in BITS[1:]:
                w = action(T, v)
                if w == v or action(T, w) == v:
                    continue
                steps += 1
                assert any(action(R, v) == w for R in involutions), (G, T, v)
        rows.append({"source_index": row["index"], "determinant": row["determinant"],
                     "mixed_order_four_elements": elements, "four_cycle_steps": steps})
        total_elements += elements
        total_orbit_steps += steps
    return {"source": str(path), "sha256": hashlib.sha256(raw).hexdigest(),
            "kernel_classes": len(rows), "mixed_order_four_elements": total_elements,
            "four_cycle_steps_with_involutive_swap": total_orbit_steps,
            "all_steps_have_an_involutive_swap": True, "lattices": rows}


def quarter_turn_probe(limit=15):
    """A bounded construction probe, not a census.

    Kernel M extends diag(2A,2A,2B,2C) by (1,1,1,1)/2.
    R rotates the first two coordinates and negates the fourth. Seek
    R-related cosets whose third overlattice is even of minus type.
    """
    from rank4_level4_search import gram_of, overlattice_basis, lll_reduce
    from rank4_feasibility import isometries
    T = [[0, -1, 0, 0], [1, 0, 0, 1], [0, 0, 1, 1], [0, 0, 0, -1]]
    assert power(T, 4) == I and power(T, 2) != I
    kernels, configurations, records, seen = 0, 0, [], set()
    for A, B, C in product(range(1, limit+1, 2), repeat=3):
        d = A*A*B*C
        if (2*A+B+C) % 4 or d % 8 != 5:
            continue
        G0 = [[2*A, 0, 0, A], [0, 2*A, 0, A], [0, 0, 2*B, B],
              [A, A, B, (2*A+B+C)//2]]
        assert determinant(G0) == 4*d
        assert multiply(multiply(transpose(T), G0), T) == G0
        kernels += 1
        for a in BITS[1:]:
            c = action(T, a)
            if c == a:
                continue
            v3 = [x ^ y for x, y in zip(a, c)]
            gram4 = gram_of(G0, overlattice_basis(v3))
            if any(x % 4 for row in gram4 for x in row) or any(gram4[i][i] % 8 for i in range(4)):
                continue
            configurations += 1
            G = [[x//4 for x in row] for row in gram4]
            _, G = lll_reduce(G)
            key = tuple(map(tuple, G))
            if key in seen:
                continue
            seen.add(key)
            group = isometries(G, G, 5_000_000)
            result = audit_lattice(G, group)
            records.append({"parameters": [A, B, C], "G0": G0, "a": a, "Ra": c,
                            "v3": v3, "G": G, **result})
    return {"universe": "odd A,B,C in [1,limit], (2A+B+C) divisible by 4, A^2 B C = 5 mod 8",
            "parameter_limit": limit, "kernel_parameter_triples": kernels,
            "configurations_with_even_third_lattice": configurations,
            "distinct_reduced_matrices_not_isometry_classes": len(records),
            "indecomposable_hits": sum(not r["decomposable"] for r in records), "lattices": records}


def extend_minus(dmin, dmax):
    """Use the existing complete generator with an explicitly narrower filter."""
    import rank4_parity_search as builder
    builder.DET_FILTERS["5_mod_8"] = lambda d: d % 8 == 5
    progress = lambda message: print(message, flush=True)
    classes = builder.class_list(dmax, even_only=True, det_filter="5_mod_8", dmin=dmin, progress=progress)
    return builder.decide_class_list(classes, progress=progress)


def replay_extension(path):
    """Independent loose-box generator, isometries, and coefficient replay."""
    import verify_rank4_parity as verifier
    verifier.DET_FILTERS["5_mod_8"] = lambda d: d % 8 == 5
    verifier.FAST_BOX[0] = True
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as source:
        doc = json.load(source)
    domain = doc["classes_doc"]
    assert domain["even_only"] and domain["determinant_filter"] == "5_mod_8"
    result = verifier.replay_decided(doc, 5_000_000, 50_000_000,
                                    progress=lambda message: print(message, flush=True))
    return {"source": str(path), "source_sha256": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            "dmin": domain["determinant_minimum"], "dmax": domain["determinant_bound"], **result}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", choices=("finite", "census", "mixed-kernels", "quarter-turn", "extend", "replay-extension"))
    parser.add_argument("--source", default="reports/rank4-parity/task-b-det201.json")
    parser.add_argument("--output", required=True)
    parser.add_argument("--limit", type=int, default=15)
    parser.add_argument("--dmin", type=int, default=202)
    parser.add_argument("--dmax", type=int, default=501)
    args = parser.parse_args()
    out = Path(args.output)
    if out.exists():
        raise FileExistsError(out)
    started = time.perf_counter()
    if args.task == "finite":
        result = finite_group_check()
    elif args.task == "census":
        result = audit_census(args.source)
    elif args.task == "mixed-kernels":
        result = mixed_kernel_check()
    elif args.task == "quarter-turn":
        result = quarter_turn_probe(args.limit)
    elif args.task == "extend":
        result = extend_minus(args.dmin, args.dmax)
    else:
        result = replay_extension(args.source)
    result["seconds"] = time.perf_counter() - started
    out.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(out).endswith(".gz") else open
    with opener(out, "wt", encoding="utf-8") as target:
        json.dump(result, target, indent=None if str(out).endswith(".gz") else 2)
        target.write("\n")
    print(json.dumps({k: v for k, v in result.items() if k not in ("lattices", "symplectic_representatives", "classes_doc")}, indent=2))


if __name__ == "__main__":
    main()
