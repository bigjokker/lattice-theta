"""Independent E6 cross-check by Euclidean sum-of-squares enumeration.

Standard library only. Does NOT import exact_theta, the builder, or the pack
verifier. Uses the explicit doubled roots in R8 and integer square bounds.
See E6-MILESTONE.md for the inverse-coordinate and completeness proofs.
"""

import argparse
from collections import Counter, defaultdict
import hashlib
from itertools import product
import json
from math import isqrt
from pathlib import Path

# Columns are 2 alpha_i, ordered as Bourbaki roots (1,3,4,5,6,2).
DOUBLED_ROOTS = [
    [1, -2, 0, 0, 0, 2],
    [-1, 2, -2, 0, 0, 2],
    [-1, 0, 2, -2, 0, 0],
    [-1, 0, 0, 2, -2, 0],
    [-1, 0, 0, 0, 2, 0],
    [-1, 0, 0, 0, 0, 0],
    [-1, 0, 0, 0, 0, 0],
    [1, 0, 0, 0, 0, 0],
]


def check(condition, message):
    if not condition:
        raise ValueError(message)


def gram_from_roots():
    raw = [[sum(row[i] * row[j] for row in DOUBLED_ROOTS) for j in range(6)] for i in range(6)]
    check(all(x % 4 == 0 for row in raw for x in row), "root Gram not integral")
    return [[x // 4 for x in row] for row in raw]


def recover_coordinates(z, h):
    """Recover integral y from (z0,...,z4,-h,-h,h)=C y, or reject."""
    y = [h, 0, 0, 0, 0, 0]
    for j in (4, 3, 2):
        numerator = z[j] + h + (2 * y[j + 1] if j < 4 else 0)
        if numerator % 2:
            return None
        y[j] = numerator // 2
    numerator = z[0] + z[1] + 2 * y[2]
    if numerator % 4:
        return None
    y[5] = numerator // 4
    numerator = h + 2 * y[5] - z[0]
    if numerator % 2:
        return None
    y[1] = numerator // 2
    ambient = z + [-h, -h, h]
    check([sum(c * v for c, v in zip(row, y)) for row in DOUBLED_ROOTS] == ambient,
          "inverse-coordinate error")
    return y


def ambient_shells(norm4_bound=4):
    check(type(norm4_bound) is int and norm4_bound >= 0, "invalid enumeration bound")
    histograms = defaultdict(Counter)
    ambient_leaf_count = accepted = 0
    seen = set()
    # ||C y||^2 = 4 y^t G y = z0^2+...+z4^2+3h^2.
    budget = 4 * norm4_bound
    for h in range(-isqrt(budget // 3), isqrt(budget // 3) + 1):
        prefix = []

        def visit(remaining):
            nonlocal ambient_leaf_count, accepted
            if len(prefix) == 5:
                ambient_leaf_count += 1
                y = recover_coordinates(prefix, h)
                if y is None:
                    return
                check(tuple(y) not in seen, "duplicate lattice vector")
                seen.add(tuple(y))
                norm_scaled = sum(v * v for v in prefix) + 3 * h * h
                check(norm_scaled % 4 == 0, "nonintegral norm")
                N = norm_scaled // 4
                a = tuple(v % 2 for v in y)
                x_mask = sum((((v - ai) // 2) % 2) << j for j, (v, ai) in enumerate(zip(y, a)))
                histograms[(a, N)][x_mask] += 1
                accepted += 1
                return
            radius = isqrt(remaining)
            lo = -radius + ((h + radius) % 2)
            for z in range(lo, radius + 1, 2):
                prefix.append(z)
                visit(remaining - z * z)
                prefix.pop()

        visit(budget - 3 * h * h)
    return histograms, {"enumeration_bound_norm4": norm4_bound,
                        "ambient_leaves": ambient_leaf_count, "accepted_lattice_vectors": accepted}


def signed_sum(histogram, b):
    mask = sum(x << i for i, x in enumerate(b))
    return sum(count * (-1 if (xmask & mask).bit_count() % 2 else 1)
               for xmask, count in histogram.items())


def verify_ambient(pack):
    check(isinstance(pack, dict) and pack.get("schema") == "work7-even-theta-pack-v1", "unknown pack")
    check(pack.get("representatives") == "binary-primal-and-dual", "unsupported representatives")
    G = pack.get("G")
    check(isinstance(G, list) and len(G) == 6 and all(isinstance(row, list) and len(row) == 6 for row in G),
          "Gram must be 6 by 6")
    check(all(type(x) is int for row in G for x in row), "Gram entries must be integers")
    check(G == gram_from_roots(), "input Gram is not the specified E6 realization")
    representatives = list(product(range(2), repeat=6))
    expected = {(a, b) for a in representatives for b in representatives
                if sum(ai * bi for ai, bi in zip(a, b)) % 2 == 0}
    entries = pack.get("certificates")
    check(isinstance(entries, list), "certificates must be a list")
    actual = set()
    # Validate all coverage and evidence before doing the enumeration.
    for entry in entries:
        check(isinstance(entry, dict), "entry must be an object")
        for label in ("a", "b"):
            v = entry.get(label)
            check(isinstance(v, list) and len(v) == 6 and all(type(x) is int and x in (0, 1) for x in v),
                  "characteristic representatives must be six binary integers")
        key = (tuple(entry["a"]), tuple(entry["b"]))
        check(key in expected, "odd class in even pack")
        check(key not in actual, "duplicate characteristic")
        actual.add(key)
        evidence = entry.get("evidence")
        check(isinstance(evidence, dict) and evidence.get("kind") == "shell", "shell evidence required")
        check(type(evidence.get("norm4")) is int and 0 <= evidence["norm4"] <= 4,
              "ambient cross-check covers norm4 from 0 through 4")
        check(type(evidence.get("coefficient")) is int and evidence["coefficient"] != 0,
              "nonzero integral coefficient required")
    check(actual == expected, f"coverage incomplete: expected {len(expected)}, found {len(actual)}")
    histograms, stats = ambient_shells(4)
    counts = Counter()
    for entry in entries:
        evidence = entry["evidence"]
        histogram = histograms.get((tuple(entry["a"]), evidence["norm4"]), {})
        c = signed_sum(histogram, entry["b"])
        check(c == evidence["coefficient"], f"coefficient mismatch for a={entry['a']}, b={entry['b']}")
        check(c != 0, "nonvanishing not certified")
        counts[evidence["norm4"]] += 1
    # This also checks the stronger minimal-shell assertion without relying on the table.
    minimal_counts = Counter()
    for a in representatives:
        shells = sorted(N for aa, N in histograms if aa == a)
        check(bool(shells), "coset missing from the complete sphere")
        N = shells[0]
        vector_count = sum(histograms[(a, N)].values())
        if any(a):
            check(vector_count % 4 == 2, "minimal antipodal-pair count must be odd")
        else:
            check(N == 0 and vector_count == 1, "zero coset must have unique zero vector")
        minimal_counts[(N, vector_count)] += 1
    odd_checked = 0
    for a in representatives:
        for b in representatives:
            if sum(ai * bi for ai, bi in zip(a, b)) % 2 == 0:
                continue
            for N in range(5):
                check(signed_sum(histograms.get((a, N), {}), b) == 0,
                      "odd characteristic has a nonzero coefficient")
            odd_checked += 1
    return {"schema": "work7-e6-ambient-result-v1", "verdict": "parity_converse_verified",
            "algorithm": "integer-euclidean-sum-of-squares-v1; no shared enumeration code",
            "coverage_complete": True, "even_characteristics_verified": len(actual),
            "odd_characteristics_checked_through_bound": odd_checked,
            "witness_norm4_counts": {str(k): counts[k] for k in sorted(counts)},
            "minimal_coset_shell_counts": [{"norm4": N, "vector_count": count, "cosets": multiplicity}
                                           for (N, count), multiplicity in sorted(minimal_counts.items())],
            **stats}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pack", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        raw = args.pack.read_bytes()
        result = verify_ambient(json.loads(raw))
        result["pack_sha256"] = hashlib.sha256(raw).hexdigest()
    except (ValueError, OSError, json.JSONDecodeError) as error:
        result = {"verdict": "invalid_certificate", "reason": str(error)}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result))
    return 0 if result["verdict"] == "parity_converse_verified" else 1


if __name__ == "__main__":
    raise SystemExit(main())
