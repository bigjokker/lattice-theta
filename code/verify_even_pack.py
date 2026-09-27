"""Check complete nonvanishing coverage of all even characteristic classes.

Each entry is passed to the existing exact_theta.verify function. Coverage is
checked separately, without automorphism orbits, using binary representatives.
"""

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import product
import json
from pathlib import Path

from exact_theta import InvalidCertificate, dot, matrix, require, vector, verify


def verify_pack(pack, max_nodes=1_000_000):
    require(isinstance(pack, dict), "pack must be an object")
    require(pack.get("schema") == "work7-even-theta-pack-v1", "unknown pack schema")
    require(pack.get("representatives") == "binary-primal-and-dual", "unsupported representatives")
    G = matrix(pack.get("G"), "G")
    n = len(G)
    representatives = list(product((0, 1), repeat=n))
    expected = {(a, b) for a in representatives for b in representatives if dot(a, b) % 2 == 0}
    entries = pack.get("certificates")
    require(isinstance(entries, list), "certificates must be a list")
    seen = set()
    for entry in entries:
        require(isinstance(entry, dict), "each certificate must be an object")
        a = vector(entry.get("a"), "a", n)
        b = vector(entry.get("b"), "b", n)
        require(all(x in (0, 1) for x in a + b), "representatives must be binary")
        key = (tuple(a), tuple(b))
        require(key in expected, "pack includes an odd characteristic")
        require(key not in seen, "duplicate characteristic")
        seen.add(key)
        evidence = entry.get("evidence")
        require(isinstance(evidence, dict) and evidence.get("kind") == "shell",
                "every class must have a nonzero-shell certificate")
    require(seen == expected, f"coverage incomplete: expected {len(expected)}, found {len(seen)}")
    results, unresolved, histogram = [], [], Counter()
    gram_hash = determinant = None
    for entry in entries:
        doc = {"schema": "work7-theta-v1", "G": G, "a": entry["a"], "b": entry["b"],
               "evidence": entry["evidence"]}
        result = verify(doc, max_nodes)
        key = {"a_mask": sum(x << i for i, x in enumerate(entry["a"])),
               "b_mask": sum(x << i for i, x in enumerate(entry["b"]))}
        gram_hash, determinant = result["gram_sha256"], result["determinant"]
        if result["verdict"] == "unresolved":
            unresolved.append({**key, "reason": result["reason"]})
            continue
        require(result["verdict"] == "does_not_vanish", "class did not prove nonvanishing")
        N = result["nonzero_norm4"]
        shell = next(s for s in result["enumeration"]["shells"] if s["norm4"] == N)
        results.append({**key, "norm4": N, "signed_coefficient": result["signed_coefficient"],
                        "vector_count": shell["vector_count"],
                        "visited_nodes": result["enumeration"]["visited_nodes"]})
        histogram[N] += 1
    return {"schema": "work7-even-theta-pack-result-v1", "name": pack.get("name", "unnamed"),
            "verdict": "unresolved" if unresolved else "parity_converse_verified",
            "rank": n, "determinant": determinant, "gram_sha256": gram_hash,
            "coverage_complete": True, "total_characteristics": 4 ** n,
            "even_characteristics": len(expected), "odd_characteristics": 4 ** n - len(expected),
            "nonvanishing_verified": len(results), "unresolved": unresolved,
            "algorithm": "exact-rational-ldl-v1; each class independently passed to exact_theta.verify",
            "witness_norm4_counts": {str(k): histogram[k] for k in sorted(histogram)},
            "maximum_witness_norm": str(Fraction(max(histogram), 4)) if histogram else None,
            "class_results": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pack", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-nodes", type=int, default=1_000_000)
    args = parser.parse_args()
    try:
        raw = args.pack.read_bytes()
        result = verify_pack(json.loads(raw), args.max_nodes)
        result["pack_sha256"] = hashlib.sha256(raw).hexdigest()
    except (InvalidCertificate, OSError, json.JSONDecodeError) as error:
        result = {"verdict": "invalid_certificate", "reason": str(error)}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "class_results"}))
    return {"parity_converse_verified": 0, "invalid_certificate": 1, "unresolved": 2}[result["verdict"]]


if __name__ == "__main__":
    raise SystemExit(main())
