"""Verify every even E7/E8 class using batch shells and exact symmetry proofs.

The ambient backend uses different shell enumeration and direct sign summation,
sharing only certificate validation with the primary LDL/Walsh backend.
"""

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from time import perf_counter

from exact_theta import EnumerationLimit, InvalidCertificate, require
from exceptional_batch import (histogram_hash, ldl_histograms, structural_checks,
                                validate_pack, walsh)


def verify_classification(pack, backend="ldl", max_nodes=1_000_000):
    total_start = start = perf_counter()
    G, rank, bound, determinant = validate_pack(pack)
    proof_seconds = perf_counter() - start
    start = perf_counter()
    if backend == "ldl":
        histograms, stats = ldl_histograms(G, bound, max_nodes)
        coefficients = {key: walsh(histogram, rank) for key, histogram in histograms.items()}
        algorithm = "exact-rational-ldl-v1 + integer Walsh transform"
    elif backend == "ambient":
        from exceptional_ambient import ambient_histograms, direct_coefficients, gram
        require(gram(rank) == G, "ambient Gram does not match the certificate")
        histograms, stats = ambient_histograms(rank, bound)
        coefficients = direct_coefficients(histograms, rank)
        algorithm = "integer Euclidean sum-of-squares + direct sign summation"
    else:
        raise InvalidCertificate("unknown backend")
    shell_seconds = perf_counter() - start
    zero, nonzero, witness_counts = 0, 0, Counter()
    minimum_norm4 = {}
    for am, N in histograms:
        minimum_norm4[am] = min(N, minimum_norm4.get(am, N))
    minimal_shell_detects_nonvanishing = True
    for entry in pack["classes"]:
        am, bm, evidence = entry["a_mask"], entry["b_mask"], entry["evidence"]
        if evidence["kind"] == "symmetry":
            # This is a consistency check. The exact seed and transport proofs
            # checked above are what prove identical vanishing.
            require(all(values[bm] == 0 for (aa, _), values in coefficients.items() if aa == am),
                    "symmetry certificate contradicts a complete shell")
            zero += 1
        else:
            N = evidence["norm4"]
            actual = coefficients.get((am, N), [0] * (1 << rank))[bm]
            require(actual == evidence["coefficient"] and actual != 0, "nonzero coefficient mismatch")
            minimal_shell_detects_nonvanishing &= coefficients[(am, minimum_norm4[am])][bm] != 0
            nonzero += 1
            witness_counts[N] += 1
    for (am, _), values in coefficients.items():
        require(all(c == 0 for bm, c in enumerate(values) if (am & bm).bit_count() % 2),
                "odd characteristic has a nonzero coefficient")
    structure = structural_checks(pack)
    raw = json.dumps(G, separators=(",", ":")).encode("ascii")
    even = zero + nonzero
    return {"schema": "work7-exceptional-result-v1", "verdict": "classification_verified",
            "name": pack.get("name", f"E{rank}"), "rank": rank, "determinant": determinant,
            "gram_sha256": hashlib.sha256(raw).hexdigest(), "backend": backend, "algorithm": algorithm,
            "coverage_complete": True, "total_characteristics": 4 ** rank,
            "even_characteristics": even, "odd_characteristics_vanish_by_parity": 4 ** rank - even,
            "even_vanishing_by_symmetry": zero, "even_nonvanishing_by_complete_shell": nonzero,
            "unresolved": 0, "seed_witnesses": len(pack["seeds"]), "transport_nodes": len(pack["zero_nodes"]),
            "full_automorphism_group_claimed": False,
            "witness_norm4_counts": {str(N): witness_counts[N] for N in sorted(witness_counts)},
            "maximum_nonzero_witness_norm": str(Fraction(max(witness_counts), 4)),
            "minimal_shell_detects_all_even_nonvanishing": minimal_shell_detects_nonvanishing,
            "shell_histogram_sha256": histogram_hash(histograms), "enumeration": stats,
            "structural_checks": structure,
            "timings_seconds": {"symmetry_and_coverage": round(proof_seconds, 4),
                                "shells_and_coefficients": round(shell_seconds, 4),
                                "total": round(perf_counter() - total_start, 4)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pack", type=Path)
    parser.add_argument("--backend", choices=("ldl", "ambient"), default="ldl")
    parser.add_argument("--max-nodes", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        raw = args.pack.read_bytes()
        result = verify_classification(json.loads(raw), args.backend, args.max_nodes)
        result["pack_sha256"] = hashlib.sha256(raw).hexdigest()
    except EnumerationLimit as error:
        result = {"verdict": "unresolved", "reason": str(error)}
    except (InvalidCertificate, OSError, json.JSONDecodeError) as error:
        result = {"verdict": "invalid_certificate", "reason": str(error)}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result))
    return {"classification_verified": 0, "invalid_certificate": 1, "unresolved": 2}[result["verdict"]]


if __name__ == "__main__":
    raise SystemExit(main())
