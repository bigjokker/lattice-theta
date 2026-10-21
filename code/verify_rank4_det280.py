"""Independent replay of the determinant-280 level-four zero.

The lattice G is half the Gram matrix of the determinant-4480 level-four
certificate. Scaling a Gram matrix by a positive rational changes neither
identical vanishing (Theta_{G/2}[a,b](tau) = Theta_G[a,b](tau/2)), nor the
sign kernel and its two coset lattices, nor isometry between them, so the
determinant-4480 certificate already implies the result. This script checks
that certificate with the unmodified level-four verifier, checks the scaling,
and then re-proves every claim directly for the halved Gram:

- the full automorphism group (complete box search) and the signs of its
  stabilizer, all zero, so the zero has no symmetry witness;
- the zero itself, by the modular cutoff and complete signed shells;
- level four, by a complete isometry search between the two coset lattices.

Only verifier modules are imported; nothing comes from a builder.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

from exact_theta import dot, matvec, require, symmetry_data, transpose, inverse
from verify_rank4_feasibility import box_isometries, det4, modular_check
from verify_rank4_level4 import box_signed_shells, gram_of, replay_certificate

ROOT = Path(__file__).resolve().parent.parent


def load_source(derived):
    path = ROOT / derived["file"]
    raw = path.read_bytes() if path.exists() else gzip.decompress((ROOT / (derived["file"] + ".gz")).read_bytes())
    require(hashlib.sha256(raw).hexdigest() == derived["sha256"], "source certificate hash")
    doc = json.loads(raw)
    found = [c for c in doc["certificates"] if c["determinant"] == derived["certificate_determinant"]]
    require(len(found) == 1 and found[0]["pair"] == derived["pair"], "source certificate")
    return found[0]


def replay(doc, max_nodes, enumeration_nodes):
    require(doc["schema"] == "work7-rank4-level4-rescaled-certificate-v1", "schema")
    derived = doc["derived_from"]
    source = load_source(derived)
    source_result = replay_certificate(source, max_nodes, enumeration_nodes)
    require(source_result["verdict"] == "verified_level_4", "source certificate does not verify")

    k = derived["gram_divisor"]
    require(k == 2 and all(x % k == 0 for row in source["G"] for x in row), "Gram not divisible")
    G = [[x // k for x in row] for row in source["G"]]
    require(G == doc["G"] and doc["a"] == source["a"] and doc["b"] == source["b"], "rescaled data")
    a, b = doc["a"], doc["b"]
    det = det4(G)
    require(det == doc["determinant"] == 280 and det * 16 == source["determinant"], "determinant")
    require(any(G[i][i] % 2 for i in range(4)), "expected an odd lattice")
    require(any(a) and any(b) and dot(a, b) % 2 == 0, "even nontrivial characteristic")

    # Full group and stabilizer signs (no symmetry witness).
    group = box_isometries(G, G, max_nodes)
    require(group == source["automorphisms"], "group differs from the source certificate")
    exp = doc["expected"]
    require(len(group) == exp["full_group_order"], "group order")
    stabilizer = []
    for T in group:
        if any((x - y) % 2 for x, y in zip(matvec(T, a), a)):
            continue
        if any((x - y) % 2 for x, y in zip(matvec(transpose(inverse(T)), b), b)):
            continue
        stabilizer.append(symmetry_data(G, a, b, T)["epsilon"])
    require(stabilizer and all(e == 0 for e in stabilizer), "an automorphism is a witness")

    # The zero: modular cutoff, then complete signed shells to the cutoff.
    B = modular_check(G, a, b, doc["modular"])
    require(B == exp["cutoff"], "cutoff")
    shells, nodes = box_signed_shells(G, a, b, B, enumeration_nodes)
    require(all(s["signed_coefficient"] == 0 for s in shells), "a signed coefficient is nonzero")
    require(len(shells) == exp["shells"] and sum(s["vector_count"] for s in shells) == exp["vectors"]
            and hashlib.sha256(json.dumps(shells).encode()).hexdigest() == exp["shells_sha256"], "shell record")

    # Level four: the two coset lattices (bases from the source certificate,
    # whose structure the source replay checked) are not isometric.
    reduced = {}
    for name in ("O_v1", "O_v2"):
        entry = source["coset_lattices"][name]
        H = gram_of(G, entry["doubled_basis"])  # 4 x the coset Gram for G
        require(H == [[x // k for x in row] for row in gram_of(source["G"], entry["doubled_basis"])], "coset scaling")
        reduced[name] = gram_of(H, entry["reduced"]["basis"])
        require(abs(det4(entry["reduced"]["basis"])) == 1, "coset reduction")
    require(box_isometries(reduced["O_v1"], reduced["O_v2"], max_nodes) == [], "coset lattices are isometric")

    return {"verdict": "verified_level_4", "determinant": det, "G": G, "a": a, "b": b,
            "source": {"pair": derived["pair"], "determinant": source["determinant"],
                       "verdict": source_result["verdict"]},
            "full_group_order": len(group), "stabilizer_epsilons": stabilizer,
            "cutoff": B, "shells": len(shells), "vectors": exp["vectors"], "box_nodes": nodes,
            "coset_lattices_isometric": False, "level": 4}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", type=Path, default=ROOT / "data/rank4-level4/det280-level4.json")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-nodes", type=int, default=5_000_000)
    parser.add_argument("--enumeration-nodes", type=int, default=100_000_000)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    doc = json.loads(args.certificate.read_text(encoding="utf-8"))
    result = replay(doc, args.max_nodes, args.enumeration_nodes)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k not in ("G",)}, indent=2))


if __name__ == "__main__":
    main()
