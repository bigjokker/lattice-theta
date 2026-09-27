"""Exact Bn audit and F4/G2 certificate packs; see ROOT-LATTICES.md."""

from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from dn_exact import characteristic_tables, expand_record, gram as dn_gram, header, zero_recipe
from exact_theta import (complete_shells, identity, inverse, matvec, multiply,
                         require, symmetry_data, transpose, verify)
from exceptional_ambient import doubled_roots
from exceptional_batch import bits, export_certificate, walsh

ROOT = Path(__file__).resolve().parent.parent
F4_G = [[2, 1, 1, 1], [1, 2, 0, 0], [1, 0, 2, 0], [1, 0, 0, 2]]
F4_TO_D4 = [[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, -1, 0], [0, 0, 0, -1]]
G2_G = [[2, -3], [-3, 6]]
A2_IN_G2 = [[1, 1], [0, 1]]


def even_pairs(rank):
    for am in range(1 << rank):
        for bm in range(1 << rank):
            if (am & bm).bit_count() % 2 == 0:
                yield am, bm


def bn_formula(rank):
    require(type(rank) is int and rank >= 1, "positive integer rank required")
    return 2 ** (2 * rank - 1) + 2 ** (rank - 1) - 3 ** rank


def bn_audit(rank):
    """Complete minimal shells by LDL; zeros have explicit coordinate flips."""
    G = [[2 * int(i == j) for j in range(rank)] for i in range(rank)]
    digest = hashlib.sha256()
    count = zero = nodes = vectors = 0
    for am in range(1 << rank):
        a = bits(am, rank)
        N = 2 * am.bit_count()
        histogram = Counter()

        def observe(y, norm4):
            require(norm4 == N, "unexpected subminimal vector")
            xm = sum((((v - ai) // 2) % 2) << i for i, (v, ai) in enumerate(zip(y, a)))
            histogram[xm] += 1

        enumeration = complete_shells(G, a, [0] * rank, N, vector_observer=observe)
        nodes += enumeration["visited_nodes"]
        vectors += sum(histogram.values())
        require(sum(histogram.values()) == 2 ** am.bit_count(), "minimal shell size mismatch")
        coefficients = walsh(histogram, rank)
        for bm, coefficient in enumerate(coefficients):
            common = am & bm
            expected = 0 if common else 2 ** am.bit_count()
            require(coefficient == expected, "product criterion disagrees with complete shell")
            if common.bit_count() % 2:
                continue
            count += 1
            if common:
                i = (common & -common).bit_length() - 1
                # T flips i. The primal/dual shifts have respective basis
                # coordinates -a_i e_i, -b_i e_i; epsilon=-a_i*b_i.
                require(a[i] == bits(bm, rank)[i] == 1, "bad flip witness")
                row = [am, bm, 1, i]
                zero += 1
            else:
                row = [am, bm, 0, N, coefficient]
            digest.update((json.dumps(row, separators=(",", ":")) + "\n").encode("ascii"))
    require(zero == bn_formula(rank) and count - zero == 3 ** rank, "Bn count mismatch")
    return {"family": "B", "rank": rank, "G": G, "coverage_complete": True,
            "even_pairs": count, "even_zero": zero, "even_nonzero": count - zero,
            "odd_zero_by_parity": 4 ** rank - count, "unresolved": 0,
            "minimal_shell_vectors": vectors, "ldl_nodes": nodes,
            "records_sha256": digest.hexdigest(), "zero_proof": "coordinate flip with epsilon one"}


def build_root_pack(family):
    require(family in ("F4", "G2"), "unknown family")
    G = F4_G if family == "F4" else G2_G
    rank = len(G)
    entries = []
    U = F4_TO_D4
    Ui = [[int(v) for v in row] for row in inverse(U)]
    tables = characteristic_tables(4)
    for am, bm in even_pairs(rank):
        a, b = bits(am, rank), bits(bm, rank)
        recipe = None
        if family == "F4":
            da, db = matvec(U, a), matvec(transpose(Ui), b)
            dam = sum((v % 2) << i for i, v in enumerate(da))
            dbm = sum((v % 2) << i for i, v in enumerate(db))
            recipe = zero_recipe(dam, dbm, 4, tables)
        if recipe:
            source = expand_record(header(4), [dam, dbm, *recipe])
            T = multiply(multiply(Ui, source["evidence"]["T"]), U)
            evidence = {"kind": "symmetry", "T": T, "order": 2}
        else:
            shell = next(s for s in complete_shells(G, a, b, 4)["shells"] if s["signed_coefficient"])
            evidence = {"kind": "shell", "norm4": shell["norm4"], "coefficient": shell["signed_coefficient"]}
        entries.append({"a": a, "b": b, "evidence": evidence})
    return {"schema": "work7-root-pack-v1", "family": family, "G": G,
            "representatives": "binary-primal-and-dual", "certificates": entries}


def verify_root_pack(pack):
    require(isinstance(pack, dict) and pack.get("schema") == "work7-root-pack-v1", "unknown pack")
    family = pack.get("family")
    require(family in ("F4", "G2"), "unknown family")
    G = F4_G if family == "F4" else G2_G
    require(pack.get("G") == G and all(type(v) is int for row in pack["G"] for v in row), "wrong Gram matrix")
    require(pack.get("representatives") == "binary-primal-and-dual", "wrong representatives")
    expected = list(even_pairs(len(G)))
    entries = pack.get("certificates")
    require(isinstance(entries, list) and len(entries) == len(expected), "incomplete coverage")
    counts = Counter()
    for entry, (am, bm) in zip(entries, expected):
        require(isinstance(entry, dict), "entry must be an object")
        require(entry.get("a") == bits(am, len(G)) and entry.get("b") == bits(bm, len(G)), "wrong class order")
        document = {**entry, "schema": "work7-theta-v1", "G": G}
        result = verify(document)
        require(result["verdict"] != "unresolved", "unresolved entry")
        counts[result["verdict"]] += 1
    return {"family": family, "rank": len(G), "coverage_complete": True,
            "even_pairs": len(expected), "even_zero": counts["vanishes_by_symmetry"],
            "even_nonzero": counts["does_not_vanish"],
            "odd_zero_by_parity": 4 ** len(G) - len(expected), "unresolved": 0}


def e8_scope_example(pack):
    """Export a verified quarter-coordinate xi illustrating the missing sector."""
    C = doubled_roots(8)
    A2, B2 = [5, 1, 1, 1, 1, 1, 1, 1], [2, 2, 2, 2, 0, 0, 0, 0]
    raw_a = matvec(inverse(C), A2)
    raw_b = [Fraction(v, 4) for v in matvec(transpose(C), B2)]
    require(all(v.denominator == 1 for v in raw_a), "nonintegral primal coordinates")
    require(all(v.denominator == 1 for v in raw_b), "nonintegral dual coordinates")
    a, b = [int(v) for v in raw_a], [int(v) for v in raw_b]
    am = sum((v % 2) << i for i, v in enumerate(a))
    bm = sum((v % 2) << i for i, v in enumerate(b))
    entry = next(e for e in pack["classes"] if e["a_mask"] == am and e["b_mask"] == bm)
    document = export_certificate(pack, entry)
    document.update(name="E8 quarter-coordinate xi: 2xi=s+2e1", a=a, b=b)
    require(verify(document)["verdict"] == "vanishes_by_symmetry", "example not certified")
    return document


def build():
    require(multiply(multiply(transpose(F4_TO_D4), dn_gram(4)), F4_TO_D4) == F4_G, "F4 isometry failed")
    require(multiply(multiply(transpose(A2_IN_G2), G2_G), A2_IN_G2) == [[2, -1], [-1, 2]], "G2 isometry failed")
    path = ROOT / "data" / "root-lattices"
    path.mkdir(parents=True, exist_ok=True)
    results = [bn_audit(n) for n in range(1, 9)]
    for family in ("F4", "G2"):
        pack = build_root_pack(family)
        result = verify_root_pack(pack)
        require(result["even_zero"] == (9 if family == "F4" else 0), "family count mismatch")
        (path / f"{family.lower()}-classes.json").write_text(json.dumps(pack, indent=2) + "\n", encoding="utf-8")
        results.append(result)
    for rank in (2, 8):
        a = [1, 1] + [0] * (rank - 2)
        T = identity(rank)
        T[0][0] = -1
        doc = {"schema": "work7-theta-v1", "name": f"B{rank} even zero", "G": [[2 * int(i == j) for j in range(rank)] for i in range(rank)],
               "a": a, "b": a, "evidence": {"kind": "symmetry", "T": T, "order": 2}}
        verify(doc)
        (path / f"b{rank}-zero.json").write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    pack = json.loads((ROOT / "data" / "e8-classification.json").read_text(encoding="utf-8"))
    doc = e8_scope_example(pack)
    (path / "e8-quarter-xi-zero.json").write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    (ROOT / "reports" / "root-lattices-audit.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(results))


if __name__ == "__main__":
    build()
