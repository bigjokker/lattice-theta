"""Construction checks accompanying docs/rank-four-minus-type.md.

The Conway-Sloane template is checked as a linear identity in all four
parameters. Its parity obstruction is proved in the note, not inferred
from sampling. The determinant-280 check records a rescaling correction.

Release layout (run from the repository root): inputs are read from data/ and
reports/ (a preserved .json.gz is read decompressed); write fresh outputs under
scratch/. This is the research-workspace script of 2026-10-21 with relocated
input paths; the preserved outputs under reports/ were produced by the
workspace copy, so paths recorded inside them are the original workspace paths.
"""

import argparse
import gzip
from fractions import Fraction
import hashlib
import json
from math import gcd
from pathlib import Path

from exact_theta import multiply, symmetry_data
from rank4_feasibility import cutoff_build, determinant, isometries
from rank4_level4_search import cs_gram, gram_of
from rank4_parity_search import kernel_configuration
from verify_rank4_feasibility import box_isometries, modular_check
from verify_rank4_level4 import box_signed_shells

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


def template(a, b, c, d):
    """3 times the third-lattice Gram in a fixed basis, for arbitrary parameters."""
    return [[9*a, 3*a, 3*a, 0],
            [3*a, a+b+d, a-b, d-b],
            [3*a, a-b, a+b+c, b-c],
            [0, d-b, b-c, b+c+d]]


def template_check():
    source = Path("reports/rank4-level4/task-a.json")
    pair = json.loads(_source_bytes(source))["pairs"][0]
    cert = json.loads(_source_bytes("data/rank4-level4/det1729-level4.json"))["certificates"][0]
    doubled_basis = pair["klein_configurations"][0]["doubled_bases"]["L"]
    B = multiply([[Fraction(x, 2) for x in row] for row in doubled_basis],
                 cert["lll_transform_from_construction"])
    i1, i2, kernel_map = pair["sub_table"]["isometric_pairs"][0]
    K1 = pair["sub_family_L1"][i1]["basis"]
    K2_map = multiply(pair["sub_family_L2"][i2]["basis"], kernel_map)
    coeffs = []
    for i in range(4):
        parameters = [int(i == j) for j in range(4)]
        computed = gram_of(cs_gram(*parameters, -1), B)
        assert gram_of(cs_gram(*parameters, 1), K1) == gram_of(cs_gram(*parameters, -1), K2_map)
        wanted = template(*parameters)
        assert [[3*x for x in row] for row in computed] == wanted
        coeffs.append(wanted)
    # The assertions below are an exhaustive F2 check of the diagonal argument.
    even_diagonal_parameters = []
    for mask in range(16):
        values = [(mask >> i) & 1 for i in range(4)]
        if all(template(*values)[i][i] % 2 == 0 for i in range(4)):
            even_diagonal_parameters.append(values)
    assert even_diagonal_parameters == [[0, 0, 0, 0]]
    return {"source": str(source), "sha256": hashlib.sha256(_source_bytes(source)).hexdigest(),
            "coefficient_matrices_for_3G": coeffs,
            "linear_identity_verified": True,
            "common_kernel_isometry_verified_for_all_parameters": True,
            "even_diagonals_force_all_parameters_even_mod2": True,
            "claim": "this fixed common-kernel construction has no even lattice of odd determinant, at any rational scaling"}


def det280_check():
    source = Path("data/rank4-level4/other-level4-certificates.json")
    original = next(c for c in json.loads(_source_bytes(source))["certificates"]
                    if c["determinant"] == 4480)
    assert all(x % 2 == 0 for row in original["G"] for x in row)
    G = [[x//2 for x in row] for row in original["G"]]
    a, b = original["a"], original["b"]
    assert determinant(G) == 280 and gcd(*(x for row in G for x in row)) == 1
    group = isometries(G, G, 5_000_000)
    assert group == box_isometries(G, G, 5_000_000) and len(group) == 2
    assert all(symmetry_data(G, a, b, T)["epsilon"] == 0 for T in group)
    cutoff = cutoff_build(G, a, b, max_nodes=50_000_000)
    assert cutoff["evidence"]["verdict"] == "proved_zero_by_modular_cutoff"
    modular_check(G, a, b, cutoff["modular"])
    shells, nodes = box_signed_shells(G, a, b, cutoff["modular"]["cutoff"], 50_000_000)
    assert shells == cutoff["evidence"]["enumeration"]["shells"]
    config = kernel_configuration(G, a, b, 5_000_000)
    assert config["level"] == 4
    H1, H2 = [config["coset_lattices"][name]["reduced"]["gram"] for name in ("O_v1", "O_v2")]
    assert not box_isometries(H1, H2, 5_000_000)
    return {"source": str(source), "sha256": hashlib.sha256(_source_bytes(source)).hexdigest(),
            "source_determinant": 4480, "gram_divisor": 2,
            "G": G, "a": a, "b": b, "determinant": 280,
            "primitive": True, "even": False, "full_group_order": len(group),
            "level": 4, "cutoff_document": cutoff, "box_enumeration_nodes": nodes,
            "independent_shells_agree": True, "independent_nonisometry": True,
            "kernel_configuration": config,
            "scope": "rescaling of an already certified example; no global minimality claim"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", choices=("template", "det280"))
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    out = Path(args.output)
    if out.exists():
        raise FileExistsError(out)
    result = template_check() if args.task == "template" else det280_check()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k not in ("coefficient_matrices_for_3G", "cutoff_document", "kernel_configuration")}, indent=2))


if __name__ == "__main__":
    main()
