"""Exact checks supporting reviews/project-assessment-2026-10-21.md.

Run from the repository root:
    python code/check_small_pairs.py

The all-level conclusion uses the classical even-lattice theta theorem and
the weight-two Sturm bound, explained in the write-up. This script checks
their numerical hypotheses and coefficients with two complete enumerators.
"""

import json
from math import isqrt, lcm
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from exact_theta import complete_shells, inverse, ldl
from rank3_cutoff import index_gamma0
from rank4_feasibility import determinant, isometries
from verify_lattice_census import box_shells
from verify_rank4_feasibility import box_isometries


PAIRS = {
    "conway_sloane_remark_x_det961": [
        [[17, -5, 5, 2], [-5, 13, 6, 5], [5, 6, 9, 5], [2, 5, 5, 5]],
        [[17, -2, 1, 7], [-2, 13, 9, 1], [1, 9, 9, 2], [7, 1, 2, 5]],
    ],
    "index_two_kernels_det3844": [
        [[4, 1, 0, 0], [1, 8, 0, 0], [0, 0, 8, 2], [0, 0, 2, 16]],
        [[4, -2, 1, 1], [-2, 8, -3, -1], [1, -3, 10, -4], [1, -1, -4, 18]],
    ],
}


def check_pair(grams):
    # Scaling both Grams by two preserves equality of their theta series.
    scale = 1 if all(G[i][i] % 2 == 0 for G in grams for i in range(4)) else 2
    even_grams = [[[scale * x for x in row] for row in G] for G in grams]
    levels = []
    for A in even_grams:
        assert all(d > 0 for d in ldl(A)[1])
        d = determinant(A)
        assert d == isqrt(d) ** 2  # Weight 2: character (det(A)/.) is trivial.
        inv = inverse(A)
        N = lcm(*(v.denominator for row in inv for v in row),
                *((inv[i][i] / 2).denominator for i in range(4)))
        assert all((N * v).denominator == 1 for row in inv for v in row)
        assert all(int(N * inv[i][i]) % 2 == 0 for i in range(4))
        levels.append(N)
    level = lcm(*levels)
    index = index_gamma0(level)
    bound = index // 6
    coefficients = []
    for A in even_grams:
        # q = exp(2*pi*i*tau); exponent is x^t A x / 2.
        # complete_shells enumerates y=2x, so norm4 = 8*exponent.
        shells = complete_shells(A, [0]*4, [0]*4, 8*bound, 5_000_000)["shells"]
        assert shells == box_shells(A, [0]*4, [0]*4, 8*bound, 5_000_000)
        assert all(s["norm4"] % 8 == 0 for s in shells)
        counts = {s["norm4"] // 8: s["vector_count"] for s in shells}
        coefficients.append([counts.get(n, 0) for n in range(bound + 1)])
    assert coefficients[0] == coefficients[1]
    assert not isometries(*grams, max_nodes=5_000_000)
    assert not box_isometries(*grams, max_nodes=5_000_000)
    return {
        "G": grams, "determinants": [determinant(G) for G in grams],
        "gram_scale_for_evenness": scale, "even_lattice_levels": levels,
        "common_level": level, "weight": 2, "character": "trivial",
        "gamma0_index": index, "inclusive_sturm_bound": bound,
        "coefficients_0_through_bound": coefficients[0],
        "independent_enumerators_agree": True,
        "nonisometry_confirmed_by_ldl_and_box": True,
        "conclusion": "equal theta series by the stated modular-form theorem",
    }


if __name__ == "__main__":
    print(json.dumps({name: check_pair(grams) for name, grams in PAIRS.items()}, indent=2))
