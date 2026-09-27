"""Regressions for the uniform binary theorem; finite tests are not its proof."""

from fractions import Fraction
import json
from pathlib import Path
import unittest

from exact_theta import (complete_shells, inverse, matvec, multiply,
                         symmetry_data, transpose)
from lattice_census import even_pairs
from verify_lattice_census import box_shells

ROOT = Path(__file__).resolve().parent.parent


class RankTwoTests(unittest.TestCase):
    def check_critical_shells(self, A, B, C):
        G = [[A, B], [B, C]]
        cases = [([1, 0], [0, 1], A, 2, 2),
                 ([0, 1], [1, 0], C, 2, 2),
                 ([1, 1], [1, 1], A + C - 2 * B, 2 if B else 4, -2 if B else 0)]
        for a, b, bound, count, coeff in cases:
            with self.subTest(G=G, a=a, b=b):
                expected = [{"norm4": bound, "norm": str(Fraction(bound, 4)),
                             "vector_count": count, "signed_coefficient": coeff}]
                self.assertEqual(complete_shells(G, a, b, bound)["shells"], expected)
                self.assertEqual(box_shells(G, a, b, bound), expected)

    def test_reduced_forms_and_reduction_boundaries(self):
        for A in range(1, 9):
            for B in range(A // 2 + 1):
                for C in (A, A + 1, 2 * A + 3):
                    self.check_critical_shells(A, B, C)

    def test_large_determinants_and_unequal_axes(self):
        # The second shell needs the asymmetric C bound; A is not a proxy for C.
        for A, B, C in ((100, 50, 100), (2, 1, 10000), (99, 49, 10001),
                         (101, 1, 10000), (1, 0, 10000), (4, 2, 4)):
            self.check_critical_shells(A, B, C)

    def test_all_even_pairs_and_the_unique_diagonal_zero(self):
        for A, B, C in ((1, 0, 24), (7, 0, 7), (7, 0, 100),
                         (2, 1, 2), (3, 1, 100), (4, 2, 4), (5, 2, 6)):
            G = [[A, B], [B, C]]
            zero = []
            for a, b in even_pairs(2):
                shells = complete_shells(G, a, b, A + C + 2 * B)["shells"]
                if B == 0 and a == b == [1, 1]:
                    zero.append((a, b))
                    self.assertTrue(all(s["signed_coefficient"] == 0 for s in shells))
                    self.assertEqual(symmetry_data(G, a, b, [[-1, 0], [0, 1]])["epsilon"], 1)
                else:
                    self.assertTrue(any(s["signed_coefficient"] for s in shells), (G, a, b))
            self.assertEqual(len(zero), int(B == 0))

    def test_factor_swap_is_not_the_even_zero_witness(self):
        G = [[7, 0], [0, 7]]
        a = b = [1, 1]
        self.assertEqual(symmetry_data(G, a, b, [[0, 1], [1, 0]])["epsilon"], 0)
        self.assertEqual(symmetry_data(G, a, b, [[-1, 0], [0, 1]])["epsilon"], 1)

    def test_basis_changes_preserve_the_prediction(self):
        for G in ([[2, 1], [1, 2]], [[7, 0], [0, 7]]):
            U = [[1, 1], [0, 1]]
            Ui = [[int(x) for x in row] for row in inverse(U)]
            H = multiply(multiply(transpose(U), G), U)
            a = matvec(Ui, [1, 1])
            b = matvec(transpose(U), [1, 1])
            expected = complete_shells(G, [1, 1], [1, 1], 2 * G[0][0] - 2 * G[0][1])["shells"]
            actual = complete_shells(H, a, b, 2 * G[0][0] - 2 * G[0][1])["shells"]
            self.assertEqual(actual, expected)
            if G[0][1] == 0:
                T = multiply(multiply(Ui, [[-1, 0], [0, 1]]), U)
                self.assertEqual(symmetry_data(H, a, b, T)["epsilon"], 1)

    def test_uniform_theorem_matches_the_preserved_census(self):
        pack = json.loads((ROOT / "data/census/rank12-det24-closed.json").read_text())
        for row in pack["lattices"]:
            G = row["G"]
            for p in row["pairs"]:
                expected_zero = len(G) == 2 and G[0][1] == 0 and p["a"] == p["b"] == [1, 1]
                self.assertEqual(p["verdict"], "proved_zero" if expected_zero else "proved_nonzero")


if __name__ == "__main__":
    unittest.main()
