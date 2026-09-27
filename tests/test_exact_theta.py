"""Correctness checks for mathematical certificates, including false positives."""

from copy import deepcopy
from itertools import product
import json
from pathlib import Path
import unittest

from exact_theta import (InvalidCertificate, complete_shells, dot, identity, inverse,
                         matvec, matrix_power, multiply, symmetry_data, transpose, verify)

ROOT = Path(__file__).resolve().parent.parent


def load(name):
    return json.loads((ROOT / "data" / f"{name}.json").read_text(encoding="utf-8"))


class CertificateTests(unittest.TestCase):
    def test_example_pack(self):
        for name in ("d4", "a3", "ternary", "order4", "order8", "order16"):
            with self.subTest(name=name):
                result = verify(load(name))
                self.assertIn(result["verdict"], ("vanishes_by_symmetry", "does_not_vanish"))
                self.assertEqual(len(result["gram_sha256"]), 64)

    def test_minimal_shell_false_positive(self):
        doc = load("ternary")
        result = verify(doc)
        self.assertEqual(result["enumeration"]["shells"], [
            {"norm4": 7, "norm": "7/4", "vector_count": 4, "signed_coefficient": 0},
            {"norm4": 11, "norm": "11/4", "vector_count": 2, "signed_coefficient": -2}])
        doc["evidence"] = {"kind": "search", "norm4_bound": 7}
        self.assertEqual(verify(doc)["verdict"], "unresolved")
        doc["evidence"]["norm4_bound"] = 11
        self.assertEqual(verify(doc)["verdict"], "does_not_vanish")

    def test_reject_wrong_symmetry_or_shell(self):
        doc = load("d4")
        doc["evidence"]["T"][0][0] += 1
        with self.assertRaises(InvalidCertificate):
            verify(doc)
        doc = load("d4")
        doc["evidence"]["T"] = identity(4)
        with self.assertRaisesRegex(InvalidCertificate, "epsilon is zero"):
            verify(doc)
        doc = load("a3")
        doc["evidence"]["coefficient"] = 1
        with self.assertRaisesRegex(InvalidCertificate, "coefficient mismatch"):
            verify(doc)
        doc["evidence"]["coefficient"] = 0
        with self.assertRaises(InvalidCertificate):
            verify(doc)

    def test_require_both_stabilizers(self):
        doc = load("order4")
        doc["a"] = [1, 0, 0, 0]
        with self.assertRaisesRegex(InvalidCertificate, "stabilize a"):
            verify(doc)
        doc = load("order4")
        doc["b"] = [1, 0, 0, 0]
        with self.assertRaisesRegex(InvalidCertificate, "stabilize b"):
            verify(doc)

    def test_positive_definite_and_integral_inputs(self):
        for G in ([[1, 2], [2, 1]], [[1, 1], [1, 1]], [[1, 1], [0, 1]], [[1.0]]):
            doc = {"schema": "work7-theta-v1", "G": G, "a": [0] * len(G),
                   "b": [0] * len(G), "evidence": {"kind": "search", "norm4_bound": 1}}
            with self.subTest(G=G), self.assertRaises(InvalidCertificate):
                verify(doc)
        doc = load("a3")
        doc["a"][0] = True
        with self.assertRaises(InvalidCertificate):
            verify(doc)

    def test_character_representatives_inverse_and_composition(self):
        for name in ("d4", "order4", "order8", "order16"):
            doc = load(name)
            G, a, b, T = (doc[k] for k in ("G", "a", "b", "evidence"))
            T = T["T"]
            n = len(G)
            invT = [[int(x) for x in row] for row in inverse(T)]
            self.assertEqual(symmetry_data(G, a, b, invT)["epsilon"], 1)
            for k in (1, 2, 3):
                self.assertEqual(symmetry_data(G, a, b, matrix_power(T, k))["epsilon"], k % 2)
            # Independently check the legacy T^t formula on a stabilizing T.
            delta = [x - y for x, y in zip(matvec(transpose(T), b), b)]
            self.assertTrue(all(x % 2 == 0 for x in delta))
            self.assertEqual(dot(a, [x // 2 for x in delta]) % 2, 1)
            shifted = deepcopy(doc)
            shifted["a"] = [x + 2 * (i - 2) for i, x in enumerate(a)]
            shifted["b"] = [x + 2 * (3 - i) for i, x in enumerate(b)]
            self.assertEqual(verify(shifted)["symmetry"]["epsilon"], 1)

    def test_order_and_unbounded_integer_size(self):
        doc = load("order8")
        doc["evidence"]["order"] = 16
        with self.assertRaisesRegex(InvalidCertificate, "not minimal"):
            verify(doc)
        doc["evidence"]["order"] = 4
        with self.assertRaises(InvalidCertificate):
            verify(doc)
        doc = load("order4")
        doc["G"] = [[10 ** 50 * x for x in row] for row in doc["G"]]
        self.assertEqual(verify(doc)["verdict"], "vanishes_by_symmetry")

    def test_resource_limit_returns_no_partial_coefficient(self):
        result = verify(load("ternary"), max_nodes=1)
        self.assertEqual(result["verdict"], "unresolved")
        self.assertFalse(result["enumeration_complete"])
        self.assertNotIn("enumeration", result)
        self.assertNotIn("signed_coefficient", result)

    def test_enumeration_against_independent_complete_boxes(self):
        # The inverse-Gram Cauchy bound proves completeness of the reference box:
        # |y_i|^2 <= (y^t G y)(G^-1)_ii. All chosen cases have |y_i| <= 8.
        for G in ([[2]], [[3, 1], [1, 2]], [[2, -1, 0], [-1, 2, -1], [0, -1, 2]]):
            n = len(G)
            bound = 12
            invG = inverse(G)
            self.assertTrue(all(bound * invG[i][i] < 81 for i in range(n)))
            for a in product((0, 1), repeat=n):
                b = [i % 2 for i in range(n)]
                expected = {}
                for y in product(range(-8, 9), repeat=n):
                    if any((yi - ai) % 2 for yi, ai in zip(y, a)):
                        continue
                    q = dot(y, matvec(G, y))
                    if q > bound:
                        continue
                    x = [(yi - ai) // 2 for yi, ai in zip(y, a)]
                    count, coefficient = expected.get(q, (0, 0))
                    expected[q] = (count + 1, coefficient + (-1 if dot(x, b) % 2 else 1))
                actual = complete_shells(G, list(a), b, bound)
                self.assertEqual({s["norm4"]: (s["vector_count"], s["signed_coefficient"])
                                  for s in actual["shells"]}, expected)

    def test_change_of_basis(self):
        doc = load("d4")
        U = [[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1], [0, 0, 0, 1]]
        Ui = [[int(x) for x in row] for row in inverse(U)]
        doc["G"] = multiply(multiply(transpose(U), doc["G"]), U)
        doc["a"] = matvec(Ui, doc["a"])
        doc["b"] = matvec(transpose(U), doc["b"])
        doc["evidence"]["T"] = multiply(multiply(Ui, doc["evidence"]["T"]), U)
        self.assertEqual(verify(doc)["verdict"], "vanishes_by_symmetry")
        doc = load("ternary")
        U = [[1, 1, 0], [0, 1, 1], [0, 0, 1]]
        Ui = [[int(x) for x in row] for row in inverse(U)]
        doc["G"] = multiply(multiply(transpose(U), doc["G"]), U)
        doc["a"] = matvec(Ui, doc["a"])
        doc["b"] = matvec(transpose(U), doc["b"])
        self.assertEqual(verify(doc)["signed_coefficient"], -2)


if __name__ == "__main__":
    unittest.main()
