"""Exact root-lattice identifications, coverage, and E8 scope regression."""

from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import unittest

from dn_exact import basis, gram as dn_gram
from exact_theta import InvalidCertificate, identity, inverse, matvec, multiply, transpose, verify
from exceptional_ambient import doubled_roots
from nonsimply_laced import (A2_IN_G2, F4_G, F4_TO_D4, G2_G,
                            bn_audit, verify_root_pack)

ROOT = Path(__file__).resolve().parent.parent


class RootLatticeTests(unittest.TestCase):
    def test_bn_complete_minimal_shell_audits(self):
        saved = json.loads((ROOT / "reports" / "root-lattices-audit.json").read_text(encoding="utf-8"))
        for rank in range(1, 9):
            self.assertEqual(bn_audit(rank), saved[rank - 1])

    def test_cn_simple_basis_generates_dn(self):
        for rank in range(2, 11):
            U = identity(rank)
            U[-2][-1] = -1
            # The last Cn simple root is 2e_n; the other columns are differences.
            C = multiply(basis(rank), U)
            self.assertEqual([row[-1] for row in C], [0] * (rank - 1) + [2])
            self.assertTrue(all(v.denominator == 1 for row in inverse(U) for v in row))
            G = multiply(multiply(transpose(U), dn_gram(rank)), U)
            self.assertEqual([G[i][i] for i in range(rank)], [2] * (rank - 1) + [4])
            self.assertEqual(G[-2][-1], -2)

    def test_f4_and_g2_unimodular_isometries(self):
        for U in (F4_TO_D4, A2_IN_G2):
            self.assertTrue(all(v.denominator == 1 for row in inverse(U) for v in row))
        self.assertEqual(multiply(multiply(transpose(F4_TO_D4), dn_gram(4)), F4_TO_D4), F4_G)
        self.assertEqual(multiply(multiply(transpose(A2_IN_G2), G2_G), A2_IN_G2), [[2, -1], [-1, 2]])
        # Bourbaki F4 simple roots in the basis (s,e2,e3,e4).
        K = [[0, 0, 0, 1], [1, 0, 0, -1], [-1, 1, 0, -1], [0, -1, 1, -1]]
        self.assertTrue(all(v.denominator == 1 for row in inverse(K) for v in row))
        self.assertEqual(multiply(multiply(transpose(K), F4_G), K),
                         [[4, -2, 0, 0], [-2, 4, -2, 0], [0, -2, 2, -1], [0, 0, -1, 2]])

    def test_saved_packs_have_complete_exact_coverage(self):
        for family, zeros, nonzeros in (("f4", 9, 127), ("g2", 0, 10)):
            pack = json.loads((ROOT / "data" / "root-lattices" / f"{family}-classes.json").read_text(encoding="utf-8"))
            result = verify_root_pack(pack)
            self.assertEqual((result["even_zero"], result["even_nonzero"]), (zeros, nonzeros))
            for change in (lambda p: p["certificates"].pop(),
                           lambda p: p["certificates"].append(p["certificates"][0]),
                           lambda p: p["certificates"].__setitem__(1, p["certificates"][0])):
                corrupt = deepcopy(pack)
                change(corrupt)
                with self.assertRaises(InvalidCertificate):
                    verify_root_pack(corrupt)

    def test_e8_quarter_coordinate_characteristic_is_a_certified_zero(self):
        doc = json.loads((ROOT / "data" / "root-lattices" / "e8-quarter-xi-zero.json").read_text(encoding="utf-8"))
        C = doubled_roots(8)
        self.assertEqual(matvec(C, doc["a"]), [5, 1, 1, 1, 1, 1, 1, 1])
        self.assertEqual(matvec(inverse(transpose(C)), [4 * v for v in doc["b"]]),
                         [2, 2, 2, 2, 0, 0, 0, 0])
        result = verify(doc)
        self.assertEqual(result["pairing"], 4)
        self.assertEqual(result["verdict"], "vanishes_by_symmetry")

    def test_bn_standalone_symmetry_certificates(self):
        for rank in (2, 8):
            doc = json.loads((ROOT / "data" / "root-lattices" / f"b{rank}-zero.json").read_text(encoding="utf-8"))
            result = verify(doc)
            self.assertEqual(result["symmetry"]["witness_order"], 2)
            self.assertEqual(result["verdict"], "vanishes_by_symmetry")


if __name__ == "__main__":
    unittest.main()
