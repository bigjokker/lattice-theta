"""Rank-three coverage, full-group, independent replay and adversarial checks."""

from copy import deepcopy
import json
from pathlib import Path
import unittest

from exact_theta import EnumerationLimit, InvalidCertificate, identity, multiply, transpose
from lattice_census import classify, stabilizer
from rank3_census import admit, build, candidates, isometries, modular_data
from verify_rank3_census import box_isometries, independent_candidates, replay

ROOT = Path(__file__).resolve().parent.parent


class RankThreeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pack = json.loads((ROOT / "data/census-rank3/det24-bound32.json").read_text(encoding="utf-8"))
        cls.small = build(determinant_bound=4)

    def test_two_generation_algorithms_all_domain_bounds(self):
        for H in range(1, 25):
            with self.subTest(H=H):
                self.assertEqual(candidates(H), independent_candidates(H))
        self.assertEqual(len(candidates()), 303)
        # Boundary shortest-norm-3 lattices are included.
        self.assertTrue(any(G[0][0] == 3 for G in candidates()))
        # Even nonprimitive scaling of Z3 is a separate domain input.
        self.assertIn([[2, 0, 0], [0, 2, 0], [0, 0, 2]], candidates())

    def test_saved_full_census_independent_replay(self):
        box, ldl = replay(self.pack), replay(self.pack, backend="ldl")
        self.assertEqual(box["verdict"], "resolved")
        self.assertTrue(box["certificate_replay_complete"])
        self.assertEqual(box["lattices"], ldl["lattices"])
        self.assertEqual(len(self.pack["lattices"]), 120)
        for row in self.pack["lattices"]:
            self.assertEqual(len(row["pairs"]), 36)
            self.assertEqual(len(row["odd_pairs"]), 28)
        self.assertEqual(self.pack["summary"]["proved_zero"], 723)
        self.assertEqual(self.pack["summary"]["proved_nonzero"], 3597)
        self.assertEqual(self.pack["summary"]["unresolved"], 0)

    def test_noncanonical_bases_and_negative_isometry(self):
        G = [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
        U = [[1, 2, 0], [0, 1, -1], [0, 0, 1]]
        H = multiply(multiply(transpose(U), G), U)
        found = isometries(H, G)
        self.assertIn(U, found)
        self.assertEqual(found, box_isometries(H, G, 1_000_000))
        # Matching determinant is not an isometry certificate.
        A3 = [[2, -1, 0], [-1, 2, -1], [0, -1, 2]]
        diagonal = [[1, 0, 0], [0, 1, 0], [0, 0, 4]]
        self.assertEqual(admit(A3), admit(diagonal))
        self.assertEqual(isometries(A3, diagonal), [])

    def test_full_group_and_decomposable_controls(self):
        for G, order in ((identity(3), 48),
                         ([[2, -1, 0], [-1, 2, -1], [0, -1, 2]], 48),
                         ([[1, 0, 0], [0, 2, 0], [0, 0, 3]], 8)):
            group = isometries(G, G)
            self.assertEqual(len(group), order)
            self.assertEqual(group, box_isometries(G, G, 1_000_000))
        group = isometries(identity(3), identity(3))
        zeros = sum(classify(identity(3), p["a"], p["b"], group, 32, 1_000_000)["verdict"] == "proved_zero"
                    for p in self.small["lattices"][0]["pairs"])
        self.assertEqual(zeros, 9)  # 36 even pairs minus 3^3 nonzero products.

    def test_stabilizing_product_missed_by_ambient_generators(self):
        G = identity(3)
        generators = [[[0, 0, 1], [1, 0, 0], [0, 1, 0]],
                      [[1, 0, 0], [0, 1, 0], [0, 0, -1]]]
        self.assertNotIn(1, stabilizer(G, [1, 1, 0], [1, 1, 0], generators)["epsilon"])
        self.assertIn(1, stabilizer(G, [1, 1, 0], [1, 1, 0], isometries(G, G))["epsilon"])

    def test_resource_limits_and_bounded_cancellation(self):
        with self.assertRaises(EnumerationLimit):
            isometries(identity(3), identity(3), max_nodes=1)
        limited = build(determinant_bound=4, max_nodes=1)
        self.assertFalse(limited["summary"]["groups_complete"])
        self.assertFalse(limited["summary"]["uniqueness_complete"])
        self.assertFalse(limited["summary"]["symmetry_converse_proved_in_domain"])
        self.assertTrue(all(p["verdict"] == "unresolved" and "enumeration" not in p
                            for r in limited["lattices"] for p in r["pairs"]))
        self.assertEqual(replay(limited)["verdict"], "unresolved")
        self.assertFalse(replay(self.small, max_nodes=1)["certificate_replay_complete"])
        shallow = build(determinant_bound=1, norm4_bound=0)
        self.assertGreater(shallow["summary"]["unresolved"], 0)
        self.assertEqual(replay(shallow)["verdict"], "unresolved")
        trap = json.loads((ROOT / "data/ternary.json").read_text(encoding="utf-8"))
        self.assertEqual(classify(trap["G"], trap["a"], trap["b"], None, 7, 1_000_000)["verdict"], "unresolved")
        self.assertEqual(classify(trap["G"], trap["a"], trap["b"], None, 11, 1_000_000)["verdict"], "proved_nonzero")
        # A first-shell trap inside the new domain (determinant 20).
        G = [[2, -1, -1], [-1, 4, -1], [-1, -1, 4]]
        group = isometries(G, G)
        self.assertEqual(classify(G, [0, 1, 1], [1, 0, 0], group, 6, 1_000_000)["verdict"], "unresolved")
        later = classify(G, [0, 1, 1], [1, 0, 0], group, 10, 1_000_000)
        self.assertEqual(later["evidence"], {"kind": "shell", "norm4": 10, "coefficient": 2})

    def test_modular_levels_not_determinants(self):
        item = modular_data(identity(3), [1, 0, 0])
        self.assertEqual(item["G0"], [[4, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.assertEqual((item["d0"], item["N0"], item["unsigned_level"]), (4, 4, 64))
        self.assertIsNone(item["cutoff"])
        self.assertEqual(modular_data([[2, 0, 0], [0, 2, 0], [0, 0, 2]], [1, 0, 0])["N0"], 8)
        self.assertIsNone(modular_data(identity(3), [0, 0, 0]))

    def test_reject_coverage_group_sign_shell_and_metadata_tampering(self):
        mutations = [lambda d: d["generation"].pop(),
                     lambda d: d["lattices"][0]["pairs"].pop(),
                     lambda d: d["lattices"][0]["odd_pairs"].pop(),
                     lambda d: d["lattices"][0]["automorphisms"].pop(),
                     lambda d: d["lattices"][0]["pairs"][0]["stabilizer"]["epsilon"].__setitem__(0, 1),
                     lambda d: d["lattices"][0]["pairs"][0]["evidence"].__setitem__("coefficient", 9),
                     lambda d: d["lattices"][0]["pairs"][0]["a"].__setitem__(0, False),
                     lambda d: d["lattices"][0]["modular_inputs"][0].__setitem__("N0", 1),
                     lambda d: d["lattices"][0]["modular_inputs"][0].__setitem__("cutoff", 1),
                     lambda d: d["generation"][0]["U"][0].__setitem__(0, 2),
                     lambda d: d["summary"].__setitem__("proved_zero", 0)]
        for mutate in mutations:
            doc = deepcopy(self.small)
            mutate(doc)
            with self.assertRaises(InvalidCertificate):
                replay(doc)
        for value in (0, 25, True, 1.0):
            with self.assertRaises(InvalidCertificate):
                candidates(value)
        for G in ([[1, 0], [0, 1]], [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
                  [[True, 0, 0], [0, 1, 0], [0, 0, 1]]):
            with self.assertRaises(InvalidCertificate):
                admit(G)


if __name__ == "__main__":
    unittest.main()
