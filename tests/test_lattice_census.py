"""Coverage, independent replay, adversarial data and resource regressions."""

from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
import unittest

from exact_theta import (EnumerationLimit, InvalidCertificate, identity,
                         multiply, transpose, verify)
from lattice_census import (admit, build, candidates, isometries,
                            stabilizer)
from verify_lattice_census import box_isometries, replay

ROOT = Path(__file__).resolve().parent.parent


class CensusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pack = build()

    def test_finite_domain_and_saved_reproducibility(self):
        self.assertEqual(len(candidates()), 94)
        self.assertEqual(self.pack, json.loads((ROOT / "data/census/rank12-det24.json").read_text()))
        self.assertEqual(self.pack["summary"]["rank_counts"], {"1": 24, "2": 70})
        self.assertEqual(self.pack["summary"]["even_pairs"], 772)
        self.assertEqual(self.pack["summary"]["proved_zero"], 44)
        self.assertEqual(self.pack["summary"]["proved_nonzero"], 695)
        self.assertEqual(self.pack["summary"]["unresolved"], 33)
        for row in self.pack["lattices"]:
            self.assertEqual(len(row["pairs"]), 2 ** (len(row["G"]) - 1) * (2 ** len(row["G"]) + 1))

    def test_independent_complete_replay(self):
        box = replay(self.pack)
        ldl = replay(self.pack, backend="ldl")
        self.assertTrue(box["certificate_replay_complete"])
        self.assertEqual(box["lattices"], ldl["lattices"])
        self.assertEqual(box["verdict"], "unresolved")
        for row in box["lattices"]:
            self.assertTrue(row["group_replayed"])
            for pair in row["pairs"]:
                if pair["verdict"] == "unresolved":
                    self.assertEqual(pair["shells"], [])

    def test_separate_closure_run_resolves_every_pair(self):
        closed = build(norm4_bound=25)
        self.assertEqual(closed, json.loads((ROOT / "data/census/rank12-det24-closed.json").read_text()))
        self.assertEqual(closed["summary"]["proved_nonzero"], 728)
        self.assertEqual(closed["summary"]["proved_zero"], 44)
        self.assertEqual(closed["summary"]["unresolved"], 0)
        self.assertEqual(closed["summary"]["parity_converse_proved"], 50)
        box, ldl = replay(closed), replay(closed, backend="ldl")
        self.assertEqual(box["verdict"], "resolved")
        self.assertTrue(box["certificate_replay_complete"])
        self.assertEqual(box["lattices"], ldl["lattices"])
        newly_resolved = 0
        for before, after in zip(self.pack["lattices"], closed["lattices"]):
            for key in ("G", "automorphisms", "minimum_squared_norm"):
                self.assertEqual(before[key], after[key])
            for old, new in zip(before["pairs"], after["pairs"]):
                self.assertEqual(old["stabilizer"], new["stabilizer"])
                if old["verdict"] != "unresolved":
                    self.assertEqual(old["verdict"], new["verdict"])
                    self.assertEqual(old["evidence"], new["evidence"])
                    continue
                newly_resolved += 1
                self.assertEqual(new["verdict"], "proved_nonzero")
                # Also replay each new claim through the original standalone verifier.
                self.assertEqual(verify({"schema": "work7-theta-v1", "G": after["G"],
                                         "a": new["a"], "b": new["b"],
                                         "evidence": new["evidence"]})["verdict"], "does_not_vanish")
        self.assertEqual(newly_resolved, 33)

    def test_noncanonical_bases_and_negative_isometry(self):
        G = [[2, 1], [1, 3]]
        U = [[1, 2], [0, 1]]
        H = multiply(multiply(transpose(U), G), U)
        images = isometries(H, G)
        self.assertIn(U, images)
        self.assertEqual(images, box_isometries(H, G, 1_000_000))
        # Same determinant is insufficient to deduplicate these lattices.
        self.assertEqual(admit(G), admit([[1, 0], [0, 5]]))
        self.assertEqual(isometries(G, [[1, 0], [0, 5]]), [])
        self.assertEqual(isometries([[2, 1], [1, 2]], [[2, -1], [-1, 2]]),
                         box_isometries([[2, 1], [1, 2]], [[2, -1], [-1, 2]], 1_000_000))

    def test_full_group_controls(self):
        for G, order in (([[7]], 2), ([[1, 0], [0, 2]], 4),
                         ([[2, 0], [0, 2]], 8), ([[2, 1], [1, 2]], 12)):
            self.assertEqual(len(isometries(G, G)), order)
        square = next(r for r in self.pack["lattices"] if r["G"] == [[2, 0], [0, 2]])
        self.assertEqual(sum(p["verdict"] == "proved_zero" for p in square["pairs"]), 1)
        a2 = next(r for r in self.pack["lattices"] if r["G"] == [[2, 1], [1, 2]])
        self.assertTrue(all(p["verdict"] == "proved_nonzero" for p in a2["pairs"]))

    def test_ambient_generators_can_miss_stabilizing_product(self):
        # Cycle and flip-third generate all signed cyclic permutations in rank 3.
        # Neither supplied stabilizing generator witnesses this even pair, but
        # a conjugate flip does. This guards the full-element stabilizer scan.
        G = identity(3)
        cycle = [[0, 0, 1], [1, 0, 0], [0, 1, 0]]
        flip = [[1, 0, 0], [0, 1, 0], [0, 0, -1]]
        generators = [cycle, flip]
        group, pending = [identity(3)], [identity(3)]
        while pending:
            T = pending.pop()
            for S in generators:
                P = multiply(T, S)
                if P not in group:
                    group.append(P)
                    pending.append(P)
        a = b = [1, 1, 0]
        self.assertNotIn(1, stabilizer(G, a, b, generators)["epsilon"])
        self.assertIn(1, stabilizer(G, a, b, group)["epsilon"])

    def test_resource_limits_retain_uncertainty(self):
        with self.assertRaises(EnumerationLimit):
            isometries([[1, 0], [0, 1]], [[1, 0], [0, 1]], max_nodes=1)
        limited = build(determinant_bound=5, max_nodes=1)
        self.assertFalse(limited["summary"]["groups_complete"])
        self.assertFalse(limited["summary"]["uniqueness_complete"])
        self.assertTrue(limited["summary"]["generation_complete"])
        self.assertTrue(all(p["verdict"] == "unresolved" and "enumeration" not in p
                            for r in limited["lattices"] for p in r["pairs"]))
        result = replay(limited)
        self.assertEqual(result["verdict"], "unresolved")
        result = replay(self.pack, max_nodes=1)
        self.assertFalse(result["certificate_replay_complete"])
        self.assertEqual(result["verdict"], "unresolved")

    def test_reject_malformed_coverage_groups_and_evidence(self):
        def check(mutate):
            doc = deepcopy(self.pack)
            mutate(doc)
            with self.assertRaises(InvalidCertificate):
                replay(doc)
        check(lambda d: d["generation"].pop())
        check(lambda d: d["lattices"][0]["pairs"].pop())
        check(lambda d: d["lattices"][0]["pairs"].append(deepcopy(d["lattices"][0]["pairs"][0])))
        check(lambda d: d["lattices"][0]["automorphisms"].pop())
        check(lambda d: d["lattices"][0]["pairs"][0]["stabilizer"]["epsilon"].__setitem__(0, 1))
        check(lambda d: d["lattices"][0]["pairs"][0]["evidence"].__setitem__("coefficient", 7))
        check(lambda d: d["lattices"][0].__setitem__("G", [[True]]))
        check(lambda d: d["summary"].__setitem__("unresolved", 0))
        check(lambda d: d.__setitem__("lattices", None))
        for G in ([[1.0]], [[True]], [[-1, 0], [0, -1]], [[1, 1], [0, 1]], [[25]]):
            with self.assertRaises(InvalidCertificate):
                admit(G)
        for value in (0, 25, True, 1.0):
            with self.assertRaises(InvalidCertificate):
                candidates(value)

    def test_ternary_first_shell_trap_remains_unresolved(self):
        doc = json.loads((ROOT / "data/ternary.json").read_text())
        doc["evidence"] = {"kind": "search", "norm4_bound": 7}
        self.assertEqual(verify(doc)["verdict"], "unresolved")
        doc["evidence"]["norm4_bound"] = 11
        self.assertEqual(verify(doc)["verdict"], "does_not_vanish")

    def test_cli_unresolved_exit_code(self):
        result = subprocess.run([sys.executable, str(ROOT / "code" / "verify_lattice_census.py"),
                                 str(ROOT / "data/census/rank12-det24.json")],
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertTrue(json.loads(result.stdout)["certificate_replay_complete"])


if __name__ == "__main__":
    unittest.main()
