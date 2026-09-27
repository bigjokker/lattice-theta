"""Adversarial census coverage and proof tests using the completed small pilot."""
from copy import deepcopy
import json
from pathlib import Path
import unittest

from exact_theta import EnumerationLimit, InvalidCertificate
from rank4_census import build
from rank4_feasibility import cutoff_build
from verify_rank4_census import replay


class CensusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pilot = json.loads((Path(__file__).resolve().parent.parent / "data/census-rank4/pilot-det4-bound32.json").read_text())
        cls.small = deepcopy(cls.pilot)
        cls.small["domain"]["determinant_bound"] = 1
        cls.small["generation"] = cls.small["generation"][:1]
        cls.small["lattices"] = cls.small["lattices"][:1]
        cls.small["summary"] = {"candidate_bases":1,"lattice_classes":1,"even_pairs":136,
                                "odd_pairs":120,"even_verdicts":{"proved_nonzero":81,"proved_zero":55}}

    def rejects(self, change, pilot=False):
        doc = deepcopy(self.pilot if pilot else self.small)
        change(doc)
        with self.assertRaises(InvalidCertificate):
            replay(doc)

    def test_complete_pilot(self):
        result = replay(self.pilot)
        self.assertEqual(result["summary"]["lattice_classes"],8)
        self.assertIs(result["symmetry_converse_in_domain"],True)

    def test_missing_candidate_and_class(self):
        self.rejects(lambda d: d["generation"].pop())
        self.rejects(lambda d: d["lattices"].pop())

    def test_transport_and_duplicate_class(self):
        self.rejects(lambda d: d["generation"][0]["T"][0].__setitem__(0,2))
        def duplicate(d):
            d["lattices"].append(deepcopy(d["lattices"][0]))
        self.rejects(duplicate)
        def actual_duplicate(d):
            item = d["generation"][-1]
            row = deepcopy(d["lattices"][item["representative"]])
            row["G"] = item["G"]
            item["representative"] = len(d["lattices"])
            item["T"] = [[int(i == j) for j in range(4)] for i in range(4)]
            d["lattices"].append(row)
        self.rejects(actual_duplicate,pilot=True)
        # Redirect the retained determinant-two class to determinant one.
        self.rejects(lambda d: d["generation"][1].__setitem__("representative",0),pilot=True)

    def test_full_group_and_stabilizer(self):
        self.rejects(lambda d: d["lattices"][0]["automorphisms"].pop())
        self.rejects(lambda d: d["lattices"][0]["pairs"][0]["stabilizer"]["indices"].pop())
        self.rejects(lambda d: d["lattices"][0]["pairs"][0]["stabilizer"]["epsilon"].__setitem__(0,1))

    def test_coefficient_and_modular_metadata(self):
        self.rejects(lambda d: d["lattices"][0]["pairs"][0]["enumeration"]["shells"][0].__setitem__("signed_coefficient",2))
        self.rejects(lambda d: d["lattices"][0]["modular_kernels"][0].__setitem__("cutoff",1))

    def test_cancelled_bound_is_unresolved(self):
        doc = build(1,bound=0)
        result = replay(doc)
        self.assertGreater(result["unresolved_pairs"],0)
        self.assertIsNone(result["symmetry_converse_in_domain"])
        target = next(p for p in doc["lattices"][0]["pairs"] if p["verdict"] == "unresolved")
        target["verdict"] = "proved_zero_without_symmetry"
        with self.assertRaises(InvalidCertificate):
            replay(doc)

    def test_modular_fallback_replayed(self):
        doc = deepcopy(self.small)
        row = doc["lattices"][0]
        p = next(p for p in row["pairs"] if p["verdict"] == "proved_nonzero" and any(p["b"]))
        del p["enumeration"]
        del p["evidence"]
        p.update(enumeration_complete=False,reason="synthetic interruption before modular fallback",
                 cutoff_certificate=cutoff_build(row["G"],p["a"],p["b"]),
                 verdict="proved_nonzero_by_modular_replay")
        doc["summary"]["even_verdicts"] = {"proved_nonzero":80,"proved_zero":55,
                                               "proved_nonzero_by_modular_replay":1}
        self.assertEqual(replay(doc)["verdict"],"verified")
        p["cutoff_certificate"]["b"] = [0]*4
        with self.assertRaises(InvalidCertificate):
            replay(doc)

    def test_limits_never_certify_partial_domain(self):
        with self.assertRaises(EnumerationLimit):
            build(4,max_nodes=1)
        limited = replay(self.small,max_nodes=1)
        self.assertEqual(limited["verdict"],"unresolved")
        self.assertIsNone(limited["symmetry_converse_in_domain"])


if __name__ == "__main__":
    unittest.main()
