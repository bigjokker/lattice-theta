"""Adversarial tests for the rank-four parity search and its independent replay."""
from copy import deepcopy
import json
from pathlib import Path
import unittest

from exact_theta import EnumerationLimit, InvalidCertificate, multiply
from rank4_feasibility import candidates as feasibility_candidates
import rank4_parity_search as S
import verify_rank4_parity as V

NODES, ENUM = 5_000_000, 50_000_000


class ParityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.small = S.task_b(29, NODES, 5_000_000)
        cls.certs = json.loads((Path(__file__).resolve().parent.parent / "data/rank4-parity/c1-counterexamples.json").read_text(encoding="utf-8"))
        cls.cert = cls.certs["certificates"][0]

    # --- genuine material ---------------------------------------------------

    def test_generator_matches_frozen_domain_and_census(self):
        self.assertEqual(S.candidate_grams(24, even_only=False, det_filter="all"), feasibility_candidates(24))
        checks = S.generator_checks(max_nodes=NODES)
        self.assertTrue(checks["bijection_with_census"])
        self.assertEqual(checks["determinants"], {5: 1, 9: 1, 13: 1, 17: 1, 21: 2})
        self.assertEqual(V.independent_candidates(24, True, "odd"), S.candidate_grams(24, True, "odd"))

    def test_small_task_b_replays(self):
        result = V.replay_decided(self.small, NODES, ENUM)
        self.assertEqual(result["verdict"], "verified")
        self.assertEqual(result["P_fails"], 0)
        self.assertEqual(result["classes"], self.small["summary"]["classes"])

    def test_counterexample_certificate_replays(self):
        self.assertEqual(len(self.certs["certificates"]), 6)
        result = V.replay_certificate(self.cert, NODES, ENUM)
        self.assertEqual(result["verdict"], "verified_c1_counterexample")
        self.assertEqual((result["determinant"], result["level"], result["witness_orders"]), (121, 1, [4]))

    def test_structure_on_data(self):
        """Constant-sign lemma (a = c), type versus det mod 8, witness orders."""
        for row in self.small["lattices"]:
            self.assertEqual(row["arf"]["type"] == "plus", row["determinant"] % 8 == 1)
            if "parity" not in row:
                continue
            for d in row["parity"]["decisions"]:
                if d["f2"].get("a_equals_some_c"):
                    self.assertEqual(d["verdict"], "nonzero")
                    self.assertTrue(all(abs(s["signed_coefficient"]) == s["vector_count"] for s in d["shells"]))
        for d in self.cert["stabilizer"]:
            if d["epsilon"] == 1:
                self.assertEqual(d["order"], 4)
                self.assertEqual(multiply(d["T"], d["T"]), S.MINUS_I4)
        for z in [d for d in self.cert["f2"].values()]:
            pass
        self.assertTrue(self.cert["f2"]["span_a_c_totally_singular"])
        self.assertEqual(self.cert["arf"]["type"], "plus")

    # --- tampering ----------------------------------------------------------

    def _reject(self, change, exception=InvalidCertificate):
        doc = deepcopy(self.small)
        change(doc)
        with self.assertRaises(exception):
            V.replay_decided(doc, NODES, ENUM)

    def test_class_removed_or_added(self):
        def remove(doc):
            doc["classes_doc"]["representatives"].pop()
            doc["classes_doc"]["classes"] -= 1
            doc["lattices"].pop()
        def add(doc):
            rep = deepcopy(doc["classes_doc"]["representatives"][0])
            doc["classes_doc"]["representatives"].append(rep)
            doc["classes_doc"]["classes"] += 1
            doc["lattices"].append(deepcopy(doc["lattices"][0]))
        def drop_candidate(doc):
            doc["classes_doc"]["candidate_bases"].pop()
            doc["classes_doc"]["candidates"] -= 1
            doc["classes_doc"]["generation"].pop()
        for change in (remove, add, drop_candidate):
            self._reject(change)

    def test_nonzero_with_wrong_bound(self):
        def wrong_bound(doc):
            row = next(r for r in doc["lattices"] if "parity" in r)
            row["parity"]["decisions"][0]["norm4_bound"] += 2
        def wrong_first(doc):
            row = next(r for r in doc["lattices"] if "parity" in r)
            row["parity"]["decisions"][0]["first_nonzero"]["norm4"] += 2
        def wrong_coefficient(doc):
            row = next(r for r in doc["lattices"] if "parity" in r)
            row["parity"]["decisions"][0]["shells"][0]["signed_coefficient"] += 1
        for change in (wrong_bound, wrong_first, wrong_coefficient):
            self._reject(change)

    def test_zero_without_cutoff_or_witness(self):
        doc = deepcopy(self.cert)
        del doc["identity"]["modular_cutoff"]["evidence"]["enumeration"]
        with self.assertRaises(InvalidCertificate):
            V.replay_certificate(doc, NODES, ENUM)
        doc = deepcopy(self.cert)
        doc["identity"]["modular_cutoff"]["run"]["norm4_bound"] -= 1
        with self.assertRaises(InvalidCertificate):
            V.replay_certificate(doc, NODES, ENUM)
        doc = deepcopy(self.cert)
        doc["identity"]["symmetry_witness"]["T"] = S.I4
        with self.assertRaises(InvalidCertificate):
            V.replay_certificate(doc, NODES, ENUM)
        doc = deepcopy(self.cert)
        for s in doc["stabilizer"]:
            s["epsilon"] = 0
        with self.assertRaises(InvalidCertificate):
            V.replay_certificate(doc, NODES, ENUM)

    def test_tampered_automorphism_list(self):
        def drop(doc):
            doc["lattices"][0]["automorphisms"].pop()
            doc["lattices"][0]["aut_order"] -= 1
        self._reject(drop)
        doc = deepcopy(self.cert)
        doc["automorphisms"].pop()
        doc["aut_order"] -= 1
        with self.assertRaises(InvalidCertificate):
            V.replay_certificate(doc, NODES, ENUM)

    def test_decomposable_presented_as_indecomposable(self):
        def flip_decomposable(doc):
            row = next(r for r in doc["lattices"] if r["decomposable"])
            row["decomposable"] = False
            row["eichler_components"] = 1
            row["parity"] = deepcopy(next(r for r in doc["lattices"] if "parity" in r)["parity"])
        def flip_indecomposable(doc):
            row = next(r for r in doc["lattices"] if not r["decomposable"])
            row["decomposable"] = True
            row["eichler_components"] = 2
        for change in (flip_decomposable, flip_indecomposable):
            self._reject(change)
        doc = deepcopy(self.cert)
        doc["indecomposable"] = False
        with self.assertRaises(InvalidCertificate):
            V.replay_certificate(doc, NODES, ENUM)

    def test_wrong_summary_and_status(self):
        def wrong_status(doc):
            row = next(r for r in doc["lattices"] if "parity" in r)
            row["parity"]["P"] = "fails"
        def wrong_summary(doc):
            doc["summary"]["P_holds"] -= 1
        for change in (wrong_status, wrong_summary):
            self._reject(change)

    def test_interrupted_replay_cannot_certify(self):
        with self.assertRaises(EnumerationLimit):
            V.replay_decided(self.small, 10, ENUM)
        with self.assertRaises(EnumerationLimit):
            V.replay_decided(self.small, NODES, 1)
        with self.assertRaises(EnumerationLimit):
            V.replay_certificate(self.cert, 10, ENUM)
        with self.assertRaises(EnumerationLimit):
            V.replay_certificate(self.cert, NODES, 1)


if __name__ == "__main__":
    unittest.main()
