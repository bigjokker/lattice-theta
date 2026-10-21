"""Adversarial tests for the Eisenstein level-2 search and its independent replay."""
from copy import deepcopy
import unittest

from exact_theta import EnumerationLimit, InvalidCertificate
import rank4_eisenstein_search as S
import verify_rank4_eisenstein as V


class EisensteinTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.task_a = S.task_a(trials=10)
        cls.task_b = S.task_b(34, 34)
        rep = next(r for r in cls.task_b["representatives"] if r["level"] == 2)
        cls.member = cls.task_b["members"][rep["member_index"]]
        cls.certificate = S.certify_member(cls.member)
        cls.task_e = S.task_e()

    # --- genuine material ---------------------------------------------------

    def test_formulas_and_control(self):
        self.assertTrue(self.task_a["all_passed"])
        self.assertEqual(self.task_a["control"]["triple_counts"], {"5..41": 498, "5..154": 5539})
        self.assertEqual(self.task_a["control"]["tally"], {"level_2": 15, "nonzero": 90})

    def test_task_b_at_34_finds_the_level_two_example_and_replays(self):
        s = self.task_b["summary"]
        self.assertEqual(s["reduced_triples"], 42)
        self.assertEqual(s["level_2_examples"], 1)
        self.assertEqual(s["smallest_level_2_primitive_determinant"], 1156)
        result = V.replay_enumeration(self.task_b, 2_000_000, 200_000_000)
        self.assertEqual(result["verdict"], "verified")
        self.assertEqual(result["level_2_examples"], 1)

    def test_certificate_is_level_two_with_three_identity_routes(self):
        c = self.certificate
        self.assertEqual(c["determinant"], 1156)
        self.assertEqual(c["G"], [[4, -2, 1, -2], [-2, 8, 2, -2], [1, 2, 6, -4], [-2, -2, -4, 12]])
        self.assertEqual(len(c["automorphisms"]), 2)
        self.assertEqual(c["kernel_automorphisms"]["elements_swapping_v1_v2"], 0)
        self.assertTrue(c["identity"]["modular_cutoff"]["attempted"])
        self.assertTrue(c["identity"]["manuscript_schema"]["available"])
        result = V.replay_certificate(c, 2_000_000, 200_000_000)
        self.assertEqual(result["verdict"], "verified_level_2")
        self.assertEqual(len(result["identity_routes"]), 3)

    def test_task_e_lemma_replays(self):
        self.assertEqual(V.replay_task_e(self.task_e, 60, 2_000_000, 200_000_000)["verdict"], "verified")

    # --- tampering ----------------------------------------------------------

    def _reject(self, change, exception=InvalidCertificate):
        doc = deepcopy(self.certificate)
        change(doc)
        with self.assertRaises(exception):
            V.replay_certificate(doc, 2_000_000, 200_000_000)

    def test_altered_gram_matrices(self):
        self._reject(lambda d: d["G"][0].__setitem__(0, 6))
        self._reject(lambda d: d["G0"][1].__setitem__(1, d["G0"][1][1] + 2))
        self._reject(lambda d: d["G0_standard"][0].__setitem__(0, d["G0_standard"][0][0] + 2))

    def test_rotation_tampering(self):
        self._reject(lambda d: d["identity"]["kernel_rotation"]["R0"][0].__setitem__(0, 1))
        self._reject(lambda d: d["identity"]["kernel_rotation"].__setitem__("maps_v1_to_v2", not d["identity"]["kernel_rotation"]["maps_v1_to_v2"]))
        self._reject(lambda d: d["classes"]["v2"].__setitem__(0, 1 - d["classes"]["v2"][0]))

    def test_forged_swap_and_group_tampering(self):
        self._reject(lambda d: d["kernel_automorphisms"].__setitem__("elements_swapping_v1_v2", 1))
        self._reject(lambda d: d["kernel_automorphisms"].__setitem__("order", d["kernel_automorphisms"]["order"] + 1))
        self._reject(lambda d: d["automorphisms"].pop())
        self._reject(lambda d: d["stabilizer"][0].__setitem__("epsilon", 1))

    def test_truncated_or_altered_enumeration(self):
        def truncate(d):
            d["identity"]["modular_cutoff"]["ldl_backend"]["evidence"]["enumeration"]["shells"].pop()
        self._reject(truncate)

        def nonzero(d):
            d["identity"]["modular_cutoff"]["ldl_backend"]["evidence"]["enumeration"]["shells"][0]["signed_coefficient"] = 2
        self._reject(nonzero)

        def shorter(d):
            d["identity"]["modular_cutoff"]["ldl_backend"]["run"]["norm4_bound"] -= 1
        self._reject(shorter)

        def incomplete(d):
            d["identity"]["modular_cutoff"]["ldl_backend"]["evidence"]["enumeration_complete"] = False
        self._reject(incomplete)

    def test_interrupted_replay_cannot_certify(self):
        with self.assertRaises(EnumerationLimit):
            V.replay_certificate(self.certificate, 2_000_000, 10)
        with self.assertRaises(EnumerationLimit):
            V.replay_certificate(self.certificate, 1, 200_000_000)

    def test_non_primitive_presented_as_primitive(self):
        def scale(d):
            d["G"] = [[2 * x for x in row] for row in d["G"]]
            d["determinant"] *= 16
        self._reject(scale)
        self._reject(lambda d: d.__setitem__("scaling_c", 2))

    def test_manuscript_schema_tampering(self):
        self._reject(lambda d: d["identity"]["manuscript_schema"]["document"]["ambient_isometry"].__setitem__(0, [1, 0, 0, 0]))
        self._reject(lambda d: d["identity"]["manuscript_schema"].__setitem__("available", False))

    def test_tampered_enumeration_table(self):
        def drop_member(d):
            d["members"].pop(0)
        def relabel(d):
            d["members"][0]["level"] = 2 if d["members"][0]["level"] == 1 else 1
        def wrong_count(d):
            d["summary"]["level_2_examples"] += 1
        def forged_isometry(d):
            d["members"][1]["isometry_to_representative"][0][0] += 1
        for change in (drop_member, relabel, wrong_count, forged_isometry):
            doc = deepcopy(self.task_b)
            change(doc)
            with self.assertRaises(InvalidCertificate):
                V.replay_enumeration(doc, 2_000_000, 200_000_000)

    def test_tampered_lemma_record(self):
        doc = deepcopy(self.task_e)
        doc["cases"][0]["swapping_element"] = doc["cases"][0]["R"]
        with self.assertRaises(InvalidCertificate):
            V.replay_task_e(doc, 60, 2_000_000, 200_000_000)


if __name__ == "__main__":
    unittest.main()
