"""Adversarial tests for the rank-four level-4 search and its independent replay."""
from copy import deepcopy
import unittest

from exact_theta import EnumerationLimit, InvalidCertificate
import rank4_level4_search as S
import verify_rank4_level4 as V


class Level4Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.task_a = S.task_a(prefix=400, filter_bound=200, names=["A1"])
        report = cls.task_a["pairs"][0]
        pair = {"name": report["name"], "source": report["source"], "G1": report["G1"], "G2": report["G2"]}
        cls.certificate = S.certify_configuration(pair, report["klein_configurations"][0])
        cls.bijection = S.cs_bijection_certificate()

    # --- genuine material ---------------------------------------------------

    def test_task_a_finds_the_configuration_and_replays(self):
        report = self.task_a["pairs"][0]
        self.assertEqual(report["verdict"], "klein_configuration_found")
        self.assertEqual(report["sub_matches"], 1)
        self.assertEqual(report["over_matches"], 1)
        self.assertEqual(report["direct_isometries_L1_L2"], 0)
        rows = V.replay_task_a(self.task_a, 1_000_000, 100_000_000)
        self.assertEqual(rows[0]["isometric_cells"], 2)

    def test_certificate_is_level_four_with_two_identity_routes(self):
        c = self.certificate
        self.assertEqual(c["determinant"], 1729)
        self.assertEqual(c["level"], 4)
        self.assertEqual(len(c["automorphisms"]), 2)
        self.assertEqual(c["kernel_automorphisms"]["elements_mapping_v1_to_v2"], 0)
        self.assertTrue(c["identity"]["bijection"]["applies"])
        self.assertTrue(c["identity"]["modular_cutoff"]["attempted"])
        self.assertEqual(c["identity"]["modular_cutoff"]["ldl_backend"]["evidence"]["verdict"],
                         "proved_zero_by_modular_cutoff")
        self.assertEqual(V.replay_bijection(self.bijection)["verdict"], "verified")
        result = V.replay_certificate(c, 1_000_000, 100_000_000)
        self.assertEqual(result["verdict"], "verified_level_4")
        self.assertEqual(len(result["identity_routes"]), 2)
        self.assertFalse(result["members_isometric"])

    def test_two_enumeration_backends_agree_with_the_box(self):
        c = self.certificate
        G, a, b = c["G"], c["a"], c["b"]
        shells = S.signed_shells(G, a, b, 600)["shells"]
        from exact_theta import complete_shells
        self.assertEqual(shells, complete_shells(G, a, b, 600)["shells"])
        self.assertEqual(shells, V.box_signed_shells(G, a, b, 600)[0])
        self.assertTrue(all(s["signed_coefficient"] == 0 for s in shells))

    # --- tampering ----------------------------------------------------------

    def _reject(self, change, exception=InvalidCertificate):
        doc = deepcopy(self.certificate)
        change(doc)
        with self.assertRaises(exception):
            V.replay_certificate(doc, 1_000_000, 100_000_000)

    def test_altered_gram(self):
        self._reject(lambda d: d["G"][0].__setitem__(0, 4))
        self._reject(lambda d: d["G0"][1].__setitem__(1, d["G0"][1][1] + 1))

    def test_altered_characteristic(self):
        self._reject(lambda d: d.__setitem__("a", [1, 1, 0, 1]))
        self._reject(lambda d: d.__setitem__("b", [1, 1, 1, 0]))
        self._reject(lambda d: d.__setitem__("construction_a", [1, 0, 0, 0]))

    def test_removed_group_element_and_wrong_sign(self):
        self._reject(lambda d: d["automorphisms"].pop())
        self._reject(lambda d: d["stabilizer"][0].__setitem__("epsilon", 1))
        self._reject(lambda d: d["kernel_reduced"]["automorphisms"].pop())
        self._reject(lambda d: d["kernel_automorphisms"].__setitem__("elements_mapping_v1_to_v2", 1))

    def test_truncated_or_altered_enumeration(self):
        def truncate(d):
            d["identity"]["modular_cutoff"]["ldl_backend"]["evidence"]["enumeration"]["shells"].pop()
        self._reject(truncate)

        def shorter_bound(d):
            d["identity"]["modular_cutoff"]["ldl_backend"]["run"]["norm4_bound"] -= 1
        self._reject(shorter_bound)

        def nonzero_coefficient(d):
            d["identity"]["modular_cutoff"]["ldl_backend"]["evidence"]["enumeration"]["shells"][0]["signed_coefficient"] = 2
        self._reject(nonzero_coefficient)

        def incomplete(d):
            d["identity"]["modular_cutoff"]["ldl_backend"]["evidence"]["enumeration_complete"] = False
        self._reject(incomplete)

    def test_interrupted_replay_cannot_certify(self):
        with self.assertRaises(EnumerationLimit):
            V.replay_certificate(self.certificate, 1_000_000, 10)
        with self.assertRaises(EnumerationLimit):
            V.replay_certificate(self.certificate, 1, 100_000_000)

    def test_tampered_links_and_false_non_isometry(self):
        self._reject(lambda d: d["links"]["O_v1"]["reduced_isometry"][0].__setitem__(0, 0))
        self._reject(lambda d: d["links"].__setitem__("O_v2", deepcopy(d["links"]["O_v1"])))

        def same_member(d):
            d["members"]["L2"] = deepcopy(d["members"]["L1"])
            d["members"]["L2_reduced"] = deepcopy(d["members"]["L1_reduced"])
        self._reject(same_member)
        self._reject(lambda d: d["non_isometry"]["ldl_sphere"].__setitem__("isometries_found", 1))

    def test_tampered_coset_lattice_and_kernel(self):
        self._reject(lambda d: d["coset_lattices"]["O_v1"]["doubled_basis"][0].__setitem__(0, 7))
        self._reject(lambda d: d["kernel_basis"][0].__setitem__(0, 4))
        self._reject(lambda d: d.__setitem__("complement_t", [0, 0, 0, 0]))

    def test_tampered_bijection(self):
        for change in (lambda d: d["cosets"][3]["image"].__setitem__(0, d["cosets"][3]["image"][0] + 3),
                       lambda d: d["cosets"][3].__setitem__("flip", (d["cosets"][3]["flip"] + 1) % 4),
                       lambda d: d["cosets"].pop(),
                       lambda d: d["M"][0].__setitem__(0, 3),
                       lambda d: d["V_plus"][0].__setitem__(0, 2)):
            doc = deepcopy(self.bijection)
            change(doc)
            with self.assertRaises(InvalidCertificate):
                V.replay_bijection(doc)

    def test_tampered_task_a_tables(self):
        def flip_cell(d):
            cell = d["pairs"][0]["sub_table"]["cells"][0]
            cell["first_theta_difference"] = None if cell["first_theta_difference"] is not None else 5

        def drop_match(d):
            d["pairs"][0]["sub_table"]["isometric_pairs"].pop()

        def wrong_gram(d):
            d["pairs"][0]["G2"][0][0] += 2

        def wrong_klein(d):
            d["pairs"][0]["klein_configurations"][0]["a"] = [1, 0, 0, 0]

        for change in (flip_cell, drop_match, wrong_gram, wrong_klein):
            doc = deepcopy(self.task_a)
            change(doc)
            with self.assertRaises(InvalidCertificate):
                V.replay_task_a(doc, 1_000_000, 100_000_000)


if __name__ == "__main__":
    unittest.main()
