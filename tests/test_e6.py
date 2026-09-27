"""E6 coverage, corrupted-evidence rejection, and independent enumeration checks."""

from copy import deepcopy
import json
from pathlib import Path
import unittest

from exact_theta import InvalidCertificate, complete_shells, verify
from verify_even_pack import verify_pack
from verify_e6_ambient import (DOUBLED_ROOTS, ambient_shells, gram_from_roots,
                               recover_coordinates, verify_ambient)

PACK = json.loads((Path(__file__).resolve().parent.parent / "data" / "e6-even.json").read_text(encoding="utf-8"))


class E6Tests(unittest.TestCase):
    def test_complete_primary_verification(self):
        result = verify_pack(PACK)
        self.assertEqual(result["verdict"], "parity_converse_verified")
        self.assertEqual(result["determinant"], 3)
        self.assertEqual(result["even_characteristics"], 2080)
        self.assertEqual(result["odd_characteristics"], 2016)
        self.assertEqual(result["nonvanishing_verified"], 2080)
        self.assertEqual(result["maximum_witness_norm"], "1")
        self.assertEqual(result["unresolved"], [])
        self.assertEqual(result["witness_norm4_counts"], {"0": 64, "2": 1152, "4": 864})

    def test_complete_independent_verification(self):
        result = verify_ambient(PACK)
        self.assertEqual(result["even_characteristics_verified"], 2080)
        self.assertEqual(result["odd_characteristics_checked_through_bound"], 2016)
        self.assertEqual(result["accepted_lattice_vectors"], 343)
        self.assertEqual(result["minimal_coset_shell_counts"], [
            {"norm4": 0, "vector_count": 1, "cosets": 1},
            {"norm4": 2, "vector_count": 2, "cosets": 36},
            {"norm4": 4, "vector_count": 10, "cosets": 27}])

    def test_missing_and_duplicate_classes_rejected(self):
        missing = deepcopy(PACK)
        missing["certificates"].pop()
        duplicate = deepcopy(PACK)
        duplicate["certificates"][-1] = duplicate["certificates"][0]
        for pack in (missing, duplicate):
            with self.subTest(pack="missing" if pack is missing else "duplicate"):
                with self.assertRaises(InvalidCertificate):
                    verify_pack(pack)
                with self.assertRaises(ValueError):
                    verify_ambient(pack)

    def test_bad_characteristic_or_evidence_rejected(self):
        mutations = []
        odd = deepcopy(PACK)
        odd["certificates"][0]["a"] = [1, 0, 0, 0, 0, 0]
        odd["certificates"][0]["b"] = [1, 0, 0, 0, 0, 0]
        mutations.append(odd)
        nonbinary = deepcopy(PACK)
        nonbinary["certificates"][0]["a"][0] = 2
        mutations.append(nonbinary)
        wrong = deepcopy(PACK)
        wrong["certificates"][0]["evidence"]["coefficient"] = 2
        mutations.append(wrong)
        zero = deepcopy(PACK)
        zero["certificates"][0]["evidence"]["coefficient"] = 0
        mutations.append(zero)
        boolean = deepcopy(PACK)
        boolean["certificates"][0]["a"][0] = False
        mutations.append(boolean)
        for i, pack in enumerate(mutations):
            with self.subTest(mutation=i):
                with self.assertRaises(InvalidCertificate):
                    verify_pack(pack)
                with self.assertRaises(ValueError):
                    verify_ambient(pack)

    def test_independent_verifier_requires_e6_input(self):
        changed = deepcopy(PACK)
        changed["G"][0][0] = 4
        with self.assertRaisesRegex(ValueError, "not the specified E6"):
            verify_ambient(changed)
        self.assertEqual(gram_from_roots(), PACK["G"])

    def test_resource_limit_cannot_prove_the_classification(self):
        # A three-class rank-one pack suffices to check the adapter's limit handling.
        pack = {"schema": "work7-even-theta-pack-v1", "G": [[2]],
                "representatives": "binary-primal-and-dual", "certificates": [
                    {"a": [0], "b": [0], "evidence": {"kind": "shell", "norm4": 0, "coefficient": 1}},
                    {"a": [0], "b": [1], "evidence": {"kind": "shell", "norm4": 0, "coefficient": 1}},
                    {"a": [1], "b": [0], "evidence": {"kind": "shell", "norm4": 2, "coefficient": 2}}]}
        result = verify_pack(pack, max_nodes=1)
        self.assertTrue(result["coverage_complete"])
        self.assertEqual(result["verdict"], "unresolved")
        self.assertEqual(result["nonvanishing_verified"], 0)
        self.assertEqual(len(result["unresolved"]), 3)

    def test_all_coset_shell_counts_agree_between_algorithms(self):
        histograms, _ = ambient_shells(4)
        by_coset = {tuple(entry["a"]) for entry in PACK["certificates"]}
        self.assertEqual(len(by_coset), 64)
        for a in by_coset:
            result = complete_shells(PACK["G"], list(a), [0] * 6, 4)
            primary = {s["norm4"]: s["vector_count"] for s in result["shells"]}
            ambient = {N: sum(hist.values()) for (aa, N), hist in histograms.items() if aa == a}
            self.assertEqual(primary, ambient)

    def test_inverse_ambient_coordinates(self):
        # Includes negative and larger coordinates beyond the verification sphere.
        for y in ([1, -3, 2, 5, -4, 7], [-2, 1, 0, -1, 3, -5], [0] * 6):
            z = [sum(c * v for c, v in zip(row, y)) for row in DOUBLED_ROOTS]
            self.assertEqual(recover_coordinates(z[:5], z[7]), y)
        self.assertIsNone(recover_coordinates([0, 1, 0, 0, 0], 0))

    def test_single_entry_is_a_standard_certificate(self):
        entry = next(entry for entry in PACK["certificates"] if entry["evidence"]["norm4"] == 4)
        doc = {"schema": "work7-theta-v1", "G": PACK["G"], **entry}
        self.assertEqual(verify(doc)["verdict"], "does_not_vanish")


if __name__ == "__main__":
    unittest.main()
