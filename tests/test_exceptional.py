"""Complete E7/E8 classifications, algorithm agreement and corruption checks."""

from copy import deepcopy
import json
from pathlib import Path
import unittest

from exact_theta import EnumerationLimit, InvalidCertificate, identity, verify
from exceptional_ambient import ambient_histograms, direct_coefficients, gram
from exceptional_batch import (export_certificate, histogram_hash, ldl_histograms,
                                path_word, root_gram, structural_checks,
                                symmetry_forest, validate_pack, walsh)
from verify_exceptional import verify_classification

ROOT = Path(__file__).resolve().parent.parent
PACKS = {rank: json.loads((ROOT / "data" / f"e{rank}-classification.json").read_text(encoding="utf-8"))
         for rank in (7, 8)}


class ExceptionalTests(unittest.TestCase):
    def test_complete_classifications_both_backends(self):
        for rank, zeros, nonzeros, odds in ((7, 1260, 6996, 8128), (8, 9450, 23446, 32640)):
            results = []
            for backend in ("ldl", "ambient"):
                with self.subTest(rank=rank, backend=backend):
                    result = verify_classification(PACKS[rank], backend)
                    self.assertEqual(result["verdict"], "classification_verified")
                    self.assertEqual(result["even_vanishing_by_symmetry"], zeros)
                    self.assertEqual(result["even_nonvanishing_by_complete_shell"], nonzeros)
                    self.assertEqual(result["odd_characteristics_vanish_by_parity"], odds)
                    self.assertEqual(result["unresolved"], 0)
                    self.assertEqual(result["seed_witnesses"], 1)
                    self.assertTrue(result["minimal_shell_detects_all_even_nonvanishing"])
                    results.append(result)
            self.assertEqual(results[0]["shell_histogram_sha256"], results[1]["shell_histogram_sha256"])
            self.assertEqual(results[0]["witness_norm4_counts"], results[1]["witness_norm4_counts"])

    def test_walsh_matches_direct_sign_sums(self):
        # Actual E8 histograms, including negative x parity and cancellations.
        histograms, _ = ambient_histograms(8, 8)
        direct = direct_coefficients(histograms, 8)
        for key, histogram in histograms.items():
            self.assertEqual(walsh(histogram, 8), direct[key])

    def test_e6_batch_compatibility(self):
        from verify_e6_ambient import ambient_shells
        old, _ = ambient_shells(4)
        new, _ = ambient_histograms(6, 4)
        canonical = {(sum(x << i for i, x in enumerate(a)), N): histogram
                     for (a, N), histogram in old.items()}
        self.assertEqual(canonical, new)
        self.assertEqual(root_gram(6), gram(6))

    def test_transport_expands_to_original_symmetry_certificate(self):
        for rank, pack in PACKS.items():
            deepest = max(range(len(pack["zero_nodes"])), key=lambda i: len(path_word(pack["zero_nodes"], i)))
            node_indices = (0, len(pack["zero_nodes"]) - 1, deepest)
            for ni in node_indices:
                entry = next(e for e in pack["classes"] if e["evidence"] == {"kind": "symmetry", "node": ni})
                with self.subTest(rank=rank, node=ni):
                    result = verify(export_certificate(pack, entry))
                    self.assertEqual(result["verdict"], "vanishes_by_symmetry")
                    self.assertEqual(result["symmetry"]["epsilon"], 1)

    def test_standalone_examples(self):
        for rank in (7, 8):
            for label, verdict in (("zero", "vanishes_by_symmetry"), ("nonzero", "does_not_vanish")):
                doc = json.loads((ROOT / "data" / f"e{rank}-{label}.json").read_text(encoding="utf-8"))
                self.assertEqual(verify(doc)["verdict"], verdict)

    def test_bad_coverage_rejected(self):
        missing = deepcopy(PACKS[7])
        missing["classes"].pop()
        duplicate = deepcopy(PACKS[7])
        duplicate["classes"][-1] = duplicate["classes"][0]
        odd = deepcopy(PACKS[7])
        odd["classes"][0].update(a_mask=1, b_mask=1)
        for pack in (missing, duplicate, odd):
            with self.assertRaises(InvalidCertificate):
                validate_pack(pack)

    def test_bad_symmetry_proof_rejected(self):
        no_sign = deepcopy(PACKS[7])
        no_sign["seeds"][0]["T"] = identity(7)
        no_sign["seeds"][0]["word"] = []
        cyclic = deepcopy(PACKS[7])
        cyclic["zero_nodes"][1]["parent"] = 1
        wrong_transport = deepcopy(PACKS[7])
        wrong_transport["zero_nodes"][1]["b_mask"] ^= 1
        nonisometry = deepcopy(PACKS[7])
        nonisometry["generators"][0][0][0] += 1
        for pack in (no_sign, cyclic, wrong_transport, nonisometry):
            with self.assertRaises(InvalidCertificate):
                validate_pack(pack)

    def test_corrupt_shell_coefficient_rejected(self):
        pack = deepcopy(PACKS[7])
        pack["classes"][0]["evidence"]["coefficient"] = 2
        with self.assertRaisesRegex(InvalidCertificate, "coefficient mismatch"):
            verify_classification(pack, "ambient")

    def test_small_search_range_or_limit_does_not_prove_nonvanishing(self):
        pack = deepcopy(PACKS[7])
        pack["norm4_bound"] = 0
        with self.assertRaises(InvalidCertificate):
            validate_pack(pack)
        with self.assertRaises(EnumerationLimit):
            ldl_histograms(root_gram(7), 8, max_nodes=1)

    def test_no_symmetry_proof_means_unresolved(self):
        # Even a=0,b=0 never vanishes; an artificial candidate must not get a verdict.
        with self.assertRaisesRegex(InvalidCertificate, "no symmetry proof; unresolved"):
            symmetry_forest(root_gram(6), {(0, 0)})

    def test_structural_characterizations(self):
        e7 = structural_checks(PACKS[7])
        self.assertEqual(e7["even_zeros_outside_primal_lattice"], 1260)
        e8 = structural_checks(PACKS[8])
        self.assertEqual(e8["totally_singular_planes"], 1575)
        self.assertEqual(e8["singular_vectors_including_zero"], 136)
        self.assertTrue(e8["quadratic_criterion_verified"])
        self.assertEqual(len(e8["totally_singular_four_space_basis_masks"]), 4)
        self.assertTrue(e8["one_zero_orbit_under_supplied_generators"])


if __name__ == "__main__":
    unittest.main()
