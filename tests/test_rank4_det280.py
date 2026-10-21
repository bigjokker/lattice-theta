"""Tests for the determinant-280 level-four certificate and its replay."""
from copy import deepcopy
import json
from pathlib import Path
import unittest

from exact_theta import InvalidCertificate
import verify_rank4_det280 as V

ROOT = Path(__file__).resolve().parent.parent


class Det280Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc = json.loads((ROOT / "data/rank4-level4/det280-level4.json").read_text(encoding="utf-8"))

    def test_certificate_replays(self):
        result = V.replay(self.doc, 5_000_000, 100_000_000)
        self.assertEqual(result["verdict"], "verified_level_4")
        self.assertEqual(result["determinant"], 280)
        self.assertEqual(result["full_group_order"], 2)
        self.assertEqual(result["stabilizer_epsilons"], [0, 0])
        self.assertEqual(result["cutoff"], 3072)

    def test_tampered_gram_is_rejected(self):
        bad = deepcopy(self.doc)
        bad["G"][3][3] += 2
        with self.assertRaises(InvalidCertificate):
            V.replay(bad, 5_000_000, 100_000_000)

    def test_tampered_characteristic_is_rejected(self):
        bad = deepcopy(self.doc)
        bad["b"] = [1, 1, 0, 0]
        with self.assertRaises(InvalidCertificate):
            V.replay(bad, 5_000_000, 100_000_000)

    def test_tampered_shell_record_is_rejected(self):
        bad = deepcopy(self.doc)
        bad["expected"]["vectors"] += 1
        with self.assertRaises(InvalidCertificate):
            V.replay(bad, 5_000_000, 100_000_000)


if __name__ == "__main__":
    unittest.main()
