"""Tests for the determinant census through 279 (release 1.2.0) and its replay."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import unittest

from exact_theta import InvalidCertificate
import rank4_gap_census as C
import verify_rank4_parity as V

ROOT = Path(__file__).resolve().parent.parent
CHUNKS = ROOT / "reports/rank4-gap-census"
SMALL = ("even-div4-65-279", 68)


def levels(summary):
    return {str(k): v for k, v in summary["zero_levels"].items()}


class GapCensusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc = C.load(C.chunk_path(CHUNKS, *SMALL))
        V.FAST_BOX[0] = True

    def test_new_regions_and_prior_lists_partition_the_gap(self):
        gap = set(range(25, 280))
        task_e = {d for d in gap if d <= 168 and d % 4}
        odd_new = set(C.admitted("odd-notdiv4-169-279")) | set(C.admitted("odd-div4-25-279"))
        self.assertEqual(task_e | odd_new, gap)
        self.assertFalse(task_e & odd_new)
        task_b = {d for d in gap if d <= 201 and d % 2}
        task_c = {d for d in gap if d <= 64 and d % 4 == 0}
        even_new = set(C.admitted("even-odd-202-279")) | set(C.admitted("even-div4-65-279"))
        # the lists' determinant filters cover every odd determinant and every multiple of
        # four; an even lattice never has det = 2 mod 4 (checked by the generator in `summary`)
        self.assertEqual(task_b | task_c | even_new, {d for d in gap if d % 2 or d % 4 == 0})
        self.assertFalse((task_b | task_c) & even_new)
        self.assertEqual(sum(len(C.admitted(r)) for r in C.CLOSURE_REGIONS), 239)

    def test_preserved_chunk_replays_and_rebuilds(self):
        result = V.replay_decided(self.doc, 5_000_000, 50_000_000)
        self.assertEqual(result["verdict"], "verified")
        self.assertEqual((result["classes"], result["indecomposable_classes"], result["P_fails"]), (7, 5, 3))
        with tempfile.TemporaryDirectory() as tmp:
            row = C.build_chunk((*SMALL, Path(tmp)))
            fresh = C.load(C.chunk_path(tmp, *SMALL))
        self.assertEqual(fresh["classes_doc"]["candidate_bases"], self.doc["classes_doc"]["candidate_bases"])
        self.assertEqual([r["G"] for r in fresh["lattices"]], [r["G"] for r in self.doc["lattices"]])
        self.assertEqual(levels(fresh["summary"]), levels(self.doc["summary"]))
        self.assertEqual(row["even_zeros"], 16)
        self.assertEqual(row["zeros_without_witness"], 0)

    def test_tampered_group_or_summary_is_rejected(self):
        bad = deepcopy(self.doc)
        bad["lattices"][0]["aut_order"] += 2
        with self.assertRaises(InvalidCertificate):
            V.replay_decided(bad, 5_000_000, 50_000_000)
        bad = deepcopy(self.doc)
        bad["summary"]["P_fails"] -= 1
        with self.assertRaises(InvalidCertificate):
            V.replay_decided(bad, 5_000_000, 50_000_000)

    def test_preserved_summaries_state_the_closure_and_the_unique_failure_at_280(self):
        closure = json.loads((CHUNKS / "closure-25-279.json").read_text(encoding="utf-8"))
        self.assertTrue(closure["symmetry_converse_holds_for_25_to_279"])
        totals = closure["new_totals"]
        self.assertEqual((totals["classes"], totals["indecomposable_classes"]), (14987, 6239))
        self.assertEqual((totals["even_zeros"], totals["zeros_without_witness"], totals["unresolved_characteristics"]),
                         (5330, 0, 0))
        for region in C.CLOSURE_REGIONS:
            s = json.loads((CHUNKS / f"{region}-summary.json").read_text(encoding="utf-8"))
            r = json.loads((CHUNKS / f"{region}-replay.json").read_text(encoding="utf-8"))
            self.assertTrue(s["complete"] and r["all_verified"])
            self.assertEqual([c["sha256"] for c in s["chunks"]], [c["sha256"] for c in r["chunks"]])
        odd = json.loads((CHUNKS / "odd-div4-280-summary.json").read_text(encoding="utf-8"))
        even = json.loads((CHUNKS / "even-div4-280-summary.json").read_text(encoding="utf-8"))
        self.assertEqual(odd["totals"]["zero_levels"], {"1": 162, "4": 1})
        self.assertEqual(odd["totals"]["zeros_without_witness"], 1)
        self.assertEqual(even["totals"]["zero_levels"], {"1": 24})
        self.assertEqual((odd["totals"]["classes"], even["totals"]["classes"]), (243, 27))


if __name__ == "__main__":
    unittest.main()
