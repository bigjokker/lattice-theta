"""Tests for release 1.3.0: Theorem 5 (minus type), the stage-two reductions and the plus-type census."""
from copy import deepcopy
import json
from pathlib import Path
import unittest

import rank4_minus_type as MT
import rank4_plus_type_census as PT
import rank4_stage2_checks as S2
import rank4_stage2_review_check as RC
import verify_rank4_parity as V
from rank4_feasibility import isometries
from rank4_level4_search import gram_of, kernel_basis

ROOT = Path(__file__).resolve().parent.parent
TASK_B = ROOT / "reports/rank4-parity/task-b-det201.json"
TASK_C = ROOT / "reports/rank4-parity/task-c-det64.json"
PLUS = ROOT / "reports/rank4-plus-type"


def preserved(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


class StageTwoTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.task_b = json.loads(MT._source_bytes(TASK_B))
        V.FAST_BOX[0] = True

    def test_order_four_symplectic_maps_fix_one_line(self):
        result = MT.finite_group_check()
        self.assertEqual(result["group_orders"], {"Sp4_F2": 720, "O4_plus_F2": 72, "O4_minus_F2": 120})
        order_four = {k: v for k, v in result["tables"]["Sp4_F2"].items() if k.startswith("order=4,")}
        self.assertEqual(order_four, {"order=4,fixed_vectors=2": 180})

    def test_theorem5_audit_of_the_odd_determinant_list(self):
        audit = MT.audit_census(TASK_B)
        kept = preserved("reports/rank4-minus-type/odd-det201-symmetry-audit.json")
        self.assertEqual(audit["sha256"], kept["sha256"])
        for key in ("classes", "indecomposable", "minus_indecomposable", "minus_indecomposable_symmetry_zeros",
                    "indecomposable_symmetry_zeros"):
            self.assertEqual(audit[key], kept[key], key)
        self.assertEqual((audit["minus_indecomposable"], audit["minus_indecomposable_symmetry_zeros"]), (131, 0))

    def test_plane_swap_transform_and_decomposable_control(self):
        self.assertEqual(S2.plane_transform_checks(self.task_b),
                         {"lattices": 300, "characteristic_swaps_and_dual_Q_checks": 31500})
        control = S2.decomposable_control()
        self.assertEqual((control["determinant"], control["Q_on_nonzero_plane_vectors"], control["witness_epsilon"]),
                         (21, [0, 1, 1], 1))

    def test_mixed_order_four_reflections_and_mapping_to_swapping(self):
        kernels = json.loads(MT._source_bytes(TASK_C))
        control = S2.decomposable_control()
        k = gram_of(control["G"], kernel_basis(control["b"]))
        data = {"lattices": kernels["lattices"] + [{"G": k, "automorphisms": isometries(k, k)}]}
        counts = S2.mixed_and_mapping_checks(data)
        self.assertEqual(counts, preserved("reports/rank4-stage2/strategy-checks.json")["mixed_kernel_and_mapping"])

    def test_independent_review_check_on_the_two_controls(self):
        rows = {r["index"]: r for r in preserved("reports/rank4-stage2/review-independent-check.json")["rows"]}
        lattices = {e["index"]: e["G"] for e in self.task_b["lattices"]}
        for index in (160, 293):  # det 21 decomposable (nonsquare), det 121 Gaussian (square)
            fresh = RC.analyse((index, lattices[index]))
            self.assertEqual(fresh, rows[index])
            self.assertGreater(fresh["maps"], 0)
            self.assertEqual(fresh["map_no_swap"], 0)

    def test_plus_type_chunk_replays_and_matches_the_h_planes(self):
        doc = PT.load(PLUS / PT.REGION / "det361.json.gz")
        result = V.replay_decided(doc, 5_000_000, 50_000_000)
        self.assertEqual(result["verdict"], "verified")
        self.assertEqual(result["P_fails"], 1)
        cmp = PT.compare(doc)
        self.assertEqual((cmp["zeros"], cmp["predicted"], cmp["unexplained_zeros"], cmp["predicted_but_nonzero"]),
                         (6, 6, 0, 0))
        bad = deepcopy(doc)
        row = next(r for r in bad["lattices"] if not r["decomposable"] and not r["has_square_root_of_minus_identity"])
        row["parity"]["decisions"][0]["verdict"] = "zero_by_symmetry"
        self.assertEqual(PT.compare(bad)["unexplained_zeros"], 1)

    def test_preserved_plus_type_outputs_state_the_combined_census(self):
        summary = preserved(f"reports/rank4-plus-type/{PT.REGION}-summary.json")
        replay = preserved(f"reports/rank4-plus-type/{PT.REGION}-replay.json")
        self.assertTrue(summary["complete"] and replay["all_verified"])
        self.assertEqual([c["sha256"] for c in summary["chunks"]], [c["sha256"] for c in replay["chunks"]])
        self.assertEqual(len(summary["chunks"]), 28)
        totals = summary["totals"]
        self.assertEqual((totals["classes"], totals["indecomposable_classes"], totals["even_zeros"],
                          totals["zeros_without_witness"], totals["unresolved_characteristics"]), (672, 623, 6, 0, 0))
        sources = preserved("reports/rank4-plus-type/planes.json")["sources"]
        self.assertEqual(sum(s["classes"] for s in sources), 300 + 269 + 787 + 672)
        self.assertTrue(all(s["unexplained_zeros"] == s["predicted_but_nonzero"] == s["unresolved"] == 0 for s in sources))


if __name__ == "__main__":
    unittest.main()
