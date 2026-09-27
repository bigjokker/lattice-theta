"""Uniform Dn formula, complete finite coverage, and witness/phase regressions."""

from copy import deepcopy
import gzip
from itertools import permutations, product
import json
from math import comb
from pathlib import Path
import unittest

from dn_ambient import ambient_histograms
from dn_exact import (basis, characteristic_tables, coefficient_tables, expand_record,
                      formula, gram, header, product_zero_count, records,
                      verify_ambient_witness, verify_records, zero_recipe)
from exact_theta import EnumerationLimit, InvalidCertificate, multiply, transpose, verify
from exceptional_batch import ldl_histograms
from verify_dn import verify_stream

ROOT = Path(__file__).resolve().parent.parent


class DnTests(unittest.TestCase):
    def test_closed_formula_matches_unsimplified_count(self):
        for rank in range(2, 61):
            balanced = comb(rank, rank // 2) // 2 if rank % 2 == 0 else 0
            self.assertEqual(formula(rank), product_zero_count(rank) + balanced)
        self.assertEqual([formula(n) for n in range(2, 11)],
                         [1, 0, 9, 60, 400, 2100, 10241, 46620, 204756])

    def test_characteristic_parametrizations_are_bijective(self):
        for rank in range(2, 11):
            primals, duals = characteristic_tables(rank)
            self.assertEqual(len({(S, e) for _, S, e in primals}), 1 << rank)
            self.assertTrue(all(S.bit_count() % 2 == 0 for _, S, _ in primals))
            full = (1 << rank) - 1
            self.assertEqual(len({(h, min(W, W ^ full)) for _, h, W in duals}), 1 << rank)

    def test_all_saved_streams_reverify_and_match_primary(self):
        total = 0
        for rank in range(2, 11):
            with self.subTest(rank=rank):
                result = verify_stream(ROOT / "data" / "dn" / f"d{rank}.jsonl.gz")
                primary = json.loads((ROOT / "reports" / f"dn-d{rank}-ldl.json").read_text(encoding="utf-8"))
                for field in ("records_sha256", "shell_histogram_sha256", "certificate_stream_sha256",
                              "even_vanishing", "even_nonvanishing", "even_characteristics"):
                    self.assertEqual(result[field], primary[field])
                self.assertEqual(result["verdict"], "classification_verified")
                self.assertEqual(result["even_vanishing"], formula(rank))
                self.assertEqual(result["unresolved"], 0)
                total += result["even_characteristics"]
        self.assertEqual(total, 700070)

    def test_ambient_and_ldl_histograms_are_identical(self):
        for rank in range(2, 7):
            bound = header(rank)["norm4_bound"]
            ambient, _ = ambient_histograms(rank, bound)
            primary, _ = ldl_histograms(gram(rank), bound)
            self.assertEqual(primary, ambient)

    def test_every_small_zero_expands_to_original_verifier(self):
        for rank in (2, 4, 5):
            doc = header(rank)
            coefficients, norms, _, _ = coefficient_tables(doc, "ambient")
            zeros = 0
            for row in records(doc, coefficients, norms):
                if row[2] == 0:
                    continue
                result = verify(expand_record(doc, row))
                self.assertEqual(result["verdict"], "vanishes_by_symmetry")
                self.assertEqual(result["symmetry"]["witness_order"], 2)
                zeros += 1
            self.assertEqual(zeros, formula(rank))

    def test_quarter_shift_even_classes_never_vanish(self):
        for rank in range(2, 7):
            document = header(rank)
            coefficients, _, _, _ = coefficient_tables(document, "ambient")
            tables = characteristic_tables(rank)
            for am, (_, S, _) in enumerate(tables[0]):
                for bm, (_, h, _) in enumerate(tables[1]):
                    if h != 1 or (am & bm).bit_count() % 2:
                        continue
                    self.assertIsNone(zero_recipe(am, bm, rank, tables))
                    self.assertEqual(abs(coefficients[(am, S.bit_count())][bm]), 2 ** (S.bit_count() // 2))

    def test_integer_sector_nonzero_coefficients(self):
        for rank in range(2, 7):
            document = header(rank)
            coefficients, _, _, _ = coefficient_tables(document, "ambient")
            tables = characteristic_tables(rank)
            for am, (_, S, e) in enumerate(tables[0]):
                for bm, (_, h, W) in enumerate(tables[1]):
                    if h or (am & bm).bit_count() % 2 or zero_recipe(am, bm, rank, tables) is not None:
                        continue
                    if S:
                        N, expected = S.bit_count(), 2 ** (S.bit_count() - 1)
                    elif e == 0:
                        N, expected = 0, 1
                    else:
                        N, expected = 4, 2 * abs(rank - 2 * W.bit_count())
                    self.assertEqual(abs(coefficients[(am, N)][bm]), expected)

    def test_missing_duplicate_corrupt_and_extra_records_rejected(self):
        document = header(4)
        coefficients, norms, histograms, stats = coefficient_tables(document, "ambient")
        rows = list(records(document, coefficients, norms))
        corrupt_coefficient = deepcopy(rows)
        corrupt_coefficient[0][4] += 1
        corrupt_witness = deepcopy(rows)
        index = next(i for i, row in enumerate(rows) if row[2] == 1)
        corrupt_witness[index][3] = corrupt_witness[index][4]
        mutations = (rows[:-1], [rows[0], *rows], [*rows, rows[0]], corrupt_coefficient, corrupt_witness)
        for changed in mutations:
            with self.assertRaises(InvalidCertificate):
                verify_records(document, changed, coefficients, norms, histograms, stats, "ambient")

    def test_resource_limit_does_not_publish_a_classification(self):
        with self.assertRaises(EnumerationLimit):
            coefficient_tables(header(4), "ldl", max_nodes=1)

    def test_d2_and_d3_isometries(self):
        self.assertEqual(gram(2), [[2, 0], [0, 2]])
        U = [[0, 1, 0], [1, 0, 0], [0, 0, 1]]
        self.assertEqual(multiply(multiply(transpose(U), gram(3)), U),
                         [[2, -1, 0], [-1, 2, -1], [0, -1, 2]])

    def test_d4_signed_permutation_witness_count(self):
        count = 0
        for permutation in permutations(range(4)):
            for signs in product((-1, 1), repeat=4):
                try:
                    verify_ambient_witness([2, 0, 0, 0], [0, 0, 2, 2], permutation, signs)
                except InvalidCertificate:
                    continue
                count += 1
        self.assertEqual(count, 64)
        with self.assertRaises(InvalidCertificate):
            verify_ambient_witness([2, 0, 0, 0], [0, 0, 2, 2], [0, 1, 2, 3], [-1, 1, -1, 1])

    def test_standalone_certificate_examples(self):
        for path in (ROOT / "data" / "dn").glob("*.json"):
            result = verify(json.loads(path.read_text(encoding="utf-8")))
            self.assertIn(result["verdict"], ("vanishes_by_symmetry", "does_not_vanish"))


if __name__ == "__main__":
    unittest.main()
