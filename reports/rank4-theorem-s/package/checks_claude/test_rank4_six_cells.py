"""Tests for rank4_secondary_fan.py and rank4_six_cells.py (exact bounds, markings, dispositions, D = Z E)."""
from fractions import Fraction as Fr
from itertools import product
import json
import random
import unittest

import rank4_secondary_fan as SF
import rank4_six_cells as C

A4STAR = [[4, -1, -1, -1], [-1, 4, -1, -1], [-1, -1, 4, -1], [-1, -1, -1, 4]]


class ExactFan(unittest.TestCase):
    def test_ellipsoid_enumeration_is_exact(self):
        rng = random.Random(3)
        for Q in (A4STAR, [[2, 1, 0, 0], [1, 9, 3, 0], [0, 3, 5, 2], [0, 0, 2, 40]]):
            for _ in range(8):
                c = [Fr(rng.randint(-15, 15), rng.randint(1, 6)) for _ in range(4)]
                r2 = Fr(rng.randint(1, 40), rng.randint(1, 4))
                got = set(SF.lattice_points_in_ellipsoid(Q, c, r2))
                base = [round(x) for x in c]
                brute = {tuple(b + d for b, d in zip(base, x)) for x in product(range(-8, 9), repeat=4)}
                self.assertEqual(got, {x for x in brute if SF.qform(Q, [Fr(a) - b for a, b in zip(x, c)]) <= r2})

    def test_principal_cone(self):
        cone = SF.secondary_cone(A4STAR)
        self.assertEqual((len(cone["facets"]), len(cone["rays"]), cone["simplices_at_0"]), (10, 10, 120))

    def test_saved_fan_has_coverage_certificate(self):
        fan = json.load(open(SF.OUT / "fan-exact.json"))
        self.assertEqual(fan["orbits"], 3)
        self.assertTrue(C.coverage_ok(fan))


class Classifier(unittest.TestCase):
    def test_markings(self):
        self.assertEqual(len(C.TAUS), 210)
        for t in C.TAUS:
            self.assertEqual(len(C.pairs_of(t)), 6)

    def test_full_principal_cone_dispositions(self):
        fan = json.load(open(SF.OUT / "fan-exact.json"))
        r = C.process_face((0, fan["cones"][0]["rays"], tuple(range(10))))
        self.assertEqual(sum(r["dispositions"].values()), 210)
        self.assertTrue(all(c["certified"] and c["D_equals_ZE"] for c in r["cells"]))

    def test_dze_rejects_a_wrong_T(self):
        fan = json.load(open(SF.OUT / "fan-exact.json"))
        rays = fan["cones"][0]["rays"]
        r = C.process_face((0, rays, tuple(range(10))))
        cell = r["cells"][0]
        I = [[int(i == j) for j in range(4)] for i in range(4)]
        swap = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
        self.assertTrue(C.dze_holds(rays, cell["equalities"], cell["T"]))
        self.assertTrue(C.dze_holds(rays, cell["equalities"], I))            # identity preserves everything
        self.assertFalse(C.dze_holds(rays, cell["equalities"], swap))        # a coordinate swap does not


if __name__ == "__main__":
    unittest.main()
