"""Tests for rank4_six_minima.py (common lattice K, six minima, symmetry certificates)."""
import gzip
import json
import random
import unittest

from exact_theta import InvalidCertificate
import rank4_pairing_search as S
import rank4_six_minima as M

G9373 = [[6, 0, 1, 2], [0, 6, -2, 1], [1, -2, 16, 0], [2, 1, 0, 18]]
G121 = [[4, -2, -1, -1], [-2, 4, 0, 1], [-1, 0, 4, 2], [-1, 1, 2, 4]]


class CommonLattice(unittest.TestCase):
    def test_structure(self):
        for G in (S.DET21, G9373, G121):
            for p in S.isotropic_planes(G):
                plane, B, GK, H, W, pairs, tau = M.common_lattice(G, p)
                self.assertEqual(abs(S.determinant(B)), 4)
                self.assertEqual(S.determinant(H), S.determinant(G))
                self.assertEqual(len(W), 4)
                self.assertEqual(sorted(k for pr in pairs.values() for k, _ in pr), sorted(set(range(16)) - set(W)))

    def test_rejects_bad_planes(self):
        with self.assertRaises(InvalidCertificate):
            M.common_lattice(S.DET21, (1, 2, 4))          # not a subspace
        bad = next((u, v) for u in range(1, 16) for v in range(u + 1, 16)
                   if sum(S.BITS[u][i] * S.DET21[i][j] * S.BITS[v][j] for i in range(4) for j in range(4)) % 2)
        with self.assertRaises(InvalidCertificate):
            M.common_lattice(S.DET21, (bad[0], bad[1], bad[0] ^ bad[1]))


class SixMinima(unittest.TestCase):
    def test_agrees_with_census_subcosets(self):
        res = json.load(gzip.open("results/rank4-joint-depth-production/per-class-1021-1997.json.gz", "rt"))
        random.seed(11)
        for r in random.sample(res, 12):
            for p in r["planes"]:
                six = {(x["shift"], x["residue_mod_8"]): (x["min_norm4"], x["plus"] + x["minus"], x["signed"])
                       for x in M.six_minima(r["G"], p["plane"])}
                for sh in p["shifts"]:
                    k = S.BITS.index(tuple(sh["a"]))
                    for cls, v in sh["subcosets"].items():
                        self.assertEqual(six[(k, int(cls))], tuple(v["first"]))

    def test_det9373_is_five_of_six_lone_pair(self):
        six = M.six_minima(G9373, (3, 12, 15))
        self.assertEqual(sum(x["balanced"] for x in six), 5)
        bad = [x for x in six if not x["balanced"]][0]
        self.assertEqual((bad["shift"], bad["min_norm4"], bad["plus"] + bad["minus"], abs(bad["signed"])), (3, 64, 2, 2))

    def test_minimizers_are_complete_and_minimal(self):
        for x in M.six_minima(S.DET21, (5, 10, 15)):
            for v in x["plus_minimizers_L"] + x["minus_minimizers_L"]:
                self.assertEqual(S.norm(S.DET21, v), x["min_norm4"])
                self.assertEqual(M.label(v), x["shift"])
            self.assertEqual(len(x["plus_minimizers_L"]), x["plus"])
            self.assertEqual(len(x["minus_minimizers_L"]), x["minus"])


class Certificates(unittest.TestCase):
    def test_det21_and_det121_have_checked_certificates(self):
        for G, plane in ((S.DET21, (5, 10, 15)), (G121, (6, 11, 13))):
            self.assertEqual(sum(x["balanced"] for x in M.six_minima(G, plane)), 6)
            cert = M.find_marked_certificate(G, plane)
            self.assertIsNotNone(cert)
            self.assertTrue(M.check_marked_certificate(cert))

    def test_near_miss_has_no_certificate(self):
        self.assertIsNone(M.find_marked_certificate(G9373, (3, 12, 15)))

    def test_forged_certificates_rejected(self):
        cert = M.find_marked_certificate(S.DET21, (5, 10, 15))
        I = [[int(i == j) for j in range(4)] for i in range(4)]
        for forged in (dict(cert, T=I),                                            # identity is not tau mod 2
                       dict(cert, T=[[2 * x for x in r] for r in cert["T"]]),       # not unimodular
                       dict(cert, B=[[2 * x for x in r] for r in cert["B"]]),       # wrong index
                       dict(cert, plane=[3, 12, 15])):                              # another plane
            with self.assertRaises(InvalidCertificate):
                M.check_marked_certificate(forged)

    def test_bad_gram_inputs_rejected(self):
        """RANK4-SIX-CELLS-QUICK-REVIEW-ASTRA.md section 3: -G, non-symmetric, non-integer and odd forms."""
        cert = M.find_marked_certificate(S.DET21, (5, 10, 15))
        neg = dict(cert, G=[[-x for x in r] for r in cert["G"]])
        asym = dict(cert, G=[r[:] for r in cert["G"]])
        asym["G"][0][1] += 2
        flt = dict(cert, G=[[float(x) for x in r] for r in cert["G"]])
        boo = dict(cert, T=[[bool(x) if x in (0, 1) else x for x in r] for r in cert["T"]])
        odd = dict(cert, G=[r[:] for r in cert["G"]])
        odd["G"][0][0] += 1
        for bad in (neg, asym, flt, boo, odd):
            with self.assertRaises(InvalidCertificate):
                M.check_marked_certificate(bad)
        from verify_rank4_feasibility import box_isometries
        import verify_rank4_parity as V
        h = V.splitting_involution(box_isometries(S.DET21, S.DET21, 5_000_000))
        with self.assertRaises(InvalidCertificate):
            M.check_splitting([[-x for x in r] for r in S.DET21], h)

    def test_splitting_witness(self):
        from verify_rank4_feasibility import box_isometries
        import verify_rank4_parity as V
        h = V.splitting_involution(box_isometries(S.DET21, S.DET21, 5_000_000))
        self.assertTrue(M.check_splitting(S.DET21, h))
        I = [[int(i == j) for j in range(4)] for i in range(4)]
        with self.assertRaises(InvalidCertificate):
            M.check_splitting(S.DET21, I)
        cert = M.find_marked_certificate(G121, (6, 11, 13))
        with self.assertRaises(InvalidCertificate):                    # h^2 = -I is not a splitting
            M.check_splitting(G121, cert["h"])


if __name__ == "__main__":
    unittest.main()
