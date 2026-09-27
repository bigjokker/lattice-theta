"""Exact inclusive cutoff, independent shell replay and adverse evidence checks."""

from copy import deepcopy
from fractions import Fraction
from math import gcd
from pathlib import Path
import json
import unittest

from exact_theta import InvalidCertificate, identity, inverse, matvec, multiply, transpose
from rank3_cutoff import build, cutoff, index_gamma0, replay


class CutoffTests(unittest.TestCase):
    def test_exact_index_and_inclusive_cutoffs(self):
        for N, expected in ((1, 3), (2, 6), (4, 12), (8, 24), (16, 48), (80, 288), (96, 384)):
            self.assertEqual(cutoff(N), expected)
            self.assertEqual(index_gamma0(16*N), 8*expected)
        # Independently count primitive pairs mod M, modulo unit scaling.
        for M in (16, 32, 48, 64, 80):
            primitive = sum(gcd(gcd(a, b), M) == 1 for a in range(M) for b in range(M))
            units = sum(gcd(a, M) == 1 for a in range(M))
            self.assertEqual(index_gamma0(M), primitive // units)
        for N in (0, -1, True, 1.0):
            with self.assertRaises(InvalidCertificate):
                cutoff(N)

    def test_modular_zero_inclusive_endpoint_and_interruption(self):
        G, a, b = identity(3), [1, 1, 0], [1, 1, 0]
        complete = build(G, a, b)
        B = complete["modular"]["cutoff"]
        self.assertEqual(B, 6)
        self.assertEqual(replay(complete)["verdict"], "proved_zero_by_modular_cutoff")
        self.assertEqual(replay(complete, backend="ldl")["shells"], replay(complete)["shells"])
        self.assertFalse(replay(complete)["full_stabilizer_checked"])
        shallow = build(G, a, b, norm4_bound=B-1)
        self.assertEqual(replay(shallow)["verdict"], "unresolved")
        interrupted = build(G, a, b, max_nodes=1)
        self.assertFalse(interrupted["evidence"]["enumeration_complete"])
        self.assertNotIn("enumeration", interrupted["evidence"])
        self.assertEqual(replay(interrupted)["verdict"], "unresolved")
        self.assertFalse(replay(complete, max_nodes=1)["certificate_replay_complete"])

    def test_levels_constant_term_nonzero_and_first_shell_trap(self):
        for G, expectedN, expectedB in ((identity(3), 4, 12),
                ([[2,-1,0],[-1,2,-1],[0,-1,2]], 16, 48),
                ([[2,0,0],[0,2,0],[0,0,2]], 8, 24)):
            cert = build(G, [0, 1, 0], [1, 0, 0])
            self.assertEqual((cert["modular"]["N0"], cert["modular"]["cutoff"]), (expectedN, expectedB))
            box, ldl = replay(cert), replay(cert, backend="ldl")
            self.assertEqual(box["verdict"], "proved_nonzero")
            self.assertEqual(box["shells"], ldl["shells"])
        cert = build(identity(3), [0, 0, 0], [1, 0, 0], norm4_bound=0)
        self.assertEqual(replay(cert)["nonzero_norm4"], 0)
        G = [[2,-1,-1],[-1,4,-1],[-1,-1,4]]
        shallow = build(G, [0, 1, 1], [1, 0, 0], norm4_bound=6)
        self.assertEqual((shallow["modular"]["N0"], shallow["modular"]["cutoff"]), (80, 288))
        self.assertEqual(replay(shallow)["verdict"], "unresolved")
        full = replay(build(G, [0, 1, 1], [1, 0, 0]))
        self.assertEqual((full["nonzero_norm4"], full["coefficient"]), (10, 2))

    def test_kernel_basis_and_lattice_basis_changes(self):
        cert = build(identity(3), [1, 1, 0], [1, 1, 0])
        V = [[1, 2, 0], [0, 1, 0], [0, 0, 1]]
        m = cert["modular"]
        m["index_two_basis"] = multiply(m["index_two_basis"], V)
        m["G0"] = multiply(multiply(transpose(m["index_two_basis"]), identity(3)), m["index_two_basis"])
        m["cleared_inverse"] = [[int(m["N0"]*v) for v in row] for row in inverse(m["G0"])]
        self.assertEqual(replay(cert)["verdict"], "proved_zero_by_modular_cutoff")
        U = [[1, 1, 0], [0, 1, 0], [0, 0, 1]]
        H = multiply(transpose(U), U)
        a = [int(v) % 2 for v in matvec(inverse(U), [1, 1, 0])]
        b = [v % 2 for v in matvec(transpose(U), [1, 1, 0])]
        self.assertEqual(replay(build(H, a, b))["verdict"], "proved_zero_by_modular_cutoff")

    def test_reject_wrong_hypotheses_and_tampered_certificates(self):
        for a, b in (([1,0,0], [1,0,0]), ([0,0,0], [0,0,0]), ([2,0,0], [1,0,0]),
                     ([False,0,0], [1,0,0])):
            with self.assertRaises(InvalidCertificate):
                build(identity(3), a, b)
        with self.assertRaises(InvalidCertificate):
            build([[1,0,0],[0,1,0],[0,0,25]], [0,0,0], [1,0,0])
        cert = build(identity(3), [1,1,0], [1,1,0])
        for key, value in (("N0", 1), ("index", 1), ("cutoff", 5), ("unsigned_level", 16),
                           ("cleared_inverse", identity(3)), ("theorem", "unknown"),
                           ("character_discriminant", -16), ("weight", "3")):
            doc = deepcopy(cert)
            doc["modular"][key] = value
            with self.assertRaises(InvalidCertificate):
                replay(doc)
        doc = deepcopy(cert)
        doc["evidence"]["enumeration"]["shells"].pop()
        with self.assertRaises(InvalidCertificate):
            replay(doc)
        doc = build(identity(3), [0,0,0], [1,0,0], norm4_bound=0)
        doc["evidence"]["enumeration"]["shells"].clear()
        doc["evidence"]["verdict"] = "unresolved"
        with self.assertRaises(InvalidCertificate):
            replay(doc)
        doc = build(identity(3), [1,1,0], [1,1,0], max_nodes=1)
        doc["evidence"]["verdict"] = "proved_zero_by_modular_cutoff"
        with self.assertRaises(InvalidCertificate):
            replay(doc)


if __name__ == "__main__":
    unittest.main()
