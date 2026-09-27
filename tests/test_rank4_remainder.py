"""Tests of projection certificates, cancellation evidence and resource limits."""
from copy import deepcopy
import unittest

from exact_theta import EnumerationLimit, InvalidCertificate
from rank4_cancellation import build, gram, replay
from rank4_remainder import direct_average_check, split_search, validate_projector


class RemainderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.controls = build()

    def test_later_shell_controls_inside_and_outside_census(self):
        self.assertEqual(replay(self.controls)["cases"],21)

    def test_initial_zero_and_endpoint_tampering(self):
        doc = deepcopy(self.controls)
        doc["cases"][0]["claim"]["coefficient"] = 0
        with self.assertRaises(InvalidCertificate):
            replay(doc)
        doc = deepcopy(self.controls)
        doc["cases"][0]["enumeration"]["norm4_bound"] -= 1
        with self.assertRaises(InvalidCertificate):
            replay(doc)
        doc = deepcopy(self.controls)
        doc["cases"][0]["enumeration"]["shells"][0]["signed_coefficient"] = 2
        with self.assertRaises(InvalidCertificate):
            replay(doc)

    def test_hypotheses_and_case_coverage(self):
        with self.assertRaises(InvalidCertificate):
            gram(2,3)
        with self.assertRaises(InvalidCertificate):
            gram(True,3)
        doc = deepcopy(self.controls)
        doc["cases"].pop()
        with self.assertRaises(InvalidCertificate):
            replay(doc)

    def test_line_plane_and_no_summand_controls(self):
        I = [[int(i==j) for j in range(4)] for i in range(4)]
        A2A2 = [[2,-1,0,0],[-1,2,0,0],[0,0,2,-1],[0,0,-1,2]]
        D4 = [[2,-1,-1,-1],[-1,2,0,0],[-1,0,2,0],[-1,0,0,2]]
        for G,bound,kind in [(I,1,"line"),(A2A2,9,"plane"),(D4,4,"indecomposable")]:
            left = split_search(G,bound,"ldl",1_000_000)
            right = split_search(G,bound,"box",1_000_000)
            self.assertEqual(left,right)
            self.assertEqual(left["kind"],kind)
            if kind != "indecomposable":
                validate_projector(G,left)

    def test_false_projector_and_interruption(self):
        I = [[int(i==j) for j in range(4)] for i in range(4)]
        item = split_search(I,1,"ldl",1_000_000)
        item["projector"][0][0] += 1
        with self.assertRaises(InvalidCertificate):
            validate_projector(I,item)
        with self.assertRaises(EnumerationLimit):
            split_search(I,1,"ldl",1)
        with self.assertRaises(EnumerationLimit):
            replay(self.controls,max_nodes=1)

    def test_direct_average_rejects_wrong_residue_weights(self):
        I = [[int(i==j) for j in range(4)] for i in range(4)]
        minus = [[-int(i==j) for j in range(4)] for i in range(4)]
        # Odd pairing gives a zero average under {I,-I}.
        direct_average_check([I,minus],[1,0,0,0],[1,0,0,0],[])
        with self.assertRaises(InvalidCertificate):
            direct_average_check([I,minus],[1,0,0,0],[1,0,0,0],
                                 [{"residue":[1,0,0,0],"numerator":1}])


if __name__ == "__main__":
    unittest.main()
