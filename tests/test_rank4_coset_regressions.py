from copy import deepcopy
import unittest
from exact_theta import EnumerationLimit,InvalidCertificate
from rank4_coset_regressions import build,replay


class CosetRegressionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.document=build()

    def test_complete_counts_and_swap(self):
        result=replay(self.document)
        self.assertEqual(result["verdict"],"verified")
        self.assertEqual(result["cases"][0]["first_difference_norm4"],14)
        self.assertEqual(result["cases"][1]["first_difference_norm4"],8)
        self.assertEqual(result["cases"][2]["verdict"],"proved_equal_by_coset_swap")

    def test_count_and_partial_group_tampering(self):
        for change in [lambda d:d["cases"][0]["coset_counts"][0].__setitem__("plus_count",0),
                       lambda d:d["cases"][0]["automorphisms"].pop()]:
            doc=deepcopy(self.document)
            change(doc)
            with self.assertRaises(InvalidCertificate):
                replay(doc)

    def test_missing_case_and_false_swap(self):
        doc=deepcopy(self.document)
        doc["cases"].pop()
        with self.assertRaises(InvalidCertificate):
            replay(doc)
        doc=deepcopy(self.document)
        doc["cases"][2]["pair"]["evidence"]["T"]=[[-int(i==j) for j in range(4)] for i in range(4)]
        with self.assertRaises(InvalidCertificate):
            replay(doc)

    def test_complement_and_interruption(self):
        doc=deepcopy(self.document)
        doc["cases"][0]["u"]=[0]*4
        with self.assertRaises(InvalidCertificate):
            replay(doc)
        with self.assertRaises(EnumerationLimit):
            build(max_nodes=1)
        with self.assertRaises(EnumerationLimit):
            replay(self.document,max_nodes=1)


if __name__=="__main__":
    unittest.main()
