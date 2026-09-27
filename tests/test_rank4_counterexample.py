from copy import deepcopy
import unittest
from exact_theta import EnumerationLimit,InvalidCertificate,identity
from rank4_counterexample import build,replay


class CounterexampleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.doc=build()

    def test_all_level_identity_and_full_group(self):
        result=replay(self.doc)
        self.assertEqual(result["determinant"],23716)
        self.assertEqual(result["full_group_order"],2)
        self.assertTrue(result["all_stabilizer_epsilon_zero"])
        self.assertFalse(result["bounded_cancellation_used_as_identity_proof"])

    def test_corrobating_bound_not_used_to_prove_identity(self):
        from exact_theta import complete_shells
        doc=deepcopy(self.doc)
        doc["corroboration"]=complete_shells(doc["G"],doc["a"],doc["b"],0)
        self.assertEqual(replay(doc)["verdict"],"verified_rank4_symmetry_converse_counterexample")

    def test_false_ambient_pairing_and_metric(self):
        for change in [lambda d:d.__setitem__("ambient_isometry",identity(4)),
                       lambda d:d["G"][0].__setitem__(0,13),
                       lambda d:d["index_two_basis"][0].__setitem__(0,4)]:
            doc=deepcopy(self.doc); change(doc)
            with self.assertRaises(InvalidCertificate):
                replay(doc)

    def test_group_coefficient_and_offset_tampering(self):
        for change in [lambda d:d["automorphisms"].pop(),
                       lambda d:d["corroboration"]["shells"][0].__setitem__("signed_coefficient",2),
                       lambda d:d["negative_coset_representative"].__setitem__(0,0)]:
            doc=deepcopy(self.doc); change(doc)
            with self.assertRaises(InvalidCertificate):
                replay(doc)

    def test_interrupted_replay_cannot_certify(self):
        with self.assertRaises(EnumerationLimit):
            replay(self.doc,max_nodes=1)


if __name__=="__main__":
    unittest.main()
