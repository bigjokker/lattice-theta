"""Rank-four generation, independent controls, boundary and tampering regressions."""

from copy import deepcopy
from itertools import permutations
import json
from math import gcd
from pathlib import Path
import unittest

from exact_theta import (EnumerationLimit, InvalidCertificate, identity, inverse,
                         matvec, multiply, transpose)
from rank4_feasibility import candidates, cutoff_build, determinant, isometries, modular_data
from verify_rank4_feasibility import (box_isometries, controls_replay, cutoff_replay,
                                      det4, independent_candidates)

ROOT=Path(__file__).resolve().parent.parent


class RankFourFeasibilityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.controls=json.loads((ROOT/"data/rank4-feasibility/controls.json").read_text(encoding="utf-8"))
        cls.cutoffs=json.loads((ROOT/"data/rank4-feasibility/cutoff-controls.json").read_text(encoding="utf-8"))

    def test_independent_candidate_domains_and_d4_coverage(self):
        for H in range(1,25):
            with self.subTest(H=H):
                self.assertEqual(candidates(H),independent_candidates(H))
        self.assertEqual(len(candidates(24)),1510)
        self.assertEqual(candidates(1),[identity(4)])
        D4=self.controls["lattices"][3]["G"]
        self.assertTrue(any(isometries(G,D4) for G in candidates(4) if determinant(G)==4))

    def test_general_determinants_and_noncanonical_full_groups(self):
        for perm in permutations(range(4)):
            T=transpose([identity(4)[i] for i in perm])
            self.assertEqual(determinant(T),det4(T))
            self.assertEqual(abs(determinant(T)),1)
        self.assertEqual(determinant([[1]*4 for _ in range(4)]),0)
        U=[[1,1,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
        H=multiply(transpose(U),U)
        group=isometries(H,H)
        self.assertEqual(group,box_isometries(H,H))
        self.assertEqual(len(group),384)

    def test_all_named_controls_full_stabilizers_and_counts(self):
        report=controls_replay(self.controls)
        self.assertEqual(report["verdict"],"verified")
        self.assertTrue(report["all_groups_and_stabilizers_replayed"])
        self.assertFalse(report["full_census_performed"])
        expected=[(384,55,81),(384,55,81),(240,0,136),(1152,9,127),(96,28,108)]
        for row,(order,zero,nonzero) in zip(report["controls"],expected):
            self.assertEqual(row["group_order"],order)
            self.assertEqual(row["even_verdicts"].get("proved_zero",0),zero)
            self.assertEqual(row["even_verdicts"].get("proved_nonzero",0),nonzero)
            self.assertEqual(row["all_pairs"],256)

    def test_cutoff_controls_endpoint_and_cancelled_first_shell(self):
        results=[cutoff_replay(cert) for cert in self.cutoffs["certificates"]]
        self.assertTrue(all(r["certificate_replay_complete"] for r in results))
        self.assertEqual(sum(r["verdict"]=="unresolved" for r in results),3)
        self.assertTrue(all(r["full_stabilizer_checked"] is False for r in results))
        self.assertEqual(results[0]["verdict"],"proved_zero_by_modular_cutoff")
        trap=self.cutoffs["certificates"][-2]
        later=self.cutoffs["certificates"][-1]
        self.assertEqual(trap["modular"]["cutoff"],384)
        self.assertEqual(trap["evidence"]["verdict"],"unresolved")
        self.assertEqual((later["evidence"]["nonzero_norm4"],later["evidence"]["coefficient"]),(10,2))

    def test_kernel_lattice_basis_changes_and_exact_indices(self):
        a,b=[1,1,0,0],[1,1,0,0]
        cert=cutoff_build(identity(4),a,b)
        V=[[1,1,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
        m=cert["modular"]
        m["index_two_basis"]=multiply(m["index_two_basis"],V)
        m["G0"]=multiply(transpose(m["index_two_basis"]),m["index_two_basis"])
        self.assertEqual(cutoff_replay(cert)["verdict"],"proved_zero_by_modular_cutoff")
        H=multiply(transpose(V),V)
        ap=[int(v)%2 for v in matvec(inverse(V),a)]
        bp=[v%2 for v in matvec(transpose(V),b)]
        self.assertEqual(cutoff_replay(cutoff_build(H,ap,bp))["verdict"],"proved_zero_by_modular_cutoff")
        for M in (16,32,48,64):
            N=M//16
            G=[[N*int(i==j) for j in range(4)] for i in range(4)]
            item=modular_data(G,[0]*4,[1,1,0,0])
            # This kernel has N0=2N, so its actual level is 2M.
            actual=2*M
            prim=sum(gcd(gcd(x,y),actual)==1 for x in range(actual) for y in range(actual))
            unit=sum(gcd(x,actual)==1 for x in range(actual))
            self.assertEqual(item["index"],prim//unit)
            self.assertEqual(item["cutoff"]*6,item["index"])

    def test_resource_limits_never_certify_partial_results(self):
        with self.assertRaises(EnumerationLimit):
            candidates(24,max_nodes=1)
        for method in (isometries,box_isometries):
            with self.assertRaises(EnumerationLimit):
                method(identity(4),identity(4),max_nodes=1)
        self.assertEqual(controls_replay(self.controls,max_nodes=1)["verdict"],"unresolved")
        cert=cutoff_build(identity(4),[1,1,0,0],[1,1,0,0],max_nodes=1)
        self.assertEqual(cutoff_replay(cert)["verdict"],"unresolved")
        self.assertNotIn("enumeration",cert["evidence"])
        complete=cutoff_build(identity(4),[1,1,0,0],[1,1,0,0])
        self.assertFalse(cutoff_replay(complete,max_nodes=1)["certificate_replay_complete"])

    def test_reject_tampered_metadata_coefficients_and_hypotheses(self):
        cert=cutoff_build(identity(4),[1,1,0,0],[1,1,0,0])
        for key,value in (("N0",1),("index",1),("cutoff",7),("weight","3/2"),
                          ("character_discriminant",-16),("unsigned_level",1),("theorem","unknown")):
            changed=deepcopy(cert)
            changed["modular"][key]=value
            with self.assertRaises(InvalidCertificate):
                cutoff_replay(changed)
        changed=deepcopy(cert)
        changed["evidence"]["enumeration"]["shells"].pop()
        with self.assertRaises(InvalidCertificate):
            cutoff_replay(changed)
        changed=deepcopy(cert)
        changed["run"]["norm4_bound"]=7
        with self.assertRaises(InvalidCertificate):
            cutoff_replay(changed)
        for a,b in (([0]*4,[0]*4),([1,0,0,0],[1,0,0,0]),([False,0,0,0],[1,0,0,0])):
            with self.assertRaises(InvalidCertificate):
                cutoff_build(identity(4),a,b)
        for H in (0,25,True):
            with self.assertRaises(InvalidCertificate):
                candidates(H)

    def test_reject_missing_group_elements_pairs_and_signs(self):
        for mutate in (lambda d:d["lattices"][0]["automorphisms"].pop(),
                       lambda d:d["lattices"][0]["pairs"].pop(),
                       lambda d:d["lattices"][0]["pairs"][0]["stabilizer"]["epsilon"].__setitem__(0,1),
                       lambda d:d["lattices"][0]["pairs"][0]["a"].__setitem__(0,False)):
            changed=deepcopy(self.controls)
            mutate(changed)
            with self.assertRaises(InvalidCertificate):
                controls_replay(changed)


if __name__=="__main__":
    unittest.main()
