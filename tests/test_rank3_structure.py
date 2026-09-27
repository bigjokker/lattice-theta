"""Independent structural replay and uniform graph-shell boundary regressions."""

from copy import deepcopy
import hashlib
from itertools import combinations, product
import json
from pathlib import Path
import unittest

from exact_theta import (EnumerationLimit, InvalidCertificate, dot, identity, inverse,
                         matvec, multiply, transpose)
from build_rank3_structure import build
from rank3_structure import (complete_primitive, decomposition, graph_type,
                             obtuse_superbase, shell_claim)
from verify_lattice_census import box_shells, box_vectors
from verify_rank3_structure import check_superbase, det, replay, splitting_vectors

ROOT = Path(__file__).resolve().parent.parent


def grounded_gram(weights):
    W = [[0]*4 for _ in range(4)]
    for (i, j), value in weights.items():
        W[i][j] = W[j][i] = value
    return [[sum(W[i]) if i == j else -W[i][j] for j in range(1,4)] for i in range(1,4)]


class StructuralTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        raw = (ROOT / "data/census-rank3/det24-bound32.json").read_bytes()
        cls.source, cls.digest = json.loads(raw), hashlib.sha256(raw).hexdigest()
        cls.pack = json.loads((ROOT / "data/census-rank3/structure-det24.json").read_text(encoding="utf-8"))

    def test_full_structural_replay_and_counts(self):
        result = replay(self.pack, self.source, self.digest)
        self.assertEqual(result["verdict"], "resolved")
        self.assertTrue(result["certificate_replay_complete"])
        s = result["summary"]
        self.assertEqual(s["structure_types"], {"three_lines": 51, "line_and_indecomposable_binary": 44, "indecomposable": 25})
        self.assertEqual((s["all_pairs"], s["even_pairs"]), (7680, 4320))
        self.assertEqual(s["positive_graphs"], {"cycle": 12, "diamond": 11, "complete": 2})
        self.assertEqual(s["symbolic_shell_cases"]["diamond_later"], 22)

    def test_saturated_completion_and_noncanonical_decompositions(self):
        for v in ([1,0,0], [-1,0,0], [0,-1,0], [0,0,1], [2,3,5], [-7,11,13]):
            U = complete_primitive(v)
            self.assertEqual(transpose(U)[0], v)
            self.assertEqual(abs(det(U)), 1)
            self.assertTrue(all(x.denominator == 1 for row in inverse(U) for x in row))
        for v in ([0,0,0], [2,4,6]):
            with self.assertRaises(InvalidCertificate):
                complete_primitive(v)
        B = [[1,2,-1],[0,1,1],[0,0,1]]
        for G, kind in ((identity(3), "three_lines"),
                        ([[1,0,0],[0,2,1],[0,1,2]], "line_and_indecomposable_binary")):
            H = multiply(multiply(transpose(B), G), B)
            result = decomposition(H)
            self.assertEqual(result["split"]["type"], kind)
            self.assertTrue(splitting_vectors(H, 1_000_000))
            U = result["split"]["basis"]
            self.assertEqual(abs(det(U)), 1)
            self.assertEqual(multiply(multiply(transpose(U), H), U), result["split"]["block_gram"])

    def test_all_graph_types_equal_conorms_and_scaling(self):
        graphs = [({(0,2):1,(0,3):1,(1,2):1,(1,3):1}, "cycle"),
                  ({(0,1):1,(0,2):1,(0,3):1,(1,2):1,(1,3):1}, "diamond"),
                  ({edge:1 for edge in combinations(range(4),2)}, "complete"),
                  ({(0,1):2,(0,2):1,(0,3):7,(1,2):3,(1,3):5}, "diamond")]
        # Isolated symbolic-proof controls, not an expanded determinant census.
        for weights, kind in graphs:
            for scale in (1, 1000):
                G = [[scale*x for x in row] for row in grounded_gram(weights)]
                self.assertEqual(decomposition(G)["verdict"], "indecomposable")
                S = obtuse_superbase(G)
                self.assertEqual(graph_type(S["conorms"]), kind)
                self.assertEqual(check_superbase(G, S, 1_000_000), kind)
                for a in product((0,1), repeat=3):
                    for b in product((0,1), repeat=3):
                        if dot(a,b) % 2:
                            continue
                        claim = shell_claim(G, list(a), list(b), S)
                        shell = next(s for s in box_shells(G, a, b, claim["norm4"])
                                     if s["norm4"] == claim["norm4"])
                        self.assertEqual(shell["signed_coefficient"], claim["coefficient"])

    def test_determinant29_first_shell_boundary_and_basis_changes(self):
        G = json.loads((ROOT / "data/ternary.json").read_text(encoding="utf-8"))["G"]
        self.assertEqual(det(G), 29)
        self.assertEqual(decomposition(G)["verdict"], "indecomposable")
        a = b = [1,0,1]
        claim = shell_claim(G, a, b, obtuse_superbase(G))
        self.assertEqual((claim["case"], claim["norm4"], claim["coefficient"]), ("diamond_later", 11, -2))
        U = [[1,2,-1],[0,1,1],[0,0,1]]
        H = multiply(multiply(transpose(U), G), U)
        ap = [int(x) % 2 for x in matvec(inverse(U), a)]
        bp = [x % 2 for x in matvec(transpose(U), b)]
        S = obtuse_superbase(H)
        self.assertEqual(check_superbase(H, S, 1_000_000), "diamond")
        claim = shell_claim(H, ap, bp, S)
        actual = next(s for s in box_shells(H, ap, bp, claim["norm4"]) if s["norm4"] == claim["norm4"])
        self.assertEqual(actual["signed_coefficient"], claim["coefficient"])

    def test_later_shell_contains_additional_cancelling_vectors(self):
        a, b = [0,1,1], [0,1,1]
        for h, N, count in ((4,20,10), (8,36,6)):
            G = grounded_gram({(0,1):h,(0,2):1,(0,3):1,(1,2):1,(1,3):1})
            claim = shell_claim(G, a, b, obtuse_superbase(G))
            self.assertEqual((claim["case"], claim["norm4"], claim["coefficient"]), ("diamond_later",N,2))
            self.assertEqual(len(claim["selected_vectors"]), 2)
            full = [y for y in box_vectors(G, N, 1_000_000)
                    if dot(y, matvec(G,y)) == N and all((yi-ai)%2 == 0 for yi,ai in zip(y,a))]
            self.assertEqual(len(full), count)
            selected = {tuple(y) for y in claim["selected_vectors"]}
            extras = [y for y in full if y not in selected]
            self.assertEqual(len(extras), count-2)
            self.assertTrue(all(y[0] == 0 for y in extras))  # Root 0 fixed; x_t is y[0].
            self.assertEqual(sum((-1)**(dot([(yi-ai)//2 for yi,ai in zip(y,a)],b)%2) for y in extras), 0)
            actual = next(s for s in box_shells(G,a,b,N) if s["norm4"] == N)
            self.assertEqual((actual["vector_count"], actual["signed_coefficient"]), (count,2))
        G = [[2,0,-1],[0,2,-1],[-1,-1,6]]
        self.assertEqual(self.source["lattices"][97]["G"], G)
        a = [1,1,0]
        for b, expected in (([1,1,0],2), ([1,1,1],-2)):
            claim = shell_claim(G,a,b,obtuse_superbase(G))
            self.assertEqual((claim["norm4"],claim["coefficient"]), (20,expected))
            shells = box_shells(G,a,b,20)
            self.assertEqual(next(s["signed_coefficient"] for s in shells if s["norm4"] == 4), 0)
            later = next(s for s in shells if s["norm4"] == 20)
            self.assertEqual((later["vector_count"],later["signed_coefficient"]), (10,expected))

    def test_limits_preserve_uncertainty(self):
        with self.assertRaises(EnumerationLimit):
            decomposition(identity(3), max_nodes=1)
        with self.assertRaises(EnumerationLimit):
            obtuse_superbase([[3,1,0],[1,3,-1],[0,-1,4]], max_steps=0)
        limited = build(self.source, self.digest, max_nodes=1)
        self.assertFalse(limited["summary"]["complete"])
        self.assertTrue(all(row["pairs"] == [] for row in limited["rows"]))
        self.assertEqual(replay(limited, self.source, self.digest)["verdict"], "unresolved")
        result = replay(self.pack, self.source, self.digest, max_nodes=1)
        self.assertFalse(result["certificate_replay_complete"])
        self.assertEqual(result["verdict"], "unresolved")

    def test_reject_incomplete_or_tampered_structure(self):
        def reject(mutate):
            doc = deepcopy(self.pack)
            mutate(doc)
            with self.assertRaises(InvalidCertificate):
                replay(doc, self.source, self.digest)
        reject(lambda d: d["rows"].pop())
        reject(lambda d: d["rows"][0]["pairs"].pop())
        reject(lambda d: d["rows"][0]["decomposition"]["split"]["basis"][0].__setitem__(1, 1))
        reject(lambda d: d["rows"][0]["decomposition"].__setitem__("verdict", "indecomposable"))
        positive = next(i for i,r in enumerate(self.pack["rows"]) if r["superbase"])
        reject(lambda d: d["rows"][positive]["superbase"]["conorms"][0].__setitem__(1, 999))
        reject(lambda d: d["rows"][positive]["pairs"][0]["shell"].__setitem__("coefficient", 2))
        reject(lambda d: d["rows"][positive]["superbase"]["steps"].append([0,1]))
        reject(lambda d: d.__setitem__("census_sha256", "wrong"))


if __name__ == "__main__":
    unittest.main()
