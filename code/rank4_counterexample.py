"""Exact all-level coset identity and independent full-group counterexample replay."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from exact_theta import (complete_shells,dot,identity,integer,inverse,ldl,matrix,
                         matvec,multiply,require,symmetry_data,transpose,vector)
from rank4_feasibility import admit,determinant,isometries
from verify_lattice_census import box_shells,box_vectors
from verify_rank4_feasibility import box_isometries,det4

M = [[8,-4,1,2],[-4,8,-3,1],[1,-3,14,-7],[2,1,-7,14]]
V = [[-1,0,0,0],[1,2,0,0],[0,0,2,0],[0,0,0,2]]
R = [[0,-1,0,0],[1,-1,0,0],[0,0,0,-1],[0,0,1,-1]]
G = [[12,12,-4,-1],[12,16,-6,2],[-4,-6,28,-14],[-1,2,-14,28]]
A,B = [0,1,0,0],[1,0,0,0]


def build():
    K = [[2,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
    W = [[-1,0,0,0],[1,1,0,0],[0,0,1,0],[0,0,0,1]]
    return {"schema":"work7-rank4-ambient-coset-counterexample-v1", "G":G,"a":A,"b":B,
            "ambient_metric":M,"embedding":V,"ambient_isometry":R,
            "index_two_basis":K,"common_coset_basis_change":W,"complementary_vector":[1,0,0,0],
            "positive_coset_representative":[0,2,0,0],"negative_coset_representative":[2,0,0,0],
            "automorphisms":isometries(G,G),
            "corroboration":complete_shells(G,A,B,64)}


def replay(doc,max_nodes=1_000_000):
    require(doc["schema"] == "work7-rank4-ambient-coset-counterexample-v1","wrong schema")
    G = matrix(doc["G"],"G",4)
    M = matrix(doc["ambient_metric"],"ambient metric",4)
    _,pivots = ldl(G)
    ldl(M)
    a,b = vector(doc["a"],"a",4),vector(doc["b"],"b",4)
    require(all(x in (0,1) for x in a+b) and any(a) and any(b) and dot(a,b)%2 == 0,
            "requires a nontrivial even characteristic")
    V = matrix(doc["embedding"],"embedding",4)
    R = matrix(doc["ambient_isometry"],"ambient isometry",4)
    K = matrix(doc["index_two_basis"],"kernel basis",4)
    W = matrix(doc["common_coset_basis_change"],"coset basis change",4)
    u = vector(doc["complementary_vector"],"complement",4)
    plus = vector(doc["positive_coset_representative"],"positive offset",4)
    minus = vector(doc["negative_coset_representative"],"negative offset",4)
    require(det4(V) != 0,"singular embedding")
    require(multiply(multiply(transpose(V),M),V) == [[2*x for x in row] for row in G],
            "wrong norm relation")
    require(abs(det4(K)) == 2 and all(x%2 == 0 for x in matvec(transpose(K),b)),
            "kernel basis is not the full index-two kernel")
    require(dot(u,b)%2 == 1,"wrong complementary vector")
    require(abs(det4(W)) == 1 and multiply(V,K) == [[2*x for x in row] for row in W],
            "common coset difference lattice is not four Z4")
    actual_plus = matvec(V,a)
    actual_minus = matvec(V,[x+2*y for x,y in zip(a,u)])
    require(all((x-y)%4 == 0 for x,y in zip(actual_plus,plus))
            and all((x-y)%4 == 0 for x,y in zip(actual_minus,minus)),"wrong coset offsets")
    require(any((x-y)%4 for x,y in zip(plus,minus)),"cosets are not distinct")
    require(abs(det4(R)) == 1 and multiply(multiply(transpose(R),M),R) == M,
            "ambient map is not a norm-preserving unimodular isometry")
    require(all((x-y)%4 == 0 for x,y in zip(matvec(R,minus),plus)),
            "ambient map does not take negative coset onto positive coset")
    # The preceding identities prove equality at ALL norms. The following box
    # check proves absence of an epsilon-one automorphism of the original lattice.
    group = box_isometries(G,G,max_nodes)
    for T in doc["automorphisms"]:
        matrix(T,"stored automorphism",4)
    require(group == doc["automorphisms"],"incorrect full automorphism group")
    stabilizer = []
    for T in group:
        if any((x-y)%2 for x,y in zip(matvec(T,a),a)):
            continue
        dual = matvec(transpose(inverse(T)),b)
        if any((x-y)%2 for x,y in zip(dual,b)):
            continue
        epsilon = symmetry_data(G,a,b,T)["epsilon"]
        require(epsilon == 0,"counterexample has a symmetry witness")
        stabilizer.append({"T":T,"epsilon":epsilon})
    require(stabilizer,"empty stabilizer")
    enum = doc["corroboration"]
    bound = integer(enum["norm4_bound"],"corroboration bound")
    require(enum["complete"] is True and bound >= 0,"incomplete corroboration")
    shells = box_shells(G,a,b,bound,max_nodes)
    require(shells == enum["shells"] and all(x["signed_coefficient"] == 0 for x in shells),
            "corroborating coefficients disagree")
    inv = inverse(G)
    from math import isqrt
    radius = [isqrt((28*inv[i][i]).numerator//(28*inv[i][i]).denominator) for i in range(4)]
    basis_bound = max(G[i][i] for i in range(4))
    vectors = box_vectors(G,basis_bound,max_nodes)
    return {"verdict":"verified_rank4_symmetry_converse_counterexample", "determinant":det4(G),
            "ldl_pivots":[str(x) for x in pivots],"identity_proof":"unimodular ambient coset bijection at all norms",
            "full_group_order":len(group),"stabilizer":stabilizer,
            "all_stabilizer_epsilon_zero":True,"bounded_cancellation_used_as_identity_proof":False,
            "basis_norm_bound":basis_bound,"box_radii_at_norm28":radius,
            "basis_image_vectors":{str(n):[list(v) for v in vectors if dot(v,matvec(G,v)) == n]
                                   for n in sorted(set(G[i][i] for i in range(4)))}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--replay",type=Path)
    args = parser.parse_args()
    require(not args.output.exists(),"choose a fresh output path")
    if args.replay:
        raw = args.replay.read_bytes()
        doc = {**replay(json.loads(raw)),"input_sha256":hashlib.sha256(raw).hexdigest()}
    else:
        doc = build()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(doc,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in doc.items() if k in ("schema","verdict","determinant","full_group_order","all_stabilizer_epsilon_zero")}))


if __name__ == "__main__":
    main()
