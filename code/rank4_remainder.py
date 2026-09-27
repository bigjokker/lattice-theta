"""Bounded orthogonal-summand and averaged-shell regressions, with box replay."""
import argparse
from collections import Counter, defaultdict
from itertools import combinations, product
import hashlib
import json
from math import gcd
from pathlib import Path
import time

from exact_theta import (EnumerationLimit, dot, integer, inverse, matrix, matvec,
                         multiply, require, transpose)
from lattice_census import sphere
from rank4_feasibility import BITS
from verify_lattice_census import box_vectors

CENSUS_SHA256 = "f3411b2ae36cda48b1d93b92c8edd16eddf659e7a3031e47f613b484507ef982"


def saturated_pair(v,w):
    index = 0
    for i,j in combinations(range(4),2):
        index = gcd(index, v[i]*w[j]-v[j]*w[i])
    return abs(index) == 1


def canonical_vectors(G,bound,backend,max_nodes):
    vectors = sphere(G,bound,max_nodes) if backend == "ldl" else box_vectors(G,bound,max_nodes)
    return sorted(v for v in vectors if any(v) and next(x for x in v if x) > 0
                  and gcd(*v) == 1)


def split_search(G,bound,backend,max_nodes):
    vectors = canonical_vectors(G,bound,backend,max_nodes)
    transforms = {v:matvec(G,v) for v in vectors}
    norms = {v:dot(v,transforms[v]) for v in vectors}
    for v in vectors:
        if all(x % norms[v] == 0 for x in transforms[v]):
            P = [[v[i]*transforms[v][j]//norms[v] for j in range(4)] for i in range(4)]
            return {"kind":"line","vectors":[list(v)],"projector":P}
    nodes = 0
    for v,w in combinations(vectors,2):
        nodes += 1
        if nodes > max_nodes:
            raise EnumerationLimit("orthogonal-plane search interrupted")
        if not saturated_pair(v,w):
            continue
        a,c,b = norms[v],norms[w],dot(v,transforms[w])
        D = a*c-b*b
        if D <= 0:
            continue
        left = [c*x-b*y for x,y in zip(transforms[v],transforms[w])]
        right = [a*y-b*x for x,y in zip(transforms[v],transforms[w])]
        if all(x % D == 0 for x in left+right):
            P = [[(v[i]*left[j]+w[i]*right[j])//D for j in range(4)] for i in range(4)]
            return {"kind":"plane","vectors":[list(v),list(w)],"projector":P}
    return {"kind":"indecomposable","enumeration_complete":True,
            "primitive_vectors_up_to_sign":len(vectors),"pairs_scanned":nodes}


def validate_projector(G,item):
    """Independent rational construction checks the stored integral projection."""
    vectors = item["vectors"]
    require(len(vectors) == (1 if item["kind"] == "line" else 2),"wrong projector rank")
    W = transpose(vectors)
    inner = multiply(multiply(transpose(W),G),W)
    P = multiply(multiply(multiply(W,inverse(inner)),transpose(W)),G)
    stored = matrix(item["projector"],"projector",4)
    require(P == stored and all(x.denominator == 1 for row in P for x in row),
            "incorrect integral projector")
    require(multiply(stored,stored) == stored
            and multiply(transpose(stored),G) == multiply(G,stored),"invalid orthogonal idempotent")
    require(sum(stored[i][i] for i in range(4)) == len(vectors),"wrong projection trace")


def direct_average_check(group,a,b,stored):
    """Check the six recorded regressions from f(T^-1 y), without label phases."""
    inverses = [[[int(x) for x in row] for row in inverse(T)] for T in group]
    expected = {tuple(item["residue"]):item["numerator"] for item in stored}
    actual = {}
    for r in product(range(4),repeat=4):
        numerator = 0
        for T in inverses:
            y = matvec(T,r)
            if all((x-z)%2 == 0 for x,z in zip(y,a)):
                numerator += (-1)**(dot([(x-z)//2 for x,z in zip(y,a)],b)%2)
        if numerator:
            actual[r] = numerator
    require(expected == actual,"direct residue average mismatch")


def averaged_regressions(row,backend,max_nodes,bound=32):
    G,group = row["G"],row["automorphisms"]
    residues = list(product(range(4),repeat=4))
    residue_index = {r:i for i,r in enumerate(residues)}
    bits_index = {tuple(b):i for i,b in enumerate(BITS)}
    walsh = [[(-1)**(dot(z,b)%2) for b in BITS] for z in BITS]
    residue_labels = [(bits_index[tuple(x%2 for x in r)],
                       bits_index[tuple(x//2 for x in r)]) for r in residues]
    actions = []
    for T in group:
        primal = [matvec(T,a) for a in BITS]
        dual = [matvec(transpose(inverse(T)),b) for b in BITS]
        actions.append((primal,dual))
    vectors = sphere(G,bound,max_nodes) if backend == "ldl" else box_vectors(G,bound,max_nodes)
    hist = defaultdict(Counter)
    for v in vectors:
        hist[dot(v,matvec(G,v))][residue_index[tuple(x%4 for x in v)]] += 1
    cancellations = []
    checked_pairs = 0
    for p in row["pairs"]:
        a,b = p["a"],p["b"]
        ai,bi = bits_index[tuple(a)],bits_index[tuple(b)]
        coefficients = Counter()
        for primal,dual in actions:
            image = primal[ai]
            ap = [x%2 for x in image]
            bp = [int(x)%2 for x in dual[bi]]
            sign = (-1)**(dot([(x-y)//2 for x,y in zip(image,ap)],bp)%2)
            coefficients[(bits_index[tuple(ap)],bits_index[tuple(bp)])] += sign
        coefficients = {k:v for k,v in coefficients.items() if v}
        has_negative = 1 in p["stabilizer"]["epsilon"]
        require((not coefficients) == has_negative,"Walsh average criterion disagreement")
        if dot(a,b)%2 or has_negative:
            continue
        checked_pairs += 1
        weights = [sum(coefficients.get((ai,bi),0)*walsh[zi][bi] for bi in range(16))
                   for ai,zi in residue_labels]
        require(any(weights),"nonzero Walsh expansion became zero function")
        shells = []
        first_support = None
        for norm,counts in sorted(hist.items()):
            if first_support is None and any(weights[i] for i in counts):
                first_support = norm
            numerator = sum(weights[i]*count for i,count in counts.items())
            require(numerator % len(group) == 0,"nonintegral averaged shell coefficient")
            shells.append((norm,numerator//len(group)))
        first = next(((n,c) for n,c in shells if c),None)
        require(p["verdict"] == "proved_nonzero" and first is not None
                and first == (p["evidence"]["norm4"],p["evidence"]["coefficient"]),
                "average does not match census nonzero")
        if dict(shells)[first_support] == 0:
            cancellation = {"a":a,"b":b,"stabilizer_order":len(p["stabilizer"]["indices"]),
                "first_support_norm4":first_support,"first_nonzero_norm4":first[0],
                "first_nonzero_coefficient":first[1],
                "average_numerators":[{"residue":list(r),"numerator":w} for r,w in zip(residues,weights) if w],
                "denominator":len(group),
                "support_shell_vectors":[{"v":list(v),"average_numerator":weights[residue_index[tuple(x%4 for x in v)]]}
                    for v in sorted(vectors) if dot(v,matvec(G,v)) == first_support
                    and weights[residue_index[tuple(x%4 for x in v)]]],
                "coefficients_through_bound":[{"norm4":n,"coefficient":c} for n,c in shells]}
            direct_average_check(group,a,b,cancellation["average_numerators"])
            cancellations.append(cancellation)
    return {"average_criterion_pairs_checked":256,"epsilon_trivial_even_pairs_checked":checked_pairs,
            "cancelled_first_support_pairs":cancellations}


def analyze(census,backend="ldl",max_nodes=3_000_000):
    require(census["schema"] == "work7-rank4-census-v1","wrong input")
    bound = integer(census["domain"]["determinant_bound"],"domain")
    require(bound == 24 and backend in ("ldl","box"),"wrong remainder domain/backend")
    integer(max_nodes,"limit")
    require(max_nodes > 0,"positive limit required")
    rows = []
    for index,row in enumerate(census["lattices"]):
        item = {"representative":index,"G":row["G"],"determinant":row["determinant"],
                "split":split_search(row["G"],bound,backend,max_nodes)}
        if item["split"]["kind"] != "indecomposable":
            validate_projector(row["G"],item["split"])
        else:
            item["group_order"] = len(row["automorphisms"])
            item["even_symmetry_zeros"] = sum(dot(p["a"],p["b"])%2 == 0 and p["verdict"] == "proved_zero" for p in row["pairs"])
            item["averages"] = averaged_regressions(row,backend,max_nodes)
        rows.append(item)
    kinds = Counter(r["split"]["kind"] for r in rows)
    return {"schema":"work7-rank4-remainder-v1","summand_norm_bound":bound,"shell_norm4_bound":32,
            "lattices":rows,"summary":{"lattice_classes":len(rows),"split_kinds":dict(kinds),
            "indecomposable_even_symmetry_zeros":sum(r.get("even_symmetry_zeros",0) for r in rows),
            "cancelled_first_support_pairs":sum(len(r.get("averages",{}).get("cancelled_first_support_pairs",[])) for r in rows)}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("census",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--replay",type=Path)
    parser.add_argument("--max-nodes",type=int,default=3_000_000)
    args = parser.parse_args()
    require(not args.output.exists(),"choose a fresh output path")
    raw = args.census.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == CENSUS_SHA256,
            "remainder analysis requires the accepted census bytes")
    start = time.perf_counter()
    # Interruptions raise before any result is saved: absence is never inferred
    # from a partially enumerated sphere or a partially scanned pair list.
    result = analyze(json.loads(raw),"box" if args.replay else "ldl",args.max_nodes)
    digest = hashlib.sha256(raw).hexdigest()
    if args.replay:
        stored_raw = args.replay.read_bytes()
        stored = json.loads(stored_raw)
        require(stored["source_sha256"] == digest,"wrong bound census input")
        require({k:v for k,v in stored.items() if k not in ("source_sha256","seconds")} == result,
                "remainder replay mismatch")
        result = {"verdict":"verified","source_sha256":digest,
                  "remainder_sha256":hashlib.sha256(stored_raw).hexdigest(),
                  "summary":result["summary"],"vector_backend":"box",
                  "orthogonal_projection_criterion_shared":True}
    else:
        result["source_sha256"] = digest
    result["seconds"] = time.perf_counter()-start
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k != "lattices"}))


if __name__ == "__main__":
    main()
