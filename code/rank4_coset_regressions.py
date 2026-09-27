"""Check the coupled-coset reformulation on two nonzeros and a symmetry zero."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path

from exact_theta import dot, integer, matvec, require, symmetry_data
from lattice_census import sphere
from rank4_feasibility import action_tables, admit, characteristic, isometries, modular_data
from verify_lattice_census import box_vectors
from verify_rank4_feasibility import box_isometries, modular_check

CASES = [
    ("cancelled first shell", [[2,-1,-1,-1],[-1,2,0,0],[-1,0,3,0],[-1,0,0,3]],
     [0,0,1,1],[0,1,1,1],14),
    ("noncancelled first shell", [[2,-1,-1,-1],[-1,3,0,0],[-1,0,4,0],[-1,0,0,4]],
     [0,0,1,1],[0,1,1,1],8),
    ("D4 even symmetry zero", [[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]],
     [0,0,1,1],[0,1,0,0],16),
]


def counts(G,a,b,bound,backend,max_nodes):
    vectors = sphere(G,bound,max_nodes) if backend == "ldl" else box_vectors(G,bound,max_nodes)
    hist = defaultdict(lambda:[0,0])
    for y in vectors:
        if any((v-w)%2 for v,w in zip(y,a)):
            continue
        x = [(v-w)//2 for v,w in zip(y,a)]
        hist[dot(y,matvec(G,y))][dot(x,b)%2] += 1
    return [{"norm4":N,"plus_count":v[0],"minus_count":v[1],"difference":v[0]-v[1]}
            for N,v in sorted(hist.items())]


def build(max_nodes=1_000_000):
    rows = []
    for name,G,a,b,bound in CASES:
        group = isometries(G,G,max_nodes)
        pair = characteristic(G,a,b,group,action_tables(group),bound,max_nodes)
        j = next(i for i,v in enumerate(b) if v)
        u = [int(i==j) for i in range(4)]
        row = {"name":name,"G":G,"determinant":admit(G),"a":a,"b":b,"u":u,
               "bound":bound,"modular":modular_data(G,a,b),"automorphisms":group,
               "pair":pair,"coset_counts":counts(G,a,b,bound,"ldl",max_nodes)}
        rows.append(row)
    return {"schema":"work7-rank4-coupled-coset-controls-v1","cases":rows}


def replay(doc,max_nodes=1_000_000):
    require(doc["schema"] == "work7-rank4-coupled-coset-controls-v1","wrong controls schema")
    require([r["name"] for r in doc["cases"]] == [r[0] for r in CASES],"wrong case coverage")
    reports = []
    for row,(name,G,a,b,bound) in zip(doc["cases"],CASES):
        require((row["G"],row["a"],row["b"]) == (G,a,b),"wrong Gram or labels")
        require(integer(row["bound"],"bound") == bound
                and integer(row["determinant"],"determinant") == admit(G),"wrong bounds or determinant")
        u = row["u"]
        require(len(u) == 4 and all(type(x) is int for x in u) and dot(u,b)%2 == 1,"wrong complementary vector")
        modular_check(G,a,b,row["modular"])
        group = box_isometries(G,G,max_nodes)
        require(group == row["automorphisms"],"wrong full group")
        actual = counts(G,a,b,bound,"box",max_nodes)
        require(actual == row["coset_counts"],"changed unsigned coset counts")
        p = row["pair"]
        # Recompute each stabilizer/sign directly; do not trust the label table.
        from exact_theta import inverse,transpose
        indices,eps = [],[]
        for i,T in enumerate(group):
            if any((x-y)%2 for x,y in zip(matvec(T,a),a)):
                continue
            dual = matvec(transpose(inverse(T)),b)
            if any((x-y)%2 for x,y in zip(dual,b)):
                continue
            indices.append(i)
            eps.append(int(dot(a,[(x-y)/2 for x,y in zip(dual,b)]))%2)
        require(p["a"] == a and p["b"] == b
                and p["stabilizer"] == {"indices":indices,"epsilon":eps},"wrong stabilizer")
        first = next((r for r in actual if r["difference"]),None)
        if name == "D4 even symmetry zero":
            T = p["evidence"]["T"]
            require(p["verdict"] == "proved_zero" and 1 in eps and T in group
                    and symmetry_data(G,a,b,T)["epsilon"] == 1,"missing symmetry proof")
            # Check both coset swaps on their offsets and their difference lattice:
            # T L0=L0 follows from dual b=b mod 2; parity below locates the offsets.
            for offset,expected_parity in [(a,1),([x+2*y for x,y in zip(a,u)],0)]:
                image = matvec(T,offset)
                require(all((x-y)%2 == 0 for x,y in zip(image,a)),"changed support coset")
                require(dot([(x-y)//2 for x,y in zip(image,a)],b)%2 == expected_parity,
                        "witness does not swap both cosets")
            require(first is None,"symmetry zero has unequal counts")
            reports.append({"name":name,"verdict":"proved_equal_by_coset_swap","group_order":len(group)})
        else:
            from fractions import Fraction
            enumeration = p["enumeration"]
            expected_shells = [{"norm4":r["norm4"],"norm":str(Fraction(r["norm4"],4)),
                                "vector_count":r["plus_count"]+r["minus_count"],
                                "signed_coefficient":r["difference"]} for r in actual]
            require(enumeration["complete"] is True and enumeration["norm4_bound"] == bound
                    and enumeration["shells"] == expected_shells,"changed signed shell histogram")
            require(p["verdict"] == "proved_nonzero" and 1 not in eps and first is not None
                    and p["evidence"] == {"kind":"shell","norm4":first["norm4"],"coefficient":first["difference"]},
                    "wrong nonzero proof")
            reports.append({"name":name,"verdict":"proved_unequal","group_order":len(group),
                            "first_difference_norm4":first["norm4"],"difference":first["difference"]})
    return {"verdict":"verified","cases":reports,"unsigned_count_backend":"box"}


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
    print(json.dumps({k:v for k,v in doc.items() if k != "cases"}))


if __name__ == "__main__":
    main()
