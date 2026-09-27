"""Exact finite regressions for the proved two-parameter later-shell lemma."""
import argparse
import hashlib
import json
from pathlib import Path

from exact_theta import complete_shells, integer, require
from verify_lattice_census import box_shells

PARAMETERS = [(3,3),(3,4),(3,5),(4,4),(5,7),(8,11),(20,31)]
A = [0,0,1,1]
PAIRS = [([0,1,1,1],-2),([1,0,0,0],2),([1,1,1,1],-2)]


def gram(s,t):
    integer(s,"s")
    integer(t,"t")
    require(s >= 3 and t >= 3,"lemma requires s,t>=3")
    return [[2,-1,-1,-1],[-1,2,0,0],[-1,0,s,0],[-1,0,0,t]]


def shell_claim(s,t,coefficient,shells):
    m = s+t
    at = {x["norm4"]:x["signed_coefficient"] for x in shells}
    first = next((x for x in shells if x["signed_coefficient"]),None)
    require(m in at and at[m] == 0,"initial shell missing or nonzero")
    require(first is not None and first["norm4"] == m+8
            and first["signed_coefficient"] == coefficient,"incorrect first nonzero shell")
    return {"cancelled_minimal_norm4":m,"first_nonzero_norm4":m+8,"coefficient":coefficient}


def build(max_nodes=1_000_000):
    cases = []
    for s,t in PARAMETERS:
        G = gram(s,t)
        for b,c in PAIRS:
            enumeration = complete_shells(G,A,b,s+t+8,max_nodes)
            cases.append({"s":s,"t":t,"G":G,"a":A,"b":b,
                          "enumeration":enumeration,"claim":shell_claim(s,t,c,enumeration["shells"])})
    return {"schema":"work7-rank4-later-shell-controls-v1","cases":cases}


def replay(doc,max_nodes=1_000_000):
    require(doc["schema"] == "work7-rank4-later-shell-controls-v1","wrong controls schema")
    require([(r["s"],r["t"],r["b"]) for r in doc["cases"]]
            == [(s,t,b) for s,t in PARAMETERS for b,c in PAIRS],"wrong case coverage")
    for row,(s,t,b,c) in zip(doc["cases"],[(s,t,b,c) for s,t in PARAMETERS for b,c in PAIRS]):
        G = gram(row["s"],row["t"])
        require(row["G"] == G and row["a"] == A and row["b"] == b,"wrong Gram or characteristic")
        enumeration = row["enumeration"]
        require(enumeration["complete"] is True
                and integer(enumeration["norm4_bound"],"bound") == s+t+8,"wrong complete bound")
        shells = box_shells(G,A,b,s+t+8,max_nodes)
        require(shells == enumeration["shells"],"changed complete shell histogram")
        require(row["claim"] == shell_claim(s,t,c,shells),"changed shell claim")
    return {"verdict":"verified","cases":len(doc["cases"]),"parameter_pairs":len(PARAMETERS),
            "vector_backend":"box","uniform_proof_performed":False}


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
    print(json.dumps({k:v for k,v in doc.items() if k != "cases"} | {"cases":len(doc["cases"]) if isinstance(doc.get("cases"),list) else doc["cases"]}))


if __name__ == "__main__":
    main()
