"""Rank-four candidate sizing and exact controls; not an isometry-class census."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from math import ceil, floor, lcm, prod
from pathlib import Path
import time

from exact_theta import (EnumerationLimit, complete_shells, dot, identity, integer,
                         inverse, ldl, matrix, matvec, multiply, require, transpose, vector)
from lattice_census import sphere

RATIO = Fraction(3, 4)
BITS = [list(v) for v in product((0, 1), repeat=4)]


def determinant(A):
    work = [[Fraction(x) for x in row] for row in A]
    sign, result = 1, Fraction(1)
    for j in range(len(A)):
        pivot = next((i for i in range(j, len(A)) if work[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            work[j], work[pivot] = work[pivot], work[j]
            sign = -sign
        value = work[j][j]
        result *= value
        for i in range(j+1, len(A)):
            ratio = work[i][j] / value
            for k in range(j+1, len(A)):
                work[i][k] -= ratio * work[j][k]
    result *= sign
    require(result.denominator == 1, "nonintegral determinant")
    return int(result)


def admit(G):
    matrix(G, "G", 4)
    _, d = ldl(G)
    delta=prod(d)
    require(delta.denominator==1,"nonintegral Gram determinant")
    return int(delta)


def candidates(H=24, max_nodes=1_000_000):
    """All admitted reduced-basis candidates, via recursive rational intervals."""
    integer(H, "H")
    integer(max_nodes, "max_nodes")
    require(1 <= H <= 24 and max_nodes > 0, "invalid feasibility domain or limit")
    found, nodes = [], 0

    def tick():
        nonlocal nodes
        nodes += 1
        if nodes > max_nodes:
            raise EnumerationLimit("candidate generation interrupted; partial domain discarded")

    def extend(G):
        size = len(G)
        mu, pivots = ldl(G)
        prefix = prod(pivots)
        remaining = 4-size

        def crosses(j, entries, coeffs):
            if j == size:
                q = sum(c*c*d for c, d in zip(coeffs, pivots))
                lower = max(G[0][0], ceil(q + RATIO*pivots[-1]))
                if size == 3:
                    diagonals = []
                    for delta in range(1, H+1):
                        value = Fraction(delta, 1)/prefix + q
                        if value.denominator == 1 and value >= lower:
                            diagonals.append(int(value))
                else:
                    diagonals = []
                    diagonal = lower
                    while prefix*(Fraction(diagonal)-q)**remaining * RATIO**(remaining*(remaining-1)//2) <= H:
                        tick()
                        diagonals.append(diagonal)
                        diagonal += 1
                for diagonal in diagonals:
                    tick()
                    row = [list(old)+[entry] for old, entry in zip(G, entries)]
                    row.append(entries+[diagonal])
                    if size == 3:
                        found.append(row)
                    else:
                        extend(row)
                return
            center = sum(mu[j][k]*coeffs[k]*pivots[k] for k in range(j))
            for entry in range(ceil(center-pivots[j]/2), floor(center+pivots[j]/2)+1):
                tick()
                crosses(j+1, entries+[entry], coeffs+[(entry-center)/pivots[j]])

        crosses(0, [], [])

    a = 1
    while RATIO**6 * a**4 <= H:
        tick()
        extend([[a]])
        a += 1
    return sorted(found)


def isometries(G, H, max_nodes=1_000_000):
    """Complete LDL spheres and recursively pruned basis-image search."""
    admit(G)
    admit(H)
    integer(max_nodes, "max_nodes")
    require(max_nodes > 0, "positive group node limit required")
    vectors = sphere(H, max(G[i][i] for i in range(4)), max_nodes)
    images = [[v for v in vectors if dot(v, matvec(H, v)) == G[i][i]] for i in range(4)]
    Himages = {v: matvec(H, v) for v in vectors}
    matches, nodes = [], 0

    def extend(columns):
        nonlocal nodes
        j = len(columns)
        if j == 4:
            T = transpose(columns)
            if abs(determinant(T)) == 1:
                matches.append(T)
            return
        for v in images[j]:
            nodes += 1
            if nodes > max_nodes:
                raise EnumerationLimit("rank-four full-group search interrupted")
            if all(dot(u, Himages[v]) == G[i][j] for i, u in enumerate(columns)):
                extend(columns+[v])

    extend([])
    return sorted(matches)


def action_tables(group):
    tables = []
    for T in group:
        dual = transpose(inverse(T))
        require(all(v.denominator == 1 for row in dual for v in row), "nonintegral dual action")
        tables.append(([[v % 2 for v in matvec(T, a)] for a in BITS],
                       [[int(v) for v in matvec(dual, b)] for b in BITS]))
    return tables


def characteristic(G, a, b, group, tables, bound=32, max_nodes=1_000_000):
    ai, bi = BITS.index(a), BITS.index(b)
    indices, signs = [], []
    for i, (primal, dual) in enumerate(tables):
        if primal[ai] == a and all((v-w) % 2 == 0 for v, w in zip(dual[bi], b)):
            indices.append(i)
            signs.append(dot(a, [(v-w)//2 for v, w in zip(dual[bi], b)]) % 2)
    row = {"a": a, "b": b, "stabilizer": {"indices": indices, "epsilon": signs}}
    if 1 in signs:
        i = indices[signs.index(1)]
        return {**row, "verdict": "proved_zero", "evidence": {"kind": "symmetry", "T": group[i]}}
    try:
        enumeration = complete_shells(G, a, b, bound, max_nodes)
    except EnumerationLimit as error:
        return {**row, "verdict": "unresolved", "enumeration_complete": False, "reason": str(error)}
    nonzero = next((s for s in enumeration["shells"] if s["signed_coefficient"]), None)
    if nonzero is None:
        return {**row, "verdict": "unresolved", "enumeration": enumeration}
    return {**row, "verdict": "proved_nonzero", "enumeration": enumeration,
            "evidence": {"kind": "shell", "norm4": nonzero["norm4"], "coefficient": nonzero["signed_coefficient"]}}


def controls(bound=32, max_nodes=1_000_000):
    A4 = [[2 if i == j else -1 if abs(i-j) == 1 else 0 for j in range(4)] for i in range(4)]
    D4 = [[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]]
    mixed = [[2,-1,0,0],[-1,2,-1,0],[0,-1,2,0],[0,0,0,1]]
    inputs = [("I4", identity(4)), ("2I4", [[2*int(i==j) for j in range(4)] for i in range(4)]),
              ("A4", A4), ("D4", D4), ("A3 plus line", mixed)]
    rows = []
    for name, G in inputs:
        group = isometries(G, G, max_nodes)
        tables = action_tables(group)
        pairs = [characteristic(G, a, b, group, tables, bound, max_nodes) for a in BITS for b in BITS]
        even = Counter(p["verdict"] for p in pairs if dot(p["a"], p["b"]) % 2 == 0)
        rows.append({"name": name, "G": G, "determinant": admit(G), "automorphisms": group,
                     "pairs": pairs, "summary": {"group_order": len(group), "even_pairs": 136,
                     "odd_pairs": 120, "even_verdicts": dict(even), "all_pairs": len(pairs)}})
    return {"schema": "work7-rank4-feasibility-controls-v1", "domain": "five named controls; not a census",
            "run": {"norm4_bound": bound, "max_nodes": max_nodes}, "lattices": rows}


def modular_data(G, a, b):
    admit(G)
    vector(a, "a", 4)
    vector(b, "b", 4)
    require(all(v in (0,1) for v in a+b) and any(b) and dot(a,b)%2 == 0,
            "binary columns, nonzero b and even pairing required")
    j = next(i for i in range(4) if b[i])
    K = identity(4)
    K[j][j] = 2
    for i in range(4):
        if i != j:
            K[j][i] = -b[i]
    G0 = multiply(multiply(transpose(K), G), K)
    inv = inverse(G0)
    N0 = lcm(*(v.denominator for row in inv for v in row))
    M, remaining, p, index = 16*N0, 16*N0, 2, 1
    factors = []
    while p*p <= remaining:
        e = 0
        while remaining % p == 0:
            e += 1
            remaining //= p
        if e:
            factors.append([p,e])
            index *= p**(e-1)*(p+1)
        p += 1
    if remaining > 1:
        factors.append([remaining,1])
        index *= remaining+1
    require(index % 6 == 0, "rank-four cutoff integrality failure")
    return {"index_two_basis": K, "G0": G0, "d0": determinant(G0), "N0": N0,
            "unsigned_level": M, "level_factors": factors, "index": index, "cutoff": index//6,
            "character_discriminant": 4*determinant(G0), "weight": "2",
            "coefficient_index": "norm4", "theorem": "rank4-square-complex-identity-v1"}


def cutoff_verdict(shells, bound, B):
    first = next((s for s in shells if s["signed_coefficient"]), None)
    if first:
        return {"verdict": "proved_nonzero", "nonzero_norm4": first["norm4"], "coefficient": first["signed_coefficient"]}
    return {"verdict": "proved_zero_by_modular_cutoff" if bound >= B else "unresolved"}


def cutoff_build(G, a, b, bound=None, max_nodes=1_000_000):
    item = modular_data(G, a, b)
    bound = item["cutoff"] if bound is None else integer(bound, "bound")
    integer(max_nodes, "max_nodes")
    require(bound >= 0 and max_nodes > 0, "invalid cutoff resource limits")
    doc = {"schema": "work7-rank4-cutoff-v1", "G": G, "a": a, "b": b, "modular": item,
           "run": {"norm4_bound": bound, "max_nodes": max_nodes}}
    try:
        enumeration = complete_shells(G,a,b,bound,max_nodes)
        doc["evidence"] = {"enumeration_complete": True, "enumeration": enumeration,
                           **cutoff_verdict(enumeration["shells"],bound,item["cutoff"])}
    except EnumerationLimit as error:
        doc["evidence"] = {"enumeration_complete": False, "verdict": "unresolved", "reason": str(error)}
    return doc


def cutoff_controls():
    G = identity(4)
    a, b = [1,1,0,0], [1,1,0,0]
    zero = cutoff_build(G,a,b)
    cases = [("complete modular zero", zero),
             ("endpoint missing", cutoff_build(G,a,b,zero["modular"]["cutoff"]-1)),
             ("interrupted", cutoff_build(G,a,b,max_nodes=1)),
             ("constant coefficient", cutoff_build(G,[0]*4,[1,0,0,0],bound=0))]
    for name, H in [("I4 nonzero", G),
                    ("2I4 nonzero", [[2*int(i==j) for j in range(4)] for i in range(4)]),
                    ("A4 nonzero", [[2 if i==j else -1 if abs(i-j)==1 else 0 for j in range(4)] for i in range(4)]),
                    ("D4 nonzero", [[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]])]:
        cases.append((name,cutoff_build(H,[0,1,0,0],[1,0,0,0])))
    trap=[[2,-1,-1,0],[-1,4,-1,0],[-1,-1,4,0],[0,0,0,1]]
    for name,bound in (("cancelled first shell",6),("later nonzero",10)):
        cases.append((name,cutoff_build(trap,[0,1,1,0],[1,0,0,0],bound=bound)))
    return {"schema": "work7-rank4-cutoff-controls-v1", "certificates": [{"name": name, **cert} for name, cert in cases]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    def save(name, doc):
        (args.output_dir/name).write_text(json.dumps(doc,indent=2)+"\n",encoding="utf-8")
    sizes = []
    for H in (1,2,4,8,12,16,24):
        start = time.perf_counter()
        forms = candidates(H)
        sizes.append({"determinant_bound": H, "candidate_bases": len(forms), "seconds": time.perf_counter()-start})
    save("candidates-det24.json", {"schema": "work7-rank4-candidate-bases-v1", "determinant_bound": 24,
         "generation_complete": True, "isometry_deduplication_performed": False, "candidate_bases": forms})
    save("candidate-sizing.json", {"schema": "work7-rank4-candidate-sizing-v1", "isometry_deduplication_performed": False, "rows": sizes})
    save("controls.json", controls())
    save("cutoff-controls.json", cutoff_controls())
    print(json.dumps({"candidate_sizes": sizes, "control_lattices": 5, "full_census_started": False}))


if __name__ == "__main__":
    main()
