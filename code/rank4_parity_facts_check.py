"""Control facts for item 4 (rank-four parity property), prepared 2026-10-20
by the brief-writing session for Fable-rank4-parity-conjecture.md.

Section 1: the 27 indecomposable census classes (det <= 24) with evenness,
det mod 4, group order, parity-property status and number of vanishing
even characteristics, read from the accepted census pack.
Section 2: the 123 Task B representatives of the Eisenstein family F
(all fail P): evenness, det mod 4, |Aut(L)|, decomposability, level.
Section 3: exact nonzero proofs for all 105 even characteristics (a, b nonzero) of the
Conway-Sloane lattices L^+(1,7,13,19) and L^-(1,7,13,19) (even, det 1729):
both have the parity property.  Control: the det-1156 lattice has one zero.
Run from the workspace root (about 20 s):
    python code/rank4_parity_facts_check.py > reports/rank4-parity/control-facts.txt
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from itertools import product
from collections import Counter
from rank4_feasibility import isometries, determinant
from exact_theta import complete_shells, multiply, identity
I = identity(4)
def decomposable(group):
    for T in group:
        if T == I or T == [[-x for x in r] for r in I]: continue
        if multiply(T, T) != I: continue
        if all((x + y) % 2 == 0 for rT, rI in zip(T, I) for x, y in zip(rT, rI)): return True
    return False

print("=== Section 1: census, 27 indecomposable classes ===")
d = json.load(open("data/census-rank4/det24-bound32.json"))
rows = []
for L in d["lattices"]:
    G = L["G"]; det = L["determinant"]; even = all(G[i][i] % 2 == 0 for i in range(4))
    dec = decomposable(L["automorphisms"])
    zeros = [p for p in L["pairs"] if p.get("verdict", p.get("status")) == "proved_zero" and any(p["a"]) and any(p["b"])
             and sum(x*y for x, y in zip(p["a"], p["b"])) % 2 == 0]
    rows.append((det, det % 4, even, dec, len(L["automorphisms"]), len(zeros) == 0, len(zeros)))
indec = [r for r in rows if not r[3]]
print("classes", len(rows), "indecomposable", len(indec), "every decomposable class fails P:", all(not r[5] for r in rows if r[3]))
print("det det%4 parity |Aut| P #even-zeros")
for r in sorted(indec): print("  ", r[0], r[1], "even" if r[2] else "odd", r[4], "holds" if r[5] else "fails", r[6])
print("holding (parity, det%4):", dict(Counter(("even" if r[2] else "odd", r[1]) for r in indec if r[5])))
print("failing (parity, det%4):", dict(Counter(("even" if r[2] else "odd", r[1]) for r in indec if not r[5])))

print("=== Section 2: Eisenstein family, Task B representatives (all fail P) ===")
d = json.load(open("reports/rank4-eisenstein/task-b.json"))
reps, members = d["representatives"], d["members"]
rows = []
for r in reps:
    m = members[r["member_index"]]; G = m["L_reduced"]; det = determinant(G)
    group = isometries(G, G, 5_000_000)
    rows.append((r["D_prime"], det, det % 4, all(G[i][i] % 2 == 0 for i in range(4)), len(group), decomposable(group), r["level"], G, m["a"], m["b"]))
indec = [x for x in rows if not x[5]]
print("representatives", len(rows), "indecomposable", len(indec))
print("indecomposable, det not divisible by 4 (all odd lattices):")
for x in sorted(indec):
    if x[2] != 0: print("   D'", x[0], "det", x[1], "even" if x[3] else "odd", "|Aut|", x[4], "level", x[6], "G", x[7], "a", x[8], "b", x[9])
print("indecomposable even representatives, det%4 tally:", dict(Counter(x[2] for x in indec if x[3])))
print("indecomposable odd-lattice representatives, det%4 tally:", dict(Counter(x[2] for x in indec if not x[3])))

print("=== Section 3: exact parity-property check of the Conway-Sloane det-1729 lattices ===")
Lp = [[4,-1,-4,-4],[-1,8,0,-4],[-4,0,12,1],[-4,-4,1,16]]
Lm = [[4,2,2,5],[2,8,-3,2],[2,-3,12,-2],[5,2,-2,16]]
G1156 = [[4,-2,1,-2],[-2,8,2,-2],[1,2,6,-4],[-2,-2,-4,12]]
def probe(G, name, bound=600):
    t0 = time.time(); survivors = []; tested = 0
    for a in product((0,1), repeat=4):
        if not any(a): continue
        for b in product((0,1), repeat=4):
            if not any(b) or sum(x*y for x,y in zip(a,b)) % 2: continue
            tested += 1
            sh = complete_shells(G, list(a), list(b), bound, 50_000_000)
            if all(s["signed_coefficient"] == 0 for s in sh["shells"]): survivors.append((a, b))
    print(f"{name}: det {determinant(G)}, even={all(G[i][i]%2==0 for i in range(4))}, |Aut|={len(isometries(G,G,5_000_000))};"
          f" even characteristics without a nonzero signed coefficient through norm4 {bound}: {len(survivors)} of {tested} ({time.time()-t0:.1f}s)")
    for s in survivors: print("    unresolved a,b =", s)
probe(Lp, "CS L+(1,7,13,19)"); probe(Lm, "CS L-(1,7,13,19)"); probe(G1156, "det-1156 L (control, expect exactly one)")
