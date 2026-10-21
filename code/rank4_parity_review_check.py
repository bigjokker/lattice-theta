"""Reviewer's independent replication of the headline claims of
RANK4-PARITY-CONJECTURE.md (review record: RANK4-PARITY-REVIEW.md).
Uses only the accepted modules rank4_feasibility and exact_theta; nothing from
rank4_parity_search.py or verify_rank4_parity.py.  Every decision below is
exact: a zero is an epsilon-one element of the characteristic stabilizer
(exact_theta.symmetry_data), a nonzero is a nonzero signed coefficient of a
complete shell (exact_theta.complete_shells).
Run from the workspace root:
    python code/rank4_parity_review_check.py > reports/rank4-parity/review-independent-check.txt
"""
import sys, os, json, time
from itertools import product
from collections import Counter
from fractions import Fraction
from math import gcd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rank4_feasibility import isometries, modular_data, determinant, cutoff_build
from exact_theta import complete_shells, symmetry_data, multiply, transpose, matvec, dot, inverse, identity

BITS = [list(v) for v in product((0, 1), repeat=4)]
EVEN = [(a, b) for a in BITS if any(a) for b in BITS if any(b) and dot(a, b) % 2 == 0]
I4 = identity(4)
MI4 = [[-x for x in r] for r in I4]
T0 = time.time()


def order(T):
    P, k = T, 1
    I = identity(len(T))
    while P != I:
        P = multiply(P, T)
        k += 1
    return k


def is_even(G):
    return all(G[i][i] % 2 == 0 for i in range(len(G)))


def stabilizer_eps(G, a, b, group):
    out = []
    for T in group:
        dual = transpose([[int(x) for x in r] for r in inverse(T)])
        if all((x - y) % 2 == 0 for x, y in zip(matvec(T, a), a)) and \
           all((x - y) % 2 == 0 for x, y in zip(matvec(dual, b), b)):
            out.append((T, symmetry_data(G, a, b, T)["epsilon"]))
    return out


def decide(G, group, a, b, ladder=(100, 400, 1600), max_nodes=20_000_000):
    st = stabilizer_eps(G, a, b, group)
    wit = [T for T, e in st if e == 1]
    if wit:
        return ("zero_by_symmetry", {"witness_orders": sorted({order(T) for T in wit}),
                                     "stabilizer": len(st), "witnesses": len(wit)})
    for bound in ladder:
        sh = complete_shells(G, a, b, bound, max_nodes)["shells"]
        nz = next((s for s in sh if s["signed_coefficient"]), None)
        if nz:
            return ("nonzero", {"norm4": nz["norm4"], "coefficient": nz["signed_coefficient"], "stabilizer": len(st)})
    return ("undecided", {"stabilizer": len(st)})


def decide_all(G, group, ladder=(100, 400, 1600)):
    return {(tuple(a), tuple(b)): decide(G, group, a, b, ladder) for a, b in EVEN}


def summarise(res):
    c = Counter(v[0] for v in res.values())
    zeros = [(k, v) for k, v in res.items() if v[0] == "zero_by_symmetry"]
    firsts = [v[1]["norm4"] for v in res.values() if v[0] == "nonzero"]
    wo = Counter(o for _, v in zeros for o in v[1]["witness_orders"])
    return dict(c), zeros, (max(firsts) if firsts else None), dict(wo)


def splitting_involutions(group):
    return [T for T in group if T not in (I4, MI4) and multiply(T, T) == I4
            and all((x + y) % 2 == 0 for rT, rI in zip(T, I4) for x, y in zip(rT, rI))]


def eichler_components(G):
    """Reviewer's own implementation: indecomposable vectors of norm <= max G_ii, joined by nonzero inner product."""
    R = max(G[i][i] for i in range(4))
    vecs = []
    complete_shells(G, [0] * 4, [0] * 4, 4 * R, 20_000_000, lambda y, q: vecs.append(tuple(v // 2 for v in y)))
    vecs = [v for v in vecs if any(v)]
    norm = {v: dot(v, matvec(G, v)) for v in vecs}
    indec = []
    for v in vecs:
        if v < tuple(-x for x in v):
            continue
        Gv = matvec(G, v)
        if not any(norm[u] < norm[v] and dot(u, Gv) == norm[u] for u in vecs):
            indec.append(v)
    parent = list(range(len(indec)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for i in range(len(indec)):
        Gi = matvec(G, indec[i])
        for j in range(i + 1, len(indec)):
            if dot(indec[j], Gi) != 0:
                parent[find(i)] = find(j)
    return len({find(i) for i in range(len(indec))}), len(indec)


def decomposable(G, group):
    inv = splitting_involutions(group)
    comps, nind = eichler_components(G)
    assert (len(inv) > 0) == (comps >= 2), "decomposability criteria disagree"
    return len(inv) > 0, comps, nind


def c_of(G, b):
    sols = [c for c in BITS if all((x - y) % 2 == 0 for x, y in zip(matvec(G, c), b))]
    return sols


def singular_count(G):
    return sum(1 for x in BITS if any(x) and (dot(x, matvec(G, x)) // 2) % 2 == 0)


def fdet(M):
    M = [[Fraction(x) for x in r] for r in M]
    n, d = len(M), Fraction(1)
    for i in range(n):
        p = next((r for r in range(i, n) if M[r][i] != 0), None)
        if p is None:
            return Fraction(0)
        if p != i:
            M[i], M[p] = M[p], M[i]
            d = -d
        d *= M[i][i]
        for r in range(i + 1, n):
            f = M[r][i] / M[i][i]
            M[r] = [x - f * y for x, y in zip(M[r], M[i])]
    return d


def report(name, G, expect_zero=None, ladder=(100, 400, 1600)):
    group = isometries(G, G, 5_000_000)
    dec, comps, nind = decomposable(G, group)
    res = decide_all(G, group, ladder)
    verd, zeros, maxfirst, wo = summarise(res)
    status = "holds" if verd.get("nonzero", 0) == 105 else ("fails" if zeros else "undecided")
    print(f"{name}: det {determinant(G)} {'even' if is_even(G) else 'odd'} |Aut| {len(group)} "
          f"decomposable {dec} (Eichler components {comps}) P {status} verdicts {verd} witness orders {wo} "
          f"max first nonzero norm4 {maxfirst}")
    return group, res, status, zeros


print("# Reviewer's independent check of RANK4-PARITY-CONJECTURE.md (accepted modules only)")
print()
print("## 1. The (C1) counterexample of determinant 121")
G = [[4, -2, -1, 0], [-2, 4, 0, 1], [-1, 0, 4, -2], [0, 1, -2, 4]]
print("det", determinant(G), "even", is_even(G), "gcd of entries", gcd(*[abs(x) for r in G for x in r]),
      "det mod 8", determinant(G) % 8, "nonzero singular classes", singular_count(G), "(9 = plus type)")
group = isometries(G, G, 5_000_000)
print("|Aut(L)| =", len(group), " element orders", dict(sorted(Counter(order(T) for T in group).items())))
sq = [T for T in group if multiply(T, T) == MI4]
print("elements with g^2 = -I:", len(sq))
g = sq[0]
print("g =", g)
cols = transpose([[x + y for x, y in zip(r1, r2)] for r1, r2 in zip(g, I4)])
Pi = sorted({tuple(sum(k * c[i] for k, c in zip(coeffs, cols)) % 2 for i in range(4))
             for coeffs in product((0, 1), repeat=4)} - {(0, 0, 0, 0)})
print("plane Pi = (g+I)L/2L nonzero classes:", Pi,
      " Q=0 on all:", all((dot(list(c), matvec(G, list(c))) // 2) % 2 == 0 for c in Pi),
      " B=0 pairwise:", all(dot(list(x), matvec(G, list(y))) % 2 == 0 for x in Pi for y in Pi))
dec, comps, nind = decomposable(G, group)
print("splitting involutions:", len(splitting_involutions(group)), " Eichler components:", comps,
      " indecomposable vectors up to sign:", nind, " -> decomposable:", dec)
res = decide_all(G, group)
verd, zeros, maxfirst, wo = summarise(res)
print("105 even characteristics:", verd, " witness orders of zeros:", wo, " max first nonzero norm4:", maxfirst)
for (a, b), v in zeros:
    cs = c_of(G, list(b))
    print("   zero a", list(a), "b", list(b), "c", cs, " a,c in Pi, distinct:",
          a in Pi and tuple(cs[0]) in Pi and tuple(cs[0]) != a, " stabilizer", v[1]["stabilizer"], "witnesses", v[1]["witnesses"])
pairs_on_plane = [(a, c) for a in Pi for c in Pi if a != c]
print("ordered pairs of distinct nonzero classes of Pi:", len(pairs_on_plane),
      " all zeros:", all(res[(a, tuple(x % 2 for x in matvec(G, list(c))))][0] == "zero_by_symmetry" for a, c in pairs_on_plane),
      " epsilon(g)=1 on all:", all(symmetry_data(G, list(a), [x % 2 for x in matvec(G, list(c))], g)["epsilon"] == 1 for a, c in pairs_on_plane))
a, b = [0, 1, 1, 0], [0, 1, 1, 0]
md = modular_data(G, a, b)
print("modular data for a=b=(0,1,1,0): d0", md["d0"], "N0", md["N0"], "M", md["unsigned_level"], "mu", md["index"], "B", md["cutoff"])
doc = cutoff_build(G, a, b, max_nodes=20_000_000)
ev = doc["evidence"]
print("cutoff enumeration complete:", ev["enumeration_complete"], "verdict", ev["verdict"], "shells", len(ev["enumeration"]["shells"]),
      "vectors", sum(s["vector_count"] for s in ev["enumeration"]["shells"]))
sh2 = complete_shells(G, a, b, 2000, 50_000_000)
print("(extra) through norm4 2000: vectors", sum(s["vector_count"] for s in sh2["shells"]), "nonzero signed", sum(1 for s in sh2["shells"] if s["signed_coefficient"]))
K = [[1, 0, 0, 0], [0, 2, -1, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
G0 = multiply(multiply(transpose(K), G), K)
print("kernel basis K, G0 =", G0, "det", determinant(G0), " all columns of K in the sign kernel:", all(dot(c, b) % 2 == 0 for c in transpose(K)))
group0 = isometries(G0, G0, 5_000_000)
Kinv = inverse(K)
v1 = [int(x) % 2 for x in matvec(Kinv, a)]
v2 = [int(x) % 2 for x in matvec(Kinv, [a[0], a[1] + 2, a[2], a[3]])]
mapping = [T for T in group0 if [x % 2 for x in matvec(T, v1)] == v2]
swapping = [T for T in mapping if [x % 2 for x in matvec(T, v2)] == v1]
print("|Aut(L0)| =", len(group0), " v1", v1, "v2", v2, " mapping v1->v2:", len(mapping), " swapping:", len(swapping))
H = [[4, 0, 2, -1], [0, 4, 1, 2], [2, 1, 4, 0], [-1, 2, 0, 4]]
print("trace form Re h of Z[i]-Gram [[4,2+i],[2-i,4]] on (e1, i e1, e2, i e2):", H, "det", determinant(H),
      " isometries to G:", len(isometries(H, G, 5_000_000)))
Ji = [[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]]
print("multiplication by i on the trace form preserves H:", multiply(multiply(transpose(Ji), H), Ji) == H, " squares to -I:", multiply(Ji, Ji) == MI4)
print("elapsed %.1f s" % (time.time() - T0))
print()

print("## 2. The smallest odd-lattice failure, determinant 49, and the family G_k")
G49 = [[3, -1, -1, 0], [-1, 3, 0, -1], [-1, 0, 3, 1], [0, -1, 1, 3]]
group49, res49, status49, zeros49 = report("det-49 lattice", G49)
W = [[0, 0, 0, -1], [0, -1, 0, 0], [0, 0, 1, 0], [-1, 0, 0, 0]]
print("   quoted witness on a=b=(0,1,1,0): in Aut:", W in group49, " order", order(W), " epsilon", symmetry_data(G49, [0, 1, 1, 0], [0, 1, 1, 0], W)["epsilon"])
print("   zeros:", [(list(a), list(b), v[1]["witness_orders"]) for (a, b), v in zeros49])
for k in (4, 5, 7, 8):
    Gk = [[3, -1, -1, 0], [-1, 3, 0, -1], [-1, 0, 3, 1], [0, -1, 1, k]]
    report(f"G_{k}", Gk)
print("elapsed %.1f s" % (time.time() - T0))
print()

print("## 3. The eighteen isospectral-pair members (reports/rank4-level4/task-a.json)")
pairs = json.load(open("reports/rank4-level4/task-a.json"))["pairs"]
for p in pairs:
    for which in ("G1", "G2"):
        Gm = p[which]
        grp, resm, st, zs = report(f"{p['name']} {which}", Gm)
        if is_even(Gm):
            cnt = 0
            for (a, b), v in resm.items():
                cs = c_of(Gm, list(b))
                if cs and (dot(list(a), matvec(Gm, list(a))) // 2) % 2 == 0 and all((dot(c, matvec(Gm, c)) // 2) % 2 == 0 for c in cs) and list(a) not in cs:
                    cnt += 1
            print(f"   type: {'plus' if singular_count(Gm) == 9 else 'minus' if singular_count(Gm) == 5 else 'degenerate'}, totally singular ordered pairs (a,c): {cnt}")
        if st == "fails" and is_even(Gm) is False:
            inv = splitting_involutions(grp)
            print("   splitting involutions:", len(inv), " zeros:", len(zs))
print("elapsed %.1f s" % (time.time() - T0))
print()

print("## 4. The eight certified level-4 lattices: the other 104 characteristics")
for path in ("data/rank4-level4/det1729-level4.json", "data/rank4-level4/other-level4-certificates.json"):
    for cert in json.load(open(path))["certificates"]:
        Gl, al, bl = cert["G"], cert["a"], cert["b"]
        grp = isometries(Gl, Gl, 5_000_000)
        dec, comps, nind = decomposable(Gl, grp)
        resl = decide_all(Gl, grp, ladder=(100, 400))
        verd, zeros, maxfirst, wo = summarise(resl)
        und = [k for k, v in resl.items() if v[0] == "undecided"]
        print(f"level-4 det {determinant(Gl)} {'even' if is_even(Gl) else 'odd'} |Aut| {len(grp)} decomposable {dec}: verdicts {verd}; "
              f"undecided through norm4 400: {[(list(a), list(b)) for a, b in und]}; certified zero is {(al, bl)}: "
              f"{'the only undecided one' if und == [(tuple(al), tuple(bl))] else 'MISMATCH'}; max first nonzero {maxfirst}")
print("elapsed %.1f s" % (time.time() - T0))
print()

print("## 5. The E8 witness (reports/rank4-parity/e8-check.json)")
e8 = json.load(open("reports/rank4-parity/e8-check.json"))
G8, g8 = e8["E8_gram"], e8["g"]
I8 = identity(8)
print("det", fdet(G8), "even", is_even(G8), " g^t G g = G:", multiply(multiply(transpose(g8), G8), g8) == G8,
      " g^2 = -I:", multiply(g8, g8) == [[-x for x in r] for r in I8])
cols8 = transpose([[x + y for x, y in zip(r1, r2)] for r1, r2 in zip(g8, I8)])
F8 = sorted({tuple(sum(k * c[i] for k, c in zip(coeffs, cols8)) % 2 for i in range(8)) for coeffs in product((0, 1), repeat=8)})
nz8 = [list(c) for c in F8 if any(c)]
print("fixed space classes:", len(F8), " nonzero:", len(nz8), " all singular (norm = 0 mod 4):", all(dot(c, matvec(G8, c)) % 4 == 0 for c in nz8),
      " pairwise orthogonal mod 2:", all(dot(x, matvec(G8, y)) % 2 == 0 for x in nz8 for y in nz8))
eps1 = 0
for a8 in nz8:
    for c8 in nz8:
        if a8 != c8:
            b8 = [x % 2 for x in matvec(G8, c8)]
            eps1 += symmetry_data(G8, a8, b8, g8)["epsilon"] == 1
print("ordered pairs of distinct nonzero classes:", len(nz8) * (len(nz8) - 1), " with epsilon(g) = 1:", eps1)
print("elapsed %.1f s" % (time.time() - T0))
print()

print("## 6. The three decomposable even odd-det classes with an element squaring to -I (det 9, 49, 121)")
for name, Gd in (("A2+A2", [[2, 1, 0, 0], [1, 2, 0, 0], [0, 0, 2, 1], [0, 0, 1, 2]]),
                 ("[[2,1],[1,4]]^2", [[2, 1, 0, 0], [1, 4, 0, 0], [0, 0, 2, 1], [0, 0, 1, 4]]),
                 ("[[2,1],[1,6]]^2", [[2, 1, 0, 0], [1, 6, 0, 0], [0, 0, 2, 1], [0, 0, 1, 6]])):
    grp = isometries(Gd, Gd, 5_000_000)
    sq = [T for T in grp if multiply(T, T) == MI4]
    resd = decide_all(Gd, grp)
    verd, zeros, maxfirst, wo = summarise(resd)
    print(f"{name}: det {determinant(Gd)} |Aut| {len(grp)} elements with g^2=-I: {len(sq)} verdicts {verd} witness orders {wo}")
print("elapsed %.1f s" % (time.time() - T0))
