"""Reviewer's independent replication of the det-1729 level-4 claims of
RANK4-LEVEL4-SEARCH.md (review record: RANK4-LEVEL4-REVIEW.md).

Uses only the pre-existing accepted project modules (rank4_feasibility,
exact_theta) plus small local helpers.  Nothing from rank4_level4_search.py or
verify_rank4_level4.py is imported.  Run from the workspace root:

    python code/rank4_level4_review_check.py > reports/rank4-level4/review-independent-check.txt

Every printed line is an exact integer statement; the run takes about 15 s."""
import sys, time
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from rank4_feasibility import isometries, modular_data, determinant
from exact_theta import complete_shells, multiply, transpose, matvec, dot

def gram(G, cols):
    """Gram of the lattice spanned by the given (rational) column vectors."""
    n = len(cols)
    return [[sum(Fr(cols[i][k]) * G[k][l] * Fr(cols[j][l]) for k in range(4) for l in range(4))
             for j in range(n)] for i in range(n)]

def to_int(A):
    assert all(Fr(v).denominator == 1 for r in A for v in r), "non-integral Gram"
    return [[int(v) for v in r] for r in A]

def overlattice_basis(Kcols, v2):
    """Basis of L0 + Z v where v = v2/2 and v2 = sum c_i k_i with some c_i odd."""
    # solve 2v = v2 in the basis Kcols (integer coordinates)
    import itertools
    M = [[Fr(Kcols[j][i]) for j in range(4)] for i in range(4)]  # columns -> matrix
    # Gaussian elimination to solve M c = v2
    aug = [M[i] + [Fr(v2[i])] for i in range(4)]
    for i in range(4):
        p = next(r for r in range(i, 4) if aug[r][i] != 0)
        aug[i], aug[p] = aug[p], aug[i]
        aug[i] = [x / aug[i][i] for x in aug[i]]
        for r in range(4):
            if r != i and aug[r][i] != 0:
                f = aug[r][i]
                aug[r] = [x - f * y for x, y in zip(aug[r], aug[i])]
    c = [aug[i][4] for i in range(4)]
    assert all(ci.denominator == 1 for ci in c), "2v not in L0"
    c = [int(ci) for ci in c]
    i = next(i for i in range(4) if c[i] % 2)
    v = [Fr(x, 2) for x in v2]
    cols = [list(map(Fr, k)) for k in Kcols]
    cols[i] = v
    return cols, c

G = [[3,1,1,0],[1,9,-2,4],[1,-2,7,-2],[0,4,-2,13]]
a = [1,0,0,1]; b = [1,1,1,1]
print("det G =", determinant(G), " a.b =", dot(a,b))

# 1. signed coefficients of Theta_G[a,b] through the full modular cutoff
md = modular_data(G, a, b)
print("modular data: d0 =", md["d0"], "N0 =", md["N0"], "M =", md["unsigned_level"],
      "mu =", md["index"], "B =", md["cutoff"])
t0 = time.time()
sh = complete_shells(G, a, b, md["cutoff"], max_nodes=20_000_000)
nz = [s for s in sh["shells"] if s["signed_coefficient"] != 0]
print("cutoff enumeration: shells =", len(sh["shells"]), "vectors =", sum(s["vector_count"] for s in sh["shells"]),
      "nodes =", sh["visited_nodes"], "nonzero signed shells =", len(nz), f"({time.time()-t0:.1f}s)")

# 2. Aut(L)
A = isometries(G, G, max_nodes=5_000_000)
print("Aut(L) size =", len(A), " elements:", A)

# 3. sign kernel and its two overlattices
K = [[2,0,0,0],[-1,1,0,0],[-1,0,1,0],[-1,0,0,1]]  # columns 2e1, e2-e1, e3-e1, e4-e1
G0 = to_int(gram(G, K))
print("G0 =", G0, " det G0 =", determinant(G0))
A0 = isometries(G0, G0, max_nodes=5_000_000)
print("Aut(L0) size =", len(A0))
O1cols, c1 = overlattice_basis(K, a)                       # L0 + Z a/2
O2cols, c2 = overlattice_basis(K, [a[0]+2, a[1], a[2], a[3]])  # L0 + Z (a+2e1)/2
O1 = to_int(gram(G, O1cols)); O2 = to_int(gram(G, O2cols))
print("Gram O1 =", O1, "det", determinant(O1))
print("Gram O2 =", O2, "det", determinant(O2))

# 4. Conway-Sloane L^+/-(1,7,13,19) from my own transcription
def cs_gram(p, q, r, s_, sigma):
    cols = [[3*sigma,-1,-1,-1],[1,3*sigma,1,-1],[1,-1,3*sigma,1],[1,1,-1,3*sigma]]
    D = [Fr(p,12),Fr(q,12),Fr(r,12),Fr(s_,12)]
    return to_int([[sum(D[k]*cols[i][k]*cols[j][k] for k in range(4)) for j in range(4)] for i in range(4)])
Lp = cs_gram(1,7,13,19,+1); Lm = cs_gram(1,7,13,19,-1)
print("L+ =", Lp, "det", determinant(Lp)); print("L- =", Lm, "det", determinant(Lm))
print("isometries L+ vs L- :", len(isometries(Lp, Lm, max_nodes=5_000_000)))
print("isometries O1 vs L- :", len(isometries(O1, Lm, max_nodes=5_000_000)),
      " O1 vs L+ :", len(isometries(O1, Lp, max_nodes=5_000_000)))
print("isometries O2 vs L+ :", len(isometries(O2, Lp, max_nodes=5_000_000)),
      " O2 vs L- :", len(isometries(O2, Lm, max_nodes=5_000_000)))
print("isometries O1 vs O2 :", len(isometries(O1, O2, max_nodes=5_000_000)))

# 5. the new fact: the b=(1,1,1,1) kernels of L+ and L- are isometric
Kp = to_int(gram(Lp, K)); Km = to_int(gram(Lm, K))
print("det ker+ =", determinant(Kp), " det ker- =", determinant(Km))
U = isometries(Kp, Km, max_nodes=5_000_000)
print("isometries ker(L+) vs ker(L-):", len(U), U[:1])
print("isometries ker(L+) vs G0:", len(isometries(Kp, G0, max_nodes=5_000_000)))

# 6. side claims: Earnest-Nipp g1 vs L^-(7,13,19,49); det-3844 pair
def gram2(c):
    a11,a22,a33,a44,a12,a13,a23,a14,a24,a34 = c
    return [[2*a11,a12,a13,a14],[a12,2*a22,a23,a24],[a13,a23,2*a33,a34],[a14,a24,a34,2*a44]]
g1 = gram2((6,8,10,15,-4,0,5,1,4,8)); g2 = gram2((6,8,10,15,-1,6,1,3,-7,2))
print("EN g1 vs CS L-(7,13,19,49):", len(isometries(g1, cs_gram(7,13,19,49,-1), max_nodes=5_000_000)),
      " g2 vs L+(7,13,19,49):", len(isometries(g2, cs_gram(7,13,19,49,+1), max_nodes=5_000_000)))
S1 = [[4,1,0,0],[1,8,0,0],[0,0,8,2],[0,0,2,16]]
S2 = [[4,-2,1,1],[-2,8,-3,-1],[1,-3,10,-4],[1,-1,-4,18]]
print("det S1, S2 =", determinant(S1), determinant(S2), " isometries S1 vs S2:", len(isometries(S1, S2, max_nodes=5_000_000)))
z = [0,0,0,0]
h1 = {s["norm4"]: s["vector_count"] for s in complete_shells(S1, z, z, 4*400, max_nodes=20_000_000)["shells"]}
h2 = {s["norm4"]: s["vector_count"] for s in complete_shells(S2, z, z, 4*400, max_nodes=20_000_000)["shells"]}
print("S1/S2 theta agree through norm 400:", h1 == h2, " first terms:", [h1.get(4*n,0) for n in range(0,12)])
# are the b=(1,1,1,1) kernels of the det-961 pair really S1, S2 up to isometry?
P = [[17,-5,5,2],[-5,13,6,5],[5,6,9,5],[2,5,5,5]]; Q = [[17,-2,1,7],[-2,13,9,1],[1,9,9,2],[7,1,2,5]]
KP = to_int(gram(P, K)); KQ = to_int(gram(Q, K))
print("ker(P) vs S1:", len(isometries(KP, S1, max_nodes=5_000_000)), " ker(Q) vs S2:", len(isometries(KQ, S2, max_nodes=5_000_000)),
      " ker(P) vs S2:", len(isometries(KP, S2, max_nodes=5_000_000)), " ker(Q) vs S1:", len(isometries(KQ, S1, max_nodes=5_000_000)))
