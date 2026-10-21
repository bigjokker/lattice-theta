"""Control for the Eisenstein (order-three) family used by item 3.

Builds the sign kernel L0 of the determinant-23,716 counterexample from the
Hermitian parameters (p, q, m, n) = (8, 14, -4, 2) with the trace-form formula
of Fable-rank4-eisenstein-search.md section 3, checks the order-three
automorphism R, the determinant formula, the three overlattices of the line
Z[w]e1, their isometry with the manuscript's G, L1, L2, and the levels of all
105 configurations over L0 (all five lines level 2).  Also counts the reduced
(p0, delta0, q0) triples of the enumeration universe.  Exact arithmetic only;
uses rank4_feasibility, rank4_level4_search and exact_theta.  Run from the
workspace root; takes a few seconds."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from rank4_feasibility import isometries, determinant
from rank4_level4_search import analyse_kernel, hnf_basis, gram_of
from exact_theta import multiply, transpose, matvec, dot

def trace_gram(p, q, m, n):
    """Trace form of H' = [[p, delta/sqrt(-3)], [conj, q]], delta = m + n*omega,
    on the Z-basis (e1, w e1, e2, w e2); entry = 2 Re H'(b_i, b_j)."""
    A = [[2*p, -p], [-p, 2*p]]
    Dq = [[2*q, -q], [-q, 2*q]]
    C = [[n, -m], [m - n, n]]          # rows: e1, w e1 ; cols: e2, w e2
    G = [[0]*4 for _ in range(4)]
    for i in range(2):
        for j in range(2):
            G[i][j] = A[i][j]; G[2+i][2+j] = Dq[i][j]; G[i][2+j] = C[i][j]; G[2+j][i] = C[i][j]
    return G

R = [[0,-1,0,0],[1,-1,0,0],[0,0,0,-1],[0,0,1,-1]]
# det-23716 example: G0 = 2M with M = Tr(H), H = [[4, alpha],[conj, 7]], delta_M = -2 + w
p, q, m, n = 8, 14, -4, 2
G0 = trace_gram(p, q, m, n)
Ndelta = m*m - m*n + n*n
print("G0 =", G0)
print("R^t G0 R == G0:", multiply(multiply(transpose(R), G0), R) == G0)
print("det G0 =", determinant(G0), " (3pq - N(delta))^2 =", (3*p*q - Ndelta)**2, " D' =", (3*p*q - Ndelta)//2)
# the three overlattices of the line Z[w]e1 mod 2
def over(G0, v):
    cols = [[2*int(i==j) for i in range(4)] for j in range(4)]
    B = hnf_basis(cols + [v])
    g = gram_of(G0, B)
    assert all(x % 4 == 0 for r in g for x in r), "not integral"
    return [[x//4 for x in r] for r in g]
O1, O2, O3 = over(G0, [1,0,0,0]), over(G0, [0,1,0,0]), over(G0, [1,1,0,0])
G = [[12,12,-4,-1],[12,16,-6,2],[-4,-6,28,-14],[-1,2,-14,28]]
L1 = [[48,12,-8,-2],[12,4,-3,1],[-8,-3,28,-14],[-2,1,-14,28]]
L2 = [[48,36,-8,-2],[36,28,-7,0],[-8,-7,28,-14],[-2,0,-14,28]]
print("det O1,O2,O3 =", determinant(O1), determinant(O2), determinant(O3))
print("O3 ~ G:", len(isometries(O3, G, 5_000_000)), " O1 ~ G:", len(isometries(O1, G, 5_000_000)), " O1 ~ L1:", len(isometries(O1, L1, 5_000_000)), " O2 ~ L2:", len(isometries(O2, L2, 5_000_000)))
print("Aut(O3) =", len(isometries(O3, O3, 5_000_000)))
res = analyse_kernel(G0, prefix=60, max_nodes=5_000_000)
print("Aut(L0) order", res["aut_order"], "tally", res["tally"])
for pr in res["zero_pairs"]:
    print("  zero pair", pr["v1"], pr["v2"], "-> v3", pr["v3"], "level", pr["level"], "scaling", pr["smallest_integral_scaling"], "det L", pr["det_scaled_L"])
# size of the enumeration universe for 5 <= D' <= 41 and 5 <= D' <= 154
from math import isqrt
def count(Dmax):
    forms = members = 0
    for Dp in range(5, Dmax+1):
        E = 2*Dp
        for p0 in range(1, isqrt(Dp)+1):
            for mm in range(-2*p0, 2*p0+1):
                for nn in range(-2*p0, 2*p0+1):
                    N = mm*mm - mm*nn + nn*nn
                    if N > p0*p0: continue
                    if (E + N) % (3*p0): continue
                    q0 = (E + N)//(3*p0)
                    if q0 < p0: continue
                    forms += 1
    return forms
print("reduced (p0, delta0, q0) triples with 5<=D'<=41:", count(41), "; with 5<=D'<=154:", count(154))
