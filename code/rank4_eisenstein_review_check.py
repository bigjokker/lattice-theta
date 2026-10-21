"""Reviewer's independent replication of the det-1156 level-2 claims of
RANK4-EISENSTEIN-SEARCH.md (review record: RANK4-EISENSTEIN-REVIEW.md).
Uses only rank4_feasibility, exact_theta and rank4_level4_search.analyse_kernel;
nothing from rank4_eisenstein_search.py or verify_rank4_eisenstein.py.
Run from the workspace root:
    python code/rank4_eisenstein_review_check.py > reports/rank4-eisenstein/review-independent-check.txt
"""
import sys, os, time
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as Fr
from rank4_feasibility import isometries, modular_data, determinant
from exact_theta import complete_shells, multiply, transpose, matvec, dot, inverse
from rank4_level4_search import analyse_kernel, hnf_basis, gram_of
G = [[4,-2,1,-2],[-2,8,2,-2],[1,2,6,-4],[-2,-2,-4,12]]; a=[0,1,0,1]; b=[1,0,1,0]
print("det G =", determinant(G), " a.b =", dot(a,b), " gcd entries =", __import__('math').gcd(*[abs(x) for r in G for x in r]), " even:", all(G[i][i]%2==0 for i in range(4)))
md = modular_data(G, a, b)
print("modular: d0", md["d0"], "N0", md["N0"], "M", md["unsigned_level"], "mu", md["index"], "B", md["cutoff"])
sh = complete_shells(G, a, b, md["cutoff"], max_nodes=20_000_000)
print("cutoff shells", len(sh["shells"]), "vectors", sum(s["vector_count"] for s in sh["shells"]), "nonzero signed", sum(1 for s in sh["shells"] if s["signed_coefficient"]))
sh2 = complete_shells(G, a, b, 2000, max_nodes=50_000_000)
print("through norm4 2000: vectors", sum(s["vector_count"] for s in sh2["shells"]), "nonzero signed", sum(1 for s in sh2["shells"] if s["signed_coefficient"]))
A = isometries(G, G, 5_000_000); print("Aut(L) =", len(A))
K = [[2,0,-1,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
G0 = multiply(multiply(transpose(K), G), K); print("G0 =", G0, "det", determinant(G0))
A0 = isometries(G0, G0, 5_000_000); print("Aut(L0) =", len(A0))
Kinv = inverse(K)
v1 = [int(x)%2 for x in matvec(Kinv, a)]; v2 = [int(x)%2 for x in matvec(Kinv, [a[0]+2,a[1],a[2],a[3]])]; v3 = [(x+y)%2 for x,y in zip(v1,v2)]
print("v1 v2 v3 =", v1, v2, v3)
mappers = [T for T in A0 if [x%2 for x in matvec(T,v1)] == v2]
swappers = [T for T in A0 if [x%2 for x in matvec(T,v1)] == v2 and [x%2 for x in matvec(T,v2)] == v1 and [x%2 for x in matvec(T,v3)] == v3]
def order(T):
    P, k = T, 1
    I = [[int(i==j) for j in range(4)] for i in range(4)]
    while P != I: P = multiply(P, T); k += 1
    return k
print("elements mapping v1->v2:", len(mappers), "orders", [order(T) for T in mappers], " swappers fixing v3:", len(swappers))
def over(G0, v):
    cols = [[2*int(i==j) for i in range(4)] for j in range(4)]
    B = hnf_basis(cols + [v]); g = gram_of(G0, B)
    assert all(x % 4 == 0 for r in g for x in r); return [[x//4 for x in r] for r in g]
O1, O2, O3 = over(G0, v1), over(G0, v2), over(G0, v3)
print("det O1 O2 O3:", determinant(O1), determinant(O2), determinant(O3), " O1~O2:", len(isometries(O1,O2,5_000_000)), " O3~G:", len(isometries(O3,G,5_000_000)))
# Hermitian parametrization (a,q,beta) = (4,6,2+5w): p=8, q=6, delta = 2beta = 4+10w
def trace_gram(p,q,m,n):
    A=[[2*p,-p],[-p,2*p]]; Dq=[[2*q,-q],[-q,2*q]]; C=[[n,-m],[m-n,n]]; Gm=[[0]*4 for _ in range(4)]
    for i in range(2):
        for j in range(2):
            Gm[i][j]=A[i][j]; Gm[2+i][2+j]=Dq[i][j]; Gm[i][2+j]=C[i][j]; Gm[2+j][i]=C[i][j]
    return Gm
p,q,m,n = 8,6,4,10; Gh = trace_gram(p,q,m,n); R=[[0,-1,0,0],[1,-1,0,0],[0,0,0,-1],[0,0,1,-1]]
print("trace G0 =", Gh, "det", determinant(Gh), "(3pq-N)^2 =", (3*p*q-(m*m-m*n+n*n))**2, " R-inv:", multiply(multiply(transpose(R),Gh),R)==Gh, " D' =", (3*p*q-(m*m-m*n+n*n))//2)
print("trace G0 ~ kernel G0:", len(isometries(Gh, G0, 5_000_000)))
Oh = over(Gh, [1,1,0,0]); print("line O_v3 of trace form ~ G:", len(isometries(Oh, G, 5_000_000)))
res = analyse_kernel(Gh, prefix=400, max_nodes=5_000_000)
print("analyse_kernel(trace G0, 400): Aut order", res["aut_order"], "tally", res["tally"])
for pr in res["zero_pairs"]: print("   pair", pr["v1"], pr["v2"], "level", pr["level"], "det L", pr["det_scaled_L"], "scaling", pr["smallest_integral_scaling"])
