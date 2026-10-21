"""Reviewer's independent re-enumeration of the Eisenstein family F for
5 <= D' <= 41 (review record: RANK4-EISENSTEIN-REVIEW.md): reduced triples,
five lines, integrality, line level by the complete Aut(L0), primitive
determinant of L, isometry of every level-2 member with the det-1156 lattice.
Uses only rank4_feasibility, exact_theta and rank4_level4_search.hnf_basis/gram_of.
Run from the workspace root (about 10 s):
    python code/rank4_eisenstein_review_enum.py > reports/rank4-eisenstein/review-independent-enum.txt
"""
import sys, math, time
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from math import isqrt, gcd
from collections import Counter
from rank4_feasibility import isometries, determinant
from exact_theta import multiply, transpose, matvec
from rank4_level4_search import hnf_basis, gram_of
def trace_gram(p,q,m,n):
    A=[[2*p,-p],[-p,2*p]]; Dq=[[2*q,-q],[-q,2*q]]; C=[[n,-m],[m-n,n]]; Gm=[[0]*4 for _ in range(4)]
    for i in range(2):
        for j in range(2):
            Gm[i][j]=A[i][j]; Gm[2+i][2+j]=Dq[i][j]; Gm[i][2+j]=C[i][j]; Gm[2+j][i]=C[i][j]
    return Gm
R=[[0,-1,0,0],[1,-1,0,0],[0,0,0,-1],[0,0,1,-1]]
LIFTS=[(1,0,0,0),(0,0,1,0),(1,0,1,0),(1,0,0,1),(1,0,1,1)]   # (1,0),(0,1),(1,1),(1,w),(1,w^2)
def over(G0, v):
    cols=[[2*int(i==j) for i in range(4)] for j in range(4)]
    g=gram_of(G0, hnf_basis(cols+[list(v)])); return g
G=[[4,-2,1,-2],[-2,8,2,-2],[1,2,6,-4],[-2,-2,-4,12]]
t0=time.time(); members=Counter(); level2=[]; triples=0; lines=0; nonint=0
for Dp in range(5,42):
    E=2*Dp
    for p0 in range(1,isqrt(Dp)+1):
        for m0 in range(-2*p0,2*p0+1):
            for n0 in range(-2*p0,2*p0+1):
                N=m0*m0-m0*n0+n0*n0
                if N>p0*p0 or (E+N)%(3*p0): continue
                q0=(E+N)//(3*p0)
                if q0<p0: continue
                triples+=1
                G0=trace_gram(p0,q0,m0,n0)
                assert determinant(G0)==E*E and multiply(multiply(transpose(R),G0),R)==G0
                for v in LIFTS:
                    lines+=1
                    v1=list(v); v2=[x%2 for x in matvec(R,v1)]; v3=[(x+y)%2 for x,y in zip(v1,v2)]
                    ok=True
                    for u in (v1,v2,v3):
                        gu=matvec(G0,u)
                        if any(x%2 for x in gu) or sum(a*b for a,b in zip(u,gu))%4: ok=False
                    if not ok: nonint+=1; continue
                    members[Dp]+=1
                    A0=isometries(G0,G0,5_000_000)
                    swap=any([x%2 for x in matvec(T,v1)]==v2 and [x%2 for x in matvec(T,v2)]==v1 for T in A0)
                    if not swap:
                        g=over(G0,v3); assert all(x%4==0 for r in g for x in r); GL=[[x//4 for x in r] for r in g]
                        c=gcd(*[abs(x) for r in GL for x in r]); GLp=[[x//c for x in r] for r in GL]
                        level2.append((Dp,p0,q0,m0,n0,v,determinant(GL),c,determinant(GLp),len(isometries(GLp,G,5_000_000)),len(A0)))
print("triples",triples,"line candidates",lines,"non-integral lines",nonint,"members",sum(members.values()),f"({time.time()-t0:.0f}s)")
print("members per D':",dict(sorted(members.items())))
print("level-2 members (D',p0,q0,m0,n0,line,detL,c,det primitive,#iso to det-1156 G,|Aut L0|):")
for x in level2: print("  ",x)
