"""Independent candidate scanner, full groups and cutoff replay for rank-four controls."""

import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import permutations, product
import json
from math import floor, gcd, prod
from pathlib import Path

from exact_theta import (EnumerationLimit, InvalidCertificate, dot, identity, integer,
                         inverse, ldl, matrix, matvec, multiply, require, symmetry_data,
                         transpose, vector)
from verify_lattice_census import box_shells, box_vectors

BITS = [list(v) for v in product((0,1),repeat=4)]


def det4(A):
    total = 0
    for perm in permutations(range(4)):
        inversions = sum(perm[i] > perm[j] for i in range(4) for j in range(i+1,4))
        total += (-1)**inversions * prod(A[i][perm[i]] for i in range(4))
    return total


def independent_candidates(H):
    """Bounded cross-entry boxes and diagonal scanning, without determinant solving."""
    integer(H,"H")
    require(1 <= H <= 24,"wrong feasibility domain")
    ratio = Fraction(3,4)
    upper = [Fraction(H)/ratio**(6-i) for i in range(4)]
    found = []

    def extend(G):
        k = len(G)
        mu, d = ldl(G)
        inv = inverse(G)
        prefix = prod(d)
        radii = [floor(sum(abs(mu[j][l])*d[l] for l in range(j+1))/2) for j in range(k)]
        for cross in product(*(range(-r,r+1) for r in radii)):
            coeff = matvec(inv,cross)
            projected = [sum(coeff[i]*mu[i][j] for i in range(j,k)) for j in range(k)]
            if any(abs(x)>Fraction(1,2) for x in projected):
                continue
            q = dot(cross,coeff)
            for diagonal in range(G[0][0],floor(upper[k]+sum(upper[:k])/4)+1):
                pivot = diagonal-q
                if pivot <= 0 or pivot < ratio*d[-1]:
                    continue
                r = 4-k
                if prefix*pivot**r*ratio**(r*(r-1)//2) > H:
                    break
                row = [list(old)+[x] for old,x in zip(G,cross)] + [list(cross)+[diagonal]]
                if k == 3:
                    delta = prefix*pivot
                    if delta.denominator == 1 and 1 <= delta <= H:
                        found.append(row)
                else:
                    extend(row)

    for a in range(1,H+1):
        if ratio**6*a**4 <= H:
            extend([[a]])
    return sorted(found)


def box_isometries(G,H,max_nodes=1_000_000):
    """Independent integer boxes and joining compatible column pairs."""
    matrix(G,"G",4)
    matrix(H,"H",4)
    integer(max_nodes,"max_nodes")
    require(max_nodes>0,"positive group limit required")
    ldl(G)
    ldl(H)
    vectors=box_vectors(H,max(G[i][i] for i in range(4)),max_nodes)
    images=[[v for v in vectors if dot(v,matvec(H,v))==G[i][i]] for i in range(4)]
    transformed={v:matvec(H,v) for v in vectors}
    nodes=0
    def tick():
        nonlocal nodes
        nodes+=1
        if nodes>max_nodes:
            raise EnumerationLimit("independent rank-four group join interrupted")
    pairs=[]
    for i,j in ((0,1),(2,3)):
        compatible=[]
        for u in images[i]:
            for v in images[j]:
                tick()
                if dot(u,transformed[v])==G[i][j]:
                    compatible.append((u,v))
        pairs.append(compatible)
    matches=[]
    for left in pairs[0]:
        for right in pairs[1]:
            tick()
            if all(dot(left[i],transformed[right[j]])==G[i][j+2] for i in range(2) for j in range(2)):
                T=transpose(left+right)
                if abs(det4(T))==1:
                    require(multiply(multiply(transpose(T),H),T)==G,"group Gram mismatch")
                    matches.append(T)
    return sorted(matches)


def modular_check(G,a,b,item):
    matrix(G,"G",4)
    ldl(G)
    vector(a,"a",4)
    vector(b,"b",4)
    require(all(v in (0,1) for v in a+b) and any(b) and dot(a,b)%2==0,"wrong cutoff hypotheses")
    K=matrix(item["index_two_basis"],"kernel basis",4)
    require(abs(det4(K))==2 and all(v%2==0 for v in matvec(transpose(K),b)),"wrong saturated kernel")
    G0=multiply(multiply(transpose(K),G),K)
    matrix(item["G0"],"stored G0",4)
    integer(item["d0"],"d0")
    require(item["G0"]==G0 and item["d0"]==det4(G0)==4*det4(G),"wrong kernel Gram/discriminant")
    N=1
    for row in inverse(G0):
        for v in row:
            N=N*v.denominator//gcd(N,v.denominator)
    M=16*N
    primes=[]
    remaining=M
    p=2
    while p*p<=remaining:
        if remaining%p==0:
            primes.append(p)
            while remaining%p==0:
                remaining//=p
        p+=1
    if remaining>1:
        primes.append(remaining)
    index=Fraction(M)
    factors=[]
    for p in primes:
        index*=Fraction(p+1,p)
        e,n=0,M
        while n%p==0:
            e+=1
            n//=p
        factors.append([p,e])
    require(index.denominator==1 and int(index)%6==0,"index integrality")
    for key,expected in (("N0",N),("unsigned_level",M),("index",int(index)),
                         ("cutoff",int(index)//6),("character_discriminant",4*det4(G0))):
        require(integer(item[key],key)==expected,"wrong cutoff metadata: "+key)
    require(item["level_factors"]==factors and item["weight"]=="2"
            and item["coefficient_index"]=="norm4"
            and item["theorem"]=="rank4-square-complex-identity-v1","wrong cutoff conventions")
    for p,e in item["level_factors"]:
        integer(p,"factor prime")
        integer(e,"factor exponent")
    return int(index)//6


def _cutoff_replay(doc,max_nodes):
    require(doc["schema"]=="work7-rank4-cutoff-v1","wrong cutoff schema")
    G,a,b=doc["G"],doc["a"],doc["b"]
    B=modular_check(G,a,b,doc["modular"])
    bound=integer(doc["run"]["norm4_bound"],"bound")
    builder_nodes=integer(doc["run"]["max_nodes"],"builder nodes")
    integer(max_nodes,"max_nodes")
    require(bound>=0 and builder_nodes>0 and max_nodes>0,"wrong resource limits")
    evidence=doc["evidence"]
    require(type(evidence["enumeration_complete"]) is bool,"wrong completeness flag")
    if not evidence["enumeration_complete"]:
        require(evidence["verdict"]=="unresolved" and "enumeration" not in evidence
                and isinstance(evidence["reason"],str),"interruption claimed partial proof")
        return {"verdict":"unresolved","certificate_replay_complete":True,"cutoff":B,
                "coefficient_replay_complete":False,"full_stabilizer_checked":False}
    enumeration=evidence["enumeration"]
    require(enumeration["complete"] is True and enumeration["norm4_bound"]==bound,"wrong enumeration bound")
    try:
        shells=box_shells(G,a,b,bound,max_nodes)
    except EnumerationLimit:
        return {"verdict":"unresolved","certificate_replay_complete":False,"cutoff":B,
                "coefficient_replay_complete":False,"full_stabilizer_checked":False}
    require(enumeration["shells"]==shells,"wrong complete coefficients")
    first=next((s for s in shells if s["signed_coefficient"]),None)
    verdict="proved_nonzero" if first else "proved_zero_by_modular_cutoff" if bound>=B else "unresolved"
    require(evidence["verdict"]==verdict,"wrong cutoff verdict")
    if first:
        require(evidence["nonzero_norm4"]==first["norm4"]
                and evidence["coefficient"]==first["signed_coefficient"],"wrong nonzero claim")
    return {"verdict":verdict,"certificate_replay_complete":True,"coefficient_replay_complete":True,
            "cutoff":B,"full_stabilizer_checked":False}


def cutoff_replay(doc,max_nodes=1_000_000):
    try:
        return _cutoff_replay(doc,max_nodes)
    except (KeyError,TypeError,IndexError,ValueError,ZeroDivisionError) as error:
        raise InvalidCertificate("malformed rank-four cutoff: "+str(error)) from error


def _controls_replay(doc,max_nodes):
    require(doc["schema"]=="work7-rank4-feasibility-controls-v1","wrong controls schema")
    require(doc["domain"]=="five named controls; not a census","wrong control scope")
    bound=integer(doc["run"]["norm4_bound"],"bound")
    builder_nodes=integer(doc["run"]["max_nodes"],"builder nodes")
    integer(max_nodes,"max_nodes")
    require(bound>=0 and builder_nodes>0 and max_nodes>0,"invalid control limits")
    expected=[("I4",identity(4)),("2I4",[[2*int(i==j) for j in range(4)] for i in range(4)]),
              ("A4",[[2 if i==j else -1 if abs(i-j)==1 else 0 for j in range(4)] for i in range(4)]),
              ("D4",[[2,-1,0,0],[-1,2,-1,-1],[0,-1,2,0],[0,-1,0,2]]),
              ("A3 plus line",[[2,-1,0,0],[-1,2,-1,0],[0,-1,2,0],[0,0,0,1]])]
    require([(r["name"],r["G"]) for r in doc["lattices"]]==expected,"wrong named-control coverage")
    reports=[]
    for row in doc["lattices"]:
        G=row["G"]
        integer(row["determinant"],"determinant")
        for T in row["automorphisms"]:
            matrix(T,"stored automorphism",4)
        require(row["determinant"]==det4(G),"wrong determinant")
        try:
            group=box_isometries(G,G,max_nodes)
        except EnumerationLimit:
            reports.append({"name":row["name"],"verdict":"unresolved","reason":"full group replay interrupted"})
            continue
        require(row["automorphisms"]==group,"full automorphism group mismatch")
        duals=[transpose(inverse(T)) for T in group]
        require([(r["a"],r["b"]) for r in row["pairs"]]==[(a,b) for a in BITS for b in BITS],"characteristic coverage mismatch")
        even=Counter()
        for p in row["pairs"]:
            a,b=p["a"],p["b"]
            vector(a,"characteristic a",4)
            vector(b,"characteristic b",4)
            indices,signs=[],[]
            for i,(T,dual) in enumerate(zip(group,duals)):
                if any((x-y)%2 for x,y in zip(matvec(T,a),a)):
                    continue
                transported=matvec(dual,b)
                if any((x-y)%2 for x,y in zip(transported,b)):
                    continue
                shift=[(x-y)/2 for x,y in zip(transported,b)]
                require(all(x.denominator==1 for x in shift),"nonintegral characteristic shift")
                indices.append(i)
                signs.append(int(dot(a,shift))%2)
            require(p["stabilizer"]=={"indices":indices,"epsilon":signs},"wrong full stabilizer/sign")
            if p["verdict"]=="proved_zero":
                require(p["evidence"]["kind"]=="symmetry" and 1 in signs
                        and symmetry_data(G,a,b,p["evidence"]["T"])["epsilon"]==1,"zero without symmetry proof")
            elif p["verdict"]=="proved_nonzero":
                require(1 not in signs and p["evidence"]["kind"]=="shell","inconsistent nonzero")
                enumeration=p["enumeration"]
                require(enumeration["complete"] is True and enumeration["norm4_bound"]==bound,"incomplete nonzero")
                shells=box_shells(G,a,b,bound,max_nodes)
                require(shells==enumeration["shells"],"nonzero shell histogram mismatch")
                claim=p["evidence"]
                require(claim["coefficient"]!=0 and any(s["norm4"]==claim["norm4"]
                        and s["signed_coefficient"]==claim["coefficient"] for s in shells),"wrong nonzero coefficient")
            else:
                raise InvalidCertificate("named control has unresolved/unknown stored verdict")
            if dot(a,b)%2==0:
                even[p["verdict"]]+=1
            else:
                require(p["verdict"]=="proved_zero","odd pair not zero")
        expected_summary={"group_order":len(group),"even_pairs":136,"odd_pairs":120,
                          "even_verdicts":dict(even),"all_pairs":256}
        require(row["summary"]==expected_summary,"wrong control summary")
        reports.append({"name":row["name"],"verdict":"verified",**expected_summary})
    return {"verdict":"verified" if all(r["verdict"]=="verified" for r in reports) else "unresolved",
            "controls":reports,"all_groups_and_stabilizers_replayed":all(r["verdict"]=="verified" for r in reports),
            "full_census_performed":False}


def controls_replay(doc,max_nodes=1_000_000):
    try:
        return _controls_replay(doc,max_nodes)
    except EnumerationLimit:
        return {"verdict":"unresolved","all_groups_and_stabilizers_replayed":False,"full_census_performed":False}
    except (KeyError,TypeError,IndexError,ValueError,ZeroDivisionError) as error:
        raise InvalidCertificate("malformed rank-four controls: "+str(error)) from error


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory",type=Path)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    candidate_raw=(args.directory/"candidates-det24.json").read_bytes()
    candidate_doc=json.loads(candidate_raw)
    require(candidate_doc["schema"]=="work7-rank4-candidate-bases-v1"
            and candidate_doc["determinant_bound"]==24
            and candidate_doc["generation_complete"] is True
            and candidate_doc["isometry_deduplication_performed"] is False,"wrong candidate pack scope")
    forms=candidate_doc["candidate_bases"]
    integer(candidate_doc["determinant_bound"],"candidate domain")
    for G in forms:
        matrix(G,"candidate Gram",4)
    require(forms==independent_candidates(24),"candidate coverage mismatch")
    determinants=[int(det4(G)) for G in forms]
    domain_checks=[]
    for H in range(1,25):
        subset=[G for G,d in zip(forms,determinants) if d<=H]
        require(subset==independent_candidates(H),"nested candidate coverage mismatch")
        domain_checks.append({"determinant_bound":H,"candidate_bases":len(subset)})
    sizing=json.loads((args.directory/"candidate-sizing.json").read_text(encoding="utf-8"))
    require(sizing["schema"]=="work7-rank4-candidate-sizing-v1"
            and sizing["isometry_deduplication_performed"] is False,"wrong sizing scope")
    require([r["determinant_bound"] for r in sizing["rows"]]==[1,2,4,8,12,16,24],"missing sizing domains")
    for row in sizing["rows"]:
        H=integer(row["determinant_bound"],"sizing domain")
        count=integer(row["candidate_bases"],"saved candidate count")
        require(1<=H<=24 and count==domain_checks[H-1]["candidate_bases"],"wrong saved size")
    raw=(args.directory/"controls.json").read_bytes()
    result=controls_replay(json.loads(raw))
    result["candidate_pack_sha256"]=hashlib.sha256(candidate_raw).hexdigest()
    result["candidate_domains_verified"]=domain_checks
    result["controls_sha256"]=hashlib.sha256(raw).hexdigest()
    raw=(args.directory/"cutoff-controls.json").read_bytes()
    cutoff_doc=json.loads(raw)
    require(cutoff_doc["schema"]=="work7-rank4-cutoff-controls-v1","wrong cutoff pack")
    result["cutoff_controls_sha256"]=hashlib.sha256(raw).hexdigest()
    result["cutoff_controls"]=[{"name":cert["name"],**cutoff_replay(cert)} for cert in cutoff_doc["certificates"]]
    if not all(r["certificate_replay_complete"] for r in result["cutoff_controls"]):
        result["verdict"]="unresolved"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result))
    return 0 if result["verdict"]=="verified" else 2


if __name__=="__main__":
    main()
