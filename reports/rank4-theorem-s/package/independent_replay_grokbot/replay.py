"""Independent exact replay of the 502,740 face/marking cases of Theorem S (Grok Bot, 2026-10-23).

Independent of rank4_six_cells.py / rank4_secondary_fan.py / rank4_six_minima.py (none imported).
Inputs: fan.json (rays, as data) and, only for the final comparison, cells.json.
All deciding arithmetic is exact (int / fractions.Fraction).  No floating point anywhere.

Usage: python3 replay.py FAN_JSON CELLS_JSON OUTDIR [--faces a:b] [--seed 1]
"""
import gzip, json, sys, time, random, hashlib, itertools
from fractions import Fraction as F

# ---------------- forms ----------------
IDX = [(0, 0), (1, 1), (2, 2), (3, 3), (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]

def toQ(q):
    Q = [[F(0)] * 4 for _ in range(4)]
    for t, (i, j) in enumerate(IDX):
        Q[i][j] = Q[j][i] = F(q[t])
    return Q

def tovec(Q):
    return [Q[i][j] for (i, j) in IDX]

def qval(Q, v):
    return sum(Q[i][j] * v[i] * v[j] for i in range(4) for j in range(4))

def rayval(r, v):  # R(v) for 10-vector r
    return (r[0]*v[0]*v[0] + r[1]*v[1]*v[1] + r[2]*v[2]*v[2] + r[3]*v[3]*v[3]
            + 2*(r[4]*v[0]*v[1] + r[5]*v[0]*v[2] + r[6]*v[0]*v[3] + r[7]*v[1]*v[2] + r[8]*v[1]*v[3] + r[9]*v[2]*v[3]))

def ldl(Q):
    """Exact Q = sum_i d_i (v_i + sum_{j>i} l[i][j] v_j)^2 ; returns (d, l) or None if not PD."""
    A = [row[:] for row in Q]
    d = [None]*4; l = [[F(0)]*4 for _ in range(4)]
    for i in range(4):
        if A[i][i] <= 0:
            return None
        d[i] = A[i][i]
        for j in range(i+1, 4):
            l[i][j] = A[i][j] / d[i]
        for j in range(i+1, 4):
            for k in range(i+1, 4):
                A[j][k] -= d[i] * l[i][j] * l[i][k]
    return d, l

def is_pd(Q):
    return ldl(Q) is not None

def isqrt_floor_frac(x):  # floor(sqrt(x)) for Fraction x >= 0, exact
    p, q = x.numerator, x.denominator
    from math import isqrt
    return isqrt(p * q) // q

def short_vectors(Q, B, dl=None):
    """All v in Z^4 with Q(v) <= B (exact, certified integer bounds)."""
    d, l = dl or ldl(Q)
    out = []
    v = [0]*4
    def rec(i, rem):
        c = sum(l[i][j] * v[j] for j in range(i+1, 4))
        s = isqrt_floor_frac(rem / d[i])
        lo = -c - s - 1; hi = -c + s + 1
        import math
        for x in range(math.floor(lo), math.ceil(hi) + 1):
            t = d[i] * (x + c) ** 2
            if t <= rem:
                v[i] = x
                if i == 0:
                    out.append(tuple(v))
                else:
                    rec(i-1, rem - t)
        v[i] = 0
    rec(3, F(B))
    return out

def label(v):
    return sum((int(x) % 2) << i for i, x in enumerate(v))

def bits(k):
    return [(k >> i) & 1 for i in range(4)]

SHIFTS = [w for w in itertools.product((-1, 0, 1), repeat=4) if sum(map(abs, w)) <= 2]

def minimal_sets(Q):
    """Complete minimal-vector set of every class k=1..15 (exact)."""
    dl = ldl(Q)
    ub = {}
    for k in range(1, 16):
        b = bits(k)
        ub[k] = min(qval(Q, [b[i] + 2*w[i] for i in range(4)]) for w in SHIFTS)
    vs = short_vectors(Q, max(ub.values()), dl)
    best = {k: None for k in range(1, 16)}; sets = {k: [] for k in range(1, 16)}
    for v in vs:
        k = label(v)
        if k == 0:
            continue
        n = qval(Q, v)
        if best[k] is None or n < best[k]:
            best[k] = n; sets[k] = [v]
        elif n == best[k]:
            sets[k].append(v)
    for k in range(1, 16):
        assert best[k] is not None and best[k] <= ub[k]
        sets[k].sort()
    return sets, best

# ---------------- markings ----------------
def markings():
    out = set()
    for m in range(1 << 16):
        N = [[(m >> (4*i+j)) & 1 for j in range(4)] for i in range(4)]
        if any(sum(N[i][k]*N[k][j] for k in range(4)) % 2 for i in range(4) for j in range(4)):
            continue
        img = {label([sum(N[i][j]*bits(k)[j] for j in range(4)) for i in range(4)]) for k in range(16)}
        if len(img) != 4:
            continue
        out.add(tuple(k ^ label([sum(N[i][j]*bits(k)[j] for j in range(4)) for i in range(4)]) for k in range(16)))
    return sorted(out)

def pairs_of(tau):
    return sorted({tuple(sorted((k, tau[k]))) for k in range(16) if tau[k] != k})

# ---------------- exact simplex ----------------
def simplex(A, b, c):
    """min c.x  s.t. A x = b, x >= 0 (exact, Bland). Returns x or None if infeasible. Assumes bounded."""
    m, n = len(A), len(A[0]) if A else 0
    A = [[F(x) for x in row] for row in A]; b = [F(x) for x in b]
    for i in range(m):
        if b[i] < 0:
            A[i] = [-x for x in A[i]]; b[i] = -b[i]
    # tableau with artificials n..n+m-1
    T = [A[i] + [F(int(i == r)) for r in range(m)] + [b[i]] for i in range(m)]
    basis = list(range(n, n+m))
    def run(cost, allowed):
        while True:
            # reduced costs
            z = [cost[j] - sum(cost[basis[i]] * T[i][j] for i in range(len(T))) for j in range(n+m)]
            ent = next((j for j in range(n+m) if allowed[j] and z[j] < 0), None)
            if ent is None:
                return True
            rows = [(T[i][-1] / T[i][ent], basis[i], i) for i in range(len(T)) if T[i][ent] > 0]
            if not rows:
                return False  # unbounded
            _, _, r = min(rows)
            pivot(r, ent)
    def pivot(r, e):
        p = T[r][e]
        T[r] = [x / p for x in T[r]]
        for i in range(len(T)):
            if i != r and T[i][e] != 0:
                f = T[i][e]
                T[i] = [a - f*bb for a, bb in zip(T[i], T[r])]
        basis[r] = e
    c1 = [F(0)]*n + [F(1)]*m
    run(c1, [True]*(n+m))
    if sum(T[i][-1] for i in range(len(T)) if basis[i] >= n) != 0:
        return None
    # drive artificials out
    i = 0
    while i < len(T):
        if basis[i] >= n:
            j = next((j for j in range(n) if T[i][j] != 0), None)
            if j is None:
                del T[i]; del basis[i]; continue
            pivot(i, j)
        i += 1
    c2 = [F(x) for x in c] + [F(0)]*m
    if not run(c2, [True]*n + [False]*m):
        raise RuntimeError("unbounded")
    x = [F(0)]*n
    for i, bi in enumerate(basis):
        x[bi] = T[i][-1]
    return x

def max_support(E, n):
    """J = max support of {mu>=0, E mu = 0}; returns (J, mu) with mu>=0, E mu=0, mu_j>=1 on J."""
    m = len(E)
    # vars: mu(n), s(n), u(n), w(n)
    A = []; b = []
    for p in range(m):
        A.append(list(E[p]) + [0]*(3*n)); b.append(0)
    for j in range(n):
        row = [0]*(4*n); row[j] = 1; row[n+j] = -1; row[2*n+j] = -1; A.append(row); b.append(0)
        row = [0]*(4*n); row[n+j] = 1; row[3*n+j] = 1; A.append(row); b.append(1)
    c = [0]*n + [-1]*n + [0]*(2*n)
    x = simplex(A, b, c)
    s = x[n:2*n]
    J = [j for j in range(n) if s[j] == 1]
    assert all(s[j] in (0, 1) for j in range(n)), s
    return J, x[:n]

def stiemke(E, n, J):
    """y with E^t y >= 0 and (E^t y)_j >= 1 for j not in J."""
    m = len(E)
    # y = yp - ym ; (E^t y)_j - z_j = l_j
    A = []; b = []
    for j in range(n):
        A.append([E[p][j] for p in range(m)] + [-E[p][j] for p in range(m)] + [-int(jj == j) for jj in range(n)])
        b.append(0 if j in J else 1)
    x = simplex(A, b, [0]*(2*m+n))
    if x is None:
        return None
    return [x[p] - x[m+p] for p in range(m)]

# ---------------- linear algebra ----------------
def rref(M):
    M = [[F(x) for x in row] for row in M]; piv = []; r = 0
    cols = len(M[0]) if M else 0
    for c in range(cols):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]
        pv = M[r][c]; M[r] = [x/pv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f*bb for a, bb in zip(M[i], M[r])]
        piv.append(c); r += 1
    return M[:r], piv

def nullspace(E, n):
    if not E:
        return [[F(int(i == j)) for i in range(n)] for j in range(n)]
    R, piv = rref(E)
    free = [j for j in range(n) if j not in piv]
    out = []
    for f in free:
        v = [F(0)]*n; v[f] = F(1)
        for i, pc in enumerate(piv):
            v[pc] = -R[i][f]
        out.append(v)
    return out

def solve_left(D, E):
    """Z with D = Z E (rows of D in row space of E) or None. D: r x n, E: m x n."""
    m = len(E); n = len(D[0])
    Z = []
    for drow in D:
        if not any(drow):
            Z.append([F(0)]*m); continue
        if m == 0:
            return None
        # solve E^t z = drow
        Aug = [[E[p][j] for p in range(m)] + [drow[j]] for j in range(n)]
        R, piv = rref(Aug)
        if m in piv:
            return None
        z = [F(0)]*m
        for i, pc in enumerate(piv):
            z[pc] = R[i][m]
        Z.append(z)
    return Z

def det_int(T):
    M = [[F(x) for x in r] for r in T]; d = F(1)
    for c in range(4):
        p = next((i for i in range(c, 4) if M[i][c] != 0), None)
        if p is None: return 0
        if p != c: M[c], M[p] = M[p], M[c]; d = -d
        d *= M[c][c]
        for i in range(c+1, 4):
            f = M[i][c] / M[c][c]; M[i] = [a - f*bb for a, bb in zip(M[i], M[c])]
    return d

def mod2_perm(T):
    return tuple(label([sum(T[i][j]*bits(k)[j] for j in range(4)) for i in range(4)]) for k in range(16))

def congr(T, Q):
    return [[sum(T[a][i]*Q[a][b]*T[b][j] for a in range(4) for b in range(4)) for j in range(4)] for i in range(4)]

def isometries(Q, tau):
    """All integral T with T^t Q T = Q and column i of T in class tau(2^i) (exact, complete)."""
    diag = [Q[i][i] for i in range(4)]
    vs = short_vectors(Q, max(diag))
    cand = [[v for v in vs if qval(Q, v) == diag[i] and label(v) == tau[1 << i]] for i in range(4)]
    def ip(u, v):
        return sum(Q[a][b]*u[a]*v[b] for a in range(4) for b in range(4))
    cols = []
    def rec(i):
        if i == 4:
            T = [[cols[j][r] for j in range(4)] for r in range(4)]
            if det_int(T) in (1, -1):
                yield T
            return
        for v in cand[i]:
            if all(ip(v, cols[j]) == Q[j][i] for j in range(i)):
                cols.append(v)
                yield from rec(i+1)
                cols.pop()
    yield from rec(0)

# ---------------- main ----------------
def fs(x):
    return str(x)

def main():
    fan_path, cells_path, outdir = sys.argv[1:4]
    seed = 1
    t0 = time.time(); c0 = time.process_time()
    fan = json.load(open(fan_path))
    cones = [c["rays"] for c in fan["cones"]]
    taus = markings()
    assert len(taus) == 210, len(taus)
    claude = {}
    for c in json.load(open(cells_path))["cells"]:
        claude[(c["cone"], tuple(c["subset"]), tuple(c["tau"]))] = c
    rng = random.Random(seed)
    counts = {"multiplicity_mismatch": 0, "infeasible": 0, "lower_face": 0, "certified": 0, "uncertified": 0}
    faces = 0; pd_faces = 0; consistency_failures = []; uncert = []; mine_cert = set()
    savedT_ok = 0; savedT_tried = 0; own_T = 0
    out = gzip.open(f"{outdir}/dispositions.jsonl.gz", "wt")
    for ci, rays in enumerate(cones):
        for mask in range(1, 1 << 10):
            S = [j for j in range(10) if mask >> j & 1]; faces += 1
            R = [rays[j] for j in S]; n = len(S)
            Q1 = toQ([sum(r[t] for r in R) for t in range(10)])
            if not is_pd(Q1):
                continue
            pd_faces += 1
            sets, mins = minimal_sets(Q1)
            w2 = [rng.randint(1, 9) for _ in S]
            Q2 = toQ([sum(w*r[t] for w, r in zip(w2, R)) for t in range(10)])
            sets2, _ = minimal_sets(Q2)
            if sets2 != sets:
                consistency_failures.append({"cone": ci, "subset": S, "w2": w2})
            # per-class functional: R_j(v) must agree for all v in M_k
            fun = {}
            for k in range(1, 16):
                fk = {tuple(rayval(r, v) for r in R) for v in sets[k]}
                if len(fk) != 1:
                    consistency_failures.append({"cone": ci, "subset": S, "class": k, "why": "functional differs"})
                fun[k] = next(iter(fk))
            for tau in taus:
                P = pairs_of(tau)
                rec = {"cone": ci, "subset": S, "tau": list(tau)}
                bad = next(((a, b) for a, b in P if len(sets[a]) != len(sets[b])), None)
                if bad:
                    a, b = bad
                    rec.update(d="multiplicity_mismatch", pair=[a, b], Ma=sets[a], Mb=sets[b])
                    counts[rec["d"]] += 1; out.write(json.dumps(rec) + "\n"); continue
                E = [[fun[a][j] - fun[b][j] for j in range(n)] for a, b in P]
                J, mu = max_support(E, n)
                # exact verification of mu
                assert all(x >= 0 for x in mu) and all(mu[j] >= 1 for j in J) and all(mu[j] == 0 for j in range(n) if j not in J)
                assert all(sum(E[p][j]*mu[j] for j in range(n)) == 0 for p in range(len(E)))
                if len(J) < n:
                    y = stiemke(E, n, J)
                    assert y is not None
                    Ety = [sum(E[p][j]*y[p] for p in range(len(E))) for j in range(n)]
                    assert all(x >= 0 for x in Ety) and all(Ety[j] > 0 for j in range(n) if j not in J)
                    rec.update(d="infeasible" if not J else "lower_face", y=[fs(x) for x in y], Ety=[fs(x) for x in Ety],
                               support=[S[j] for j in J], mu=[fs(x) for x in mu])
                    counts[rec["d"]] += 1; out.write(json.dumps(rec) + "\n"); continue
                # feasible with mu > 0: generic point in the cell
                Nb = nullspace(E, n)
                for _ in range(50):
                    g = [F(rng.randint(-1000, 1000), 1000) for _ in Nb]
                    eps = F(1, 4)
                    mg = [mu[j] + eps*sum(g[i]*Nb[i][j] for i in range(len(Nb))) for j in range(n)]
                    if all(x > 0 for x in mg):
                        break
                assert all(x > 0 for x in mg) and all(sum(E[p][j]*mg[j] for j in range(n)) == 0 for p in range(len(E)))
                Hg = toQ([sum(mg[j]*R[j][t] for j in range(n)) for t in range(10)])
                def ok(T):
                    if det_int(T) not in (1, -1) or mod2_perm(T) != tau:
                        return None
                    Hmu = toQ([sum(mu[j]*R[j][t] for j in range(n)) for t in range(10)])
                    if congr(T, Hmu) != Hmu:
                        return None
                    D = [[F(0)]*n for _ in range(10)]
                    for j in range(n):
                        Rj = toQ(R[j]); dv = tovec(congr(T, Rj))
                        for t in range(10):
                            D[t][j] = dv[t] - F(R[j][t])
                    return solve_left(D, E)
                T = Z = None; src = None
                key = (ci, tuple(S), tuple(tau))
                if key in claude:
                    savedT_tried += 1
                    Z = ok(claude[key]["T"])
                    if Z is not None:
                        T = claude[key]["T"]; src = "saved"; savedT_ok += 1
                if T is None:
                    for Tc in isometries(Hg, tau):
                        Z = ok(Tc)
                        if Z is not None:
                            T = Tc; src = "own"; own_T += 1; break
                if T is None:
                    rec.update(d="uncertified", mu=[fs(x) for x in mu], E=E)
                    uncert.append(rec)
                else:
                    rec.update(d="certified", mu=[fs(x) for x in mu], T=T, Tsrc=src, E=E,
                               Z=[[fs(x) for x in r] for r in Z])
                    mine_cert.add(key)
                counts[rec["d"]] += 1; out.write(json.dumps(rec) + "\n")
        print(json.dumps({"cone": ci, "done_faces": faces, "pd": pd_faces, "counts": counts, "sec": round(time.time()-t0)}), flush=True)
    out.close()
    claude_cert = {k for k, c in claude.items() if c["certified"]}
    summ = {"faces": faces, "pd_faces": pd_faces, "markings": len(taus), "cases": sum(counts.values()),
            "counts": counts, "consistency_failures": consistency_failures, "uncertified": uncert,
            "claude_cells": len(claude), "claude_certified_not_mine": sorted(map(str, claude_cert - mine_cert)),
            "mine_not_in_claude": sorted(map(str, mine_cert - set(claude))),
            "savedT_tried": savedT_tried, "savedT_passed": savedT_ok, "own_T_found": own_T,
            "wall_s": round(time.time()-t0, 1), "cpu_s": round(time.process_time()-c0, 1)}
    json.dump(summ, open(f"{outdir}/summary.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in summ.items() if k not in ("uncertified",)})[:3000])

if __name__ == "__main__":
    main()
