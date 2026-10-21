"""Independent check of docs/rank-four-stage-two-reductions.md sections 4-5 (own code only).

Universe: the 300 even rank-4 classes with odd det <= 201 stored in
reports/rank4-parity/task-b-det201.json (only their Gram matrices are read).
For each L and each of the 15 nonzero b mod 2, M = {x : b.x even}; Aut(M) is
computed from scratch (own LLL + Fincke-Pohst + backtracking). For each of the
7 pairs {v, v+w} (w = class of L in M/2M) record whether some g maps v to v+w
and whether some g swaps them.
Claims tested:
  C5: det L nonsquare and some g maps v1 -> v2  =>  some g swaps them.
  C4: every order-4 g with g^2 != -I, every 4-cycle class v: some involution R
      in Aut(M) has Rv = gv (mod 2M).
  Sq: det L nonsquare => no g in Aut(M) with g^2 = -I.

Usage (release layout): python code/rank4_stage2_review_check.py
    reports/rank4-parity/task-b-det201.json scratch/rank4-stage2/review-independent-check.json
A .json source that is present only as .json.gz is read compressed.
"""
import gzip, json, os, sys, itertools, math
from fractions import Fraction
from math import isqrt
from multiprocessing import Pool


def mat_mul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]


def T(A):
    return [list(r) for r in zip(*A)]


def gram(G, B):  # B columns = basis vectors
    return mat_mul(mat_mul(T(B), G), B)


def lll(A):
    """LLL on Gram matrix A; returns U (columns) and reduced Gram."""
    n = len(A)
    U = [[int(i == j) for j in range(n)] for i in range(n)]
    def g():
        return mat_mul(mat_mul(T(U), A), U)
    k = 1
    while k < n:
        Gm = g()
        # Gram-Schmidt
        mu = [[Fraction(0)] * n for _ in range(n)]
        Bs = [Fraction(0)] * n
        for i in range(n):
            for j in range(i):
                mu[i][j] = (Fraction(Gm[i][j]) - sum(mu[j][l] * mu[i][l] * Bs[l] for l in range(j))) / Bs[j]
            Bs[i] = Gm[i][i] - sum(mu[i][l] ** 2 * Bs[l] for l in range(i))
        changed = False
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q:
                for r in range(n):
                    U[r][k] -= q * U[r][j]
                changed = True
                break
        if changed:
            continue
        if Bs[k] < (Fraction(3, 4) - mu[k][k - 1] ** 2) * Bs[k - 1]:
            for r in range(n):
                U[r][k], U[r][k - 1] = U[r][k - 1], U[r][k]
            k = max(k - 1, 1)
        else:
            k += 1
    return U, g()


def short_vectors(A, R):
    """All integer x != 0 with x^T A x <= R (exact check)."""
    n = len(A)
    # float Cholesky-type decomposition q_ij
    q = [[float(A[i][j]) for j in range(n)] for i in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            q[j][i] = q[i][j]
            q[i][j] = q[i][j] / q[i][i]
        for k in range(i + 1, n):
            for l in range(k, n):
                q[k][l] -= q[k][i] * q[i][l]
    out = []
    x = [0] * n
    def rec(i, rem):
        c = -sum(q[i][j] * x[j] for j in range(i + 1, n))
        r = math.sqrt(max(rem, 0) / q[i][i]) + 1e-9
        for xi in range(math.ceil(c - r), math.floor(c + r) + 1):
            x[i] = xi
            nr = rem - q[i][i] * (xi - c) ** 2
            if nr < -1e-7:
                continue
            if i == 0:
                if any(x):
                    v = list(x)
                    if sum(v[a] * A[a][b] * v[b] for a in range(n) for b in range(n)) <= R:
                        out.append(v)
            else:
                rec(i - 1, nr)
        x[i] = 0
    rec(n - 1, R + 1e-7)
    return out


def automorphisms(A):
    n = len(A)
    R = max(A[i][i] for i in range(n))
    vecs = short_vectors(A, R)
    def ip(u, v):
        return sum(u[a] * A[a][b] * v[b] for a in range(n) for b in range(n))
    by_norm = {}
    for v in vecs:
        by_norm.setdefault(ip(v, v), []).append(v)
    res = []
    imgs = []
    def rec(i):
        if i == n:
            g = T(imgs)  # columns = images of basis vectors
            res.append(g)
            return
        for v in by_norm.get(A[i][i], []):
            if all(ip(imgs[j], v) == A[j][i] for j in range(i)):
                imgs.append(v)
                rec(i + 1)
                imgs.pop()
    rec(0)
    # unimodularity check
    for g in res:
        assert abs(det4(g)) == 1
        assert mat_mul(mat_mul(T(g), A), g) == A
    return res


def det4(M):
    n = len(M)
    if n == 1:
        return M[0][0]
    return sum((-1) ** j * M[0][j] * det4([r[:j] + r[j + 1:] for r in M[1:]]) for j in range(n))


def solve_int(B, y):
    """Integer x with B x = y (B square integer, invertible)."""
    n = len(B)
    M = [[Fraction(B[i][j]) for j in range(n)] + [Fraction(y[i])] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [M[r][k] - f * M[c][k] for k in range(n + 1)]
    x = [M[i][n] / M[i][i] for i in range(n)]
    assert all(v.denominator == 1 for v in x)
    return [int(v) for v in x]


def mod2(v):
    return tuple(x % 2 for x in v)


def matvec(g, v):
    return [sum(g[i][j] * v[j] for j in range(len(v))) for i in range(len(g))]


def analyse(item):
    idx, G = item
    d = det4(G)
    square = isqrt(d) ** 2 == d
    I = [[int(i == j) for j in range(4)] for i in range(4)]
    negI = [[-x for x in r] for r in I]
    stats = dict(index=idx, det=d, square=square, pairs=0, maps=0, swaps=0,
                 map_no_swap=0, sqrt_minus_I=0, mixed4=0, cycles=0, cycle_fail=0)
    for b in itertools.product((0, 1), repeat=4):
        if not any(b):
            continue
        i0 = b.index(1)
        cols = []
        for j in range(4):
            if j == i0:
                cols.append([2 * int(r == i0) for r in range(4)])
            elif b[j] == 0:
                cols.append([int(r == j) for r in range(4)])
            else:
                cols.append([int(r == j) + int(r == i0) for r in range(4)])
        Bm = T(cols)
        A0 = gram(G, Bm)
        U, A = lll(A0)
        Bm = mat_mul(Bm, U)
        assert abs(det4(Bm)) == 2 and det4(A) == 4 * d
        auts = automorphisms(A)
        w = mod2(solve_int(Bm, [2 * int(r == i0) for r in range(4)]))
        acts = [(g, {}) for g in auts]
        classes = [c for c in itertools.product((0, 1), repeat=4) if any(c) and c != w]
        for g, act in acts:
            for c in itertools.product((0, 1), repeat=4):
                act[c] = mod2(matvec(g, list(c)))
        for g, act in acts:
            g2 = mat_mul(g, g)
            if g2 == negI:
                stats['sqrt_minus_I'] += 1
            elif g2 != I and mat_mul(g2, g2) == I:
                stats['mixed4'] += 1
                for v in classes + [w]:
                    gv = act[v]
                    if gv != v and act[gv] != v:  # 4-cycle
                        stats['cycles'] += 1
                        ok = any(mat_mul(R, R) == I and actR[v] == gv for R, actR in acts)
                        if not ok:
                            stats['cycle_fail'] += 1
        seen = set()
        for v in classes:
            v2 = tuple((x + y) % 2 for x, y in zip(v, w))
            key = frozenset((v, v2))
            if key in seen:
                continue
            seen.add(key)
            stats['pairs'] += 1
            m = any(act[v] == v2 for _, act in acts)
            s = any(act[v] == v2 and act[v2] == v for _, act in acts)
            stats['maps'] += m
            stats['swaps'] += s
            if m and not s:
                stats['map_no_swap'] += 1
    return stats


if __name__ == '__main__':
    src, out = sys.argv[1], sys.argv[2]
    path = src if os.path.exists(src) else src + '.gz'
    raw = open(path, 'rb').read()
    data = json.loads(gzip.decompress(raw) if path.endswith('.gz') else raw)
    if os.path.exists(out):
        raise FileExistsError(out)
    os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
    items = [(e['index'], e['G']) for e in data['lattices']]
    with Pool(min(12, os.cpu_count() or 1)) as p:
        rows = p.map(analyse, items, chunksize=1)
    agg = {}
    for key in ('nonsquare', 'square'):
        sel = [r for r in rows if r['square'] == (key == 'square')]
        agg[key] = {k: sum(r[k] for r in sel) for k in ('pairs', 'maps', 'swaps', 'map_no_swap', 'sqrt_minus_I', 'mixed4', 'cycles', 'cycle_fail')}
        agg[key]['lattices'] = len(sel)
    agg['square_map_no_swap_dets'] = sorted({r['det'] for r in rows if r['square'] and r['map_no_swap']})
    json.dump({'summary': agg, 'rows': rows}, open(out, 'w'), indent=1)
    print(json.dumps(agg, indent=1))
