"""Voronoi's secondary fan (L-type domains) of positive definite quaternary forms, computed exactly.

Forms are integer symmetric 4x4 matrices Q; their coordinates are q = (Q11, Q22, Q33, Q44, Q12, Q13, Q14,
Q23, Q24, Q34), and Q(v) = ev(v) . q with ev(v) = (v1^2, .., v4^2, 2 v1 v2, ..).

delaunay_star(Q):   the Delaunay simplices of Z^4 under a GENERIC rational Q that contain 0.  Candidates come
                    from a floating Delaunay triangulation of a box of lattice points (Cholesky coordinates);
                    every simplex is then checked exactly (empty circumsphere: no lattice point inside or on it
                    other than its vertices, by an exact box enumeration around the centre), and completeness
                    is checked exactly: the simplices at 0 must have total volume 5 (sum |det| = 120), since
                    every Delaunay simplex is counted once at each of its five vertices and covol Z^4 = 1.
regulators(star):   for every facet, the affine dependence of the six points (simplex + opposite vertex of the
                    neighbour) gives a linear inequality r . q > 0 that is exactly the condition that the
                    neighbour's vertex lies outside the circumsphere.  The secondary cone is {q : r . q >= 0}.
secondary_cone(Q):  irredundant inequalities and extreme rays (cdd, exact rationals).
fan():              Voronoi's algorithm: start at the principal form, cross every facet that meets the positive
                    definite cone, identify neighbours up to GL4(Z) by an exact isometry test of central forms
                    (sum of primitive rays), until no new orbit appears.  Returns the orbit representatives.

Usage: python -B rank4_secondary_fan.py fan
"""
import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from math import gcd, isqrt
from pathlib import Path
import time

import numpy as np
import cdd.gmp as cdd
from scipy.spatial import Delaunay

from exact_theta import require

OUT = Path("results/rank4-secondary-fan")
IDX = [(0, 0), (1, 1), (2, 2), (3, 3), (0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def ev(v):
    return [v[i] * v[j] * (1 if i == j else 2) for i, j in IDX]


def to_q(Q):
    return [Q[i][j] for i, j in IDX]


def to_Q(q):
    Q = [[0] * 4 for _ in range(4)]
    for (i, j), x in zip(IDX, q):
        Q[i][j] = Q[j][i] = x
    return Q


def qform(Q, v):
    return sum(Q[i][j] * v[i] * v[j] for i in range(4) for j in range(4))


def det(M):
    M = [[F(x) for x in r] for r in M]
    n, d = len(M), F(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return F(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            d = -d
        d *= M[c][c]
        for r in range(c + 1, n):
            f = M[r][c] / M[c][c]
            M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return d


def solve(A, b):
    """Exact solution of a square system."""
    n = len(A)
    M = [[F(x) for x in r] + [F(y)] for r, y in zip(A, b)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        M[c] = [x / M[c][c] for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b2 for a, b2 in zip(M[r], M[c])]
    return [M[r][n] for r in range(n)]


def positive_definite(Q):
    return all(det([r[:k] for r in Q[:k]]) > 0 for k in range(1, 5))


def primitive(vec):
    """Primitive integer vector on the ray of a rational vector."""
    fr = [F(x) for x in vec]
    L = 1
    for x in fr:
        L = L * x.denominator // gcd(L, x.denominator)
    ints = [int(x * L) for x in fr]
    g = 0
    for x in ints:
        g = gcd(g, abs(x))
    return [x // g for x in ints] if g else ints


# ----------------------------------------------------------------------------
# Delaunay star at the origin (exact)
# ----------------------------------------------------------------------------

def circumcentre(Q, verts):
    """c with Q(v - c) = Q(c) for the four nonzero vertices (0 is a vertex): 2 Q c . v = Q(v)."""
    A = [[2 * sum(Q[i][j] * v[j] for j in range(4)) for i in range(4)] for v in verts]
    return solve(A, [qform(Q, v) for v in verts])


def rational_inverse(Q):
    """Exact inverse of a nonsingular rational matrix (Gauss-Jordan over Fractions)."""
    n = len(Q)
    M = [[F(x) for x in r] + [F(int(i == j)) for j in range(n)] for i, r in enumerate(Q)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        M[c] = [x / M[c][c] for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]
                M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [r[n:] for r in M]


def _ceil_sqrt(t):
    """An integer s >= sqrt(t) for a rational t >= 0 (certified, no floating point)."""
    t = F(t)
    n = -(-t.numerator // t.denominator)          # ceil(t)
    s = isqrt(n)
    return s if s * s >= n else s + 1


def coordinate_bounds(Q, centre, r2):
    """Certified integer ranges: Q(x - c) <= r2 implies (x_i - c_i)^2 <= r2 (Q^-1)_ii (Cauchy-Schwarz)."""
    Qi = rational_inverse(Q)
    out = []
    for i in range(4):
        t = F(r2) * Qi[i][i]
        require(t >= 0, "coordinate bound: negative radius")
        s = _ceil_sqrt(t)
        c = F(centre[i])
        lo = (c - s).__floor__()
        hi = (c + s).__ceil__()
        out.append(range(lo, hi + 1))
    return out


def lattice_points_in_ellipsoid(Q, centre, r2):
    """All x in Z^4 with Q(x - centre) <= r2.  Exact: rational inverse, certified integer coordinate bounds,
    exact rational norm test (no floating point anywhere)."""
    rng = coordinate_bounds(Q, centre, r2)
    out = []
    for x0 in rng[0]:
        for x1 in rng[1]:
            for x2 in rng[2]:
                for x3 in rng[3]:
                    x = (x0, x1, x2, x3)
                    d = [F(a) - c for a, c in zip(x, centre)]
                    if qform(Q, d) <= r2:
                        out.append(x)
    return out


def delaunay_star(Q, box=2):
    """Exact Delaunay simplices containing 0 (each a frozenset of 5 lattice points), for generic Q."""
    pts = [p for p in np.ndindex(*(2 * box + 1,) * 4)]
    pts = [tuple(int(x) - box for x in p) for p in pts]
    Lc = np.linalg.cholesky(np.array(Q, dtype=float))
    tri = Delaunay(np.array(pts, dtype=float) @ Lc)
    zero = pts.index((0, 0, 0, 0))
    star = set()
    for simplex in tri.simplices:
        if zero in simplex:
            star.add(frozenset(pts[k] for k in simplex))
    good = []
    for S in star:
        verts = [v for v in S if any(v)]
        c = circumcentre(Q, verts)
        r2 = qform(Q, c)
        inside = lattice_points_in_ellipsoid(Q, c, r2)
        require(set(inside) == set(S), "delaunay: simplex not empty or Q not generic")
        good.append(S)
    vol = sum(abs(det([list(v) for v in S if any(v)])) for S in good)
    require(vol == 120, f"delaunay star incomplete (sum |det| = {vol}); enlarge the box")
    return good


def regulators(Q, star):
    """Primitive integer inequality vectors r with r . q > 0 at Q, one per facet type."""
    by_translate = {}
    for S in star:
        by_translate.setdefault(S, None)
    regs = set()
    for S in star:
        for f in combinations(sorted(S), 4):
            fs = set(f)
            p = next(v for v in S if v not in fs)
            w = f[0]
            nb = None
            for T in star:                            # neighbour across f: a translate by w of a star simplex
                T2 = frozenset(tuple(a + b for a, b in zip(t, w)) for t in T)
                if fs <= T2 and T2 != S:
                    nb = T2
                    break
            require(nb is not None, "regulators: missing neighbour")
            q = next(v for v in nb if v not in fs)
            six = list(S) + [q]
            # affine dependence: lambda over six points with sum lambda (v, 1) = 0, lambda_q = 1
            A = [[six[k][i] for k in range(5)] for i in range(4)] + [[1] * 5]
            b = [-q[i] for i in range(4)] + [-1]
            lam = solve(A, b) + [F(1)]
            r = [sum(l * e for l, e in zip(lam, col)) for col in zip(*[ev(v) for v in six])]
            r = primitive(r)
            require(sum(a * b2 for a, b2 in zip(r, to_q(Q))) > 0, "regulator must be positive at Q")
            regs.add(tuple(r))
    return sorted(regs)


def secondary_cone(Q):
    """(facet inequalities, extreme rays) of the L-type domain of generic Q, exact."""
    star = delaunay_star(Q)
    regs = regulators(Q, star)
    m = cdd.matrix_from_array([[0] + list(r) for r in regs], rep_type=cdd.RepType.INEQUALITY)
    cdd.matrix_canonicalize(m)
    facets = [primitive(row[1:]) for row in m.array]
    poly = cdd.polyhedron_from_matrix(m)
    gens = cdd.copy_generators(poly)
    rays = [primitive(row[1:]) for row in gens.array if row[0] == 0]
    require(all(row[0] == 0 for row in gens.array) or len(gens.array) == len(rays) + 1, "cone with a vertex expected")
    rays = [r for r in rays if any(r)]
    return {"facets": facets, "rays": rays, "simplices_at_0": len(star)}


def central(cone):
    s = [sum(r[k] for r in cone["rays"]) for k in range(10)]
    return to_Q(primitive(s))


def isometric(Q1, Q2):
    from verify_rank4_feasibility import box_isometries
    return bool(box_isometries(Q1, Q2, 50_000_000))


def act(U, q):
    """q-coordinates of U^t Q U."""
    Q = to_Q(q)
    return to_q([[sum(U[a][i] * Q[a][b] * U[b][j] for a in range(4) for b in range(4)) for j in range(4)] for i in range(4)])


def identify(cone, reps):
    """(index, U) with U^t R U running exactly over the rays of reps[index] as R runs over cone's rays, or None.
    U comes from an isometry of central forms; the ray-set equality is then checked exactly."""
    from verify_rank4_feasibility import box_isometries
    c = central(cone)
    for i, r in enumerate(reps):
        for U in box_isometries(r["central"], c, 50_000_000):      # U^t c U = rep central
            if sorted(tuple(primitive(act(U, x))) for x in cone["rays"]) == sorted(tuple(x) for x in r["rays"]):
                return i, U
    return None


def neighbour_point(cone, fct):
    on = [r for r in cone["rays"] if sum(a * b for a, b in zip(fct, r)) == 0]
    fpt = [sum(r[k] for r in on) for k in range(10)]
    if det(to_Q(fpt)) == 0:
        return on, fpt, None                          # facet in the boundary of the PSD cone
    cq = to_q(central(cone))
    for t in (F(1, 8), F(1, 32), F(1, 128), F(1, 512)):
        Qn = to_Q(primitive([x + t * (x - y) for x, y in zip(fpt, cq)]))
        if positive_definite(Qn) and _generic(Qn):
            return on, fpt, Qn
    raise RuntimeError("could not find a generic point across a facet")


def fan(start=None, log=print):
    """Orbit representatives of the maximal secondary cones, with a coverage certificate.

    Voronoi's algorithm.  For every facet of every representative the certificate records either
    (psd_boundary) the facet's ray sum is singular, so the facet lies in the boundary of the PSD cone; or
    (neighbour) the cone across the facet: it contains every ray of the facet, has the reversed facet
    inequality among its own facets (so it is the adjacent cone), and an integral U maps its ray set exactly
    onto the ray set of a representative.  Every facet of every representative being accounted for, the union
    of the GL4(Z)-images of the representatives is closed under crossing positive definite facets."""
    Q0 = start or [[4, -1, -1, -1], [-1, 4, -1, -1], [-1, -1, 4, -1], [-1, -1, -1, 4]]   # A4*: principal type, generic
    reps = []

    def add(cone, Q):
        cone["central"] = central(cone)
        cone["generic_point"] = Q
        reps.append(cone)
        log(f"new cone orbit {len(reps)}: {len(cone['facets'])} facets, {len(cone['rays'])} rays, central {cone['central']}")

    add(secondary_cone(Q0), Q0)
    done = 0
    while done < len(reps):
        cone = reps[done]
        cert = []
        for fct in cone["facets"]:
            on, fpt, Qn = neighbour_point(cone, fct)
            if Qn is None:
                cert.append({"facet": fct, "kind": "psd_boundary", "facet_ray_sum_det": 0})
                continue
            nb = secondary_cone(Qn)
            require(all(r in nb["rays"] for r in on), "coverage: neighbour does not contain the facet")
            require([-x for x in fct] in nb["facets"], "coverage: neighbour does not share the facet")
            hit = identify(nb, reps)
            if hit is None:
                add(nb, Qn)
                hit = identify(nb, reps)
            cert.append({"facet": fct, "kind": "neighbour", "neighbour_rays": nb["rays"], "rep": hit[0], "U": hit[1]})
        cone["coverage"] = cert
        done += 1
    return reps


def _generic(Q):
    try:
        delaunay_star(Q)
        return True
    except Exception:
        return False


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("task", choices=("fan",))
    a = ap.parse_args()
    t = time.perf_counter()
    reps = fan()
    OUT.mkdir(parents=True, exist_ok=True)
    doc = {"orbits": len(reps), "cones": reps, "seconds": round(time.perf_counter() - t, 1)}
    doc["checker"] = "exact: rational inverse, certified integer bounds, rational norms; coverage certificate per facet"
    path = OUT / "fan-exact.json"
    require(not path.exists(), f"{path} exists: refusing to overwrite a saved fan")
    path.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({"orbits": len(reps), "facets": [len(r["facets"]) for r in reps], "rays": [len(r["rays"]) for r in reps],
                      "seconds": doc["seconds"]}))


if __name__ == "__main__":
    main()
