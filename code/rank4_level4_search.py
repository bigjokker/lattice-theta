"""Rank-four level-4 search: isospectral, non-isometric index-two overlattices.

Question Q (item 2 of A_what_is_next.txt).  Let L0 be a positive definite
rank-four lattice.  Its fifteen index-two overlattices O_v = L0 + Z(v/2),
v in (L0/2L0) minus 0, come in Klein configurations: two classes v1, v2 span
the line {v1, v2, v3 = v1 + v2}, and the even characteristic on L = O_v3 whose
sign kernel is L0 vanishes identically exactly when O_v1 and O_v2 are
isospectral.  The level of such a zero is

    1  some g in Aut(L0) interchanges O_v1 and O_v2 (a symmetry witness);
    2  some g in Aut(L0) maps O_v1 onto O_v2, none interchanges them;
    3  O_v1 and O_v2 are isometric, by no element of Aut(L0);
    4  O_v1 and O_v2 are isospectral and not isometric.

Equivalently, level 4 is a pair of isospectral, non-isometric lattices that
are 2-neighbours in Kneser's sense (their intersection has index two in each).

Everything claimed is exact: integer Fincke-Pohst enumeration with exact
integer bounds derived from the rational LDL decomposition, exact LLL, and
the complete basis-image isometry search of rank4_feasibility.  Standard
library only.  See RANK4-LEVEL4-SEARCH.md for the report.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import isqrt, prod
from pathlib import Path
import time

from exact_theta import (EnumerationLimit, InvalidCertificate, complete_shells, dot,
                         identity, integer, inverse, ldl, matrix, matvec, multiply,
                         require, transpose, vector)
from rank4_feasibility import determinant, isometries as ldl_isometries, modular_data

SCHEMA = "work7-rank4-level4-search-v1"
CLASSES = [list(v) for v in product((0, 1), repeat=4)]
NONZERO = [c for c in CLASSES if any(c)]


def class_index(v):
    """Index of the class of v in Z^4/2Z^4; bit i of the index is v_i mod 2."""
    return sum((x & 1) << i for i, x in enumerate(v))


def class_vector(k):
    return [(k >> i) & 1 for i in range(4)]


# ----------------------------------------------------------------------------
# exact linear algebra helpers
# ----------------------------------------------------------------------------

def gram_of(G, B):
    """Gram matrix of the lattice spanned by the columns of B inside (Z^4, G)."""
    return multiply(multiply(transpose(B), G), B)


def is_integral_matrix(A):
    return all(Fraction(x).denominator == 1 for row in A for x in row)


def int_matrix(A):
    return [[int(Fraction(x)) for x in row] for row in A]


def lll_reduce(G):
    """Exact LLL (delta = 3/4) on a Gram matrix.  Returns (U, U^t G U)."""
    n = len(G)
    matrix(G, "G", n)
    U = identity(n)

    def gram():
        return gram_of(G, U)

    A = gram()
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            L, D = ldl(A)
            q = (L[k][j] + Fraction(1, 2)).__floor__()
            if q:
                for r in range(n):
                    U[r][k] -= q * U[r][j]
                A = gram()
        L, D = ldl(A)
        if D[k] >= (Fraction(3, 4) - L[k][k - 1] ** 2) * D[k - 1]:
            k += 1
        else:
            for r in range(n):
                U[r][k - 1], U[r][k] = U[r][k], U[r][k - 1]
            A = gram()
            k = max(k - 1, 1)
    require(abs(determinant(U)) == 1, "LLL transform not unimodular")
    return U, A


def hnf_basis(generators):
    """Basis (four integer columns) of the Z-span of integer vectors in Z^4.

    Column-style Hermite reduction by repeated Euclid; the span is unchanged.
    """
    cols = [list(g) for g in generators]
    basis = []
    for i in range(4):
        pivots = [c for c in cols if c[i] != 0]
        rest = [c for c in cols if c[i] == 0]
        while len(pivots) > 1:
            pivots.sort(key=lambda c: abs(c[i]))
            p = pivots[0]
            new = [p]
            for c in pivots[1:]:
                q = c[i] // p[i]
                c = [x - q * y for x, y in zip(c, p)]
                (new if c[i] != 0 else rest).append(c)
            pivots = new
        if not pivots:
            raise InvalidCertificate("generators do not span a rank-four lattice")
        p = pivots[0]
        if p[i] < 0:
            p = [-x for x in p]
        basis.append(p)
        cols = rest
    return transpose(basis)      # columns are the basis vectors


def solve_coordinates(B, w):
    """Rational coordinates x with B x = w."""
    return matvec(inverse(B), w)


def in_lattice(B, w):
    return all(x.denominator == 1 for x in solve_coordinates(B, w))


# ----------------------------------------------------------------------------
# exact integer Fincke-Pohst enumeration
# ----------------------------------------------------------------------------

def integer_ldl(G):
    """Leading minors Delta, integer numerators N_ji = mu_ji Delta_{i+1}, scales.

    With Delta_0 = 1 and u_i = Delta_{i+1} y_i + sum_{j>i} N_ji y_j one has
    y^t G y = sum_i u_i^2 / (Delta_i Delta_{i+1}); all quantities are integers.
    """
    L, D = ldl(G)
    n = len(G)
    Delta = [Fraction(1)]
    for d in D:
        Delta.append(Delta[-1] * d)
    require(all(x.denominator == 1 for x in Delta), "nonintegral leading minor")
    Delta = [int(x) for x in Delta]
    N = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            x = L[j][i] * Delta[i + 1]
            require(x.denominator == 1, "nonintegral LDL numerator")
            N[j][i] = int(x)
    W = prod(Delta[1:])
    scale = []
    for i in range(n):
        require(W % (Delta[i] * Delta[i + 1]) == 0, "scale integrality")
        scale.append(W // (Delta[i] * Delta[i + 1]))
    return Delta, N, W, scale


def _enumerate(G, bound, parity, inner, max_nodes):
    """Visit every y in Z^4 with y^t G y <= bound (and y = parity mod 2 if given).

    inner(lo, hi, step, c0, C, y) must process y0 in range(lo, hi+1, step) with
    y = (y0, y[1], y[2], y[3]) and norm (G00*y0 + 2*c0)*y0 + C.  Raises
    EnumerationLimit (discarding everything) past max_nodes nodes.
    """
    matrix(G, "G", 4)
    integer(bound, "bound")
    require(bound >= 0, "negative bound")
    integer(max_nodes, "max_nodes")
    require(max_nodes > 0, "positive node limit required")
    Delta, N, W, scale = integer_ldl(G)
    y = [0, 0, 0, 0]
    nodes = 0

    def level(i, Rhat):
        nonlocal nodes
        nodes += 1
        if nodes > max_nodes:
            raise EnumerationLimit(f"integer Fincke-Pohst exceeded {max_nodes} nodes")
        c = sum(N[j][i] * y[j] for j in range(i + 1, 4))
        T = Rhat // scale[i]
        if T < 0:
            return
        s = isqrt(T)
        Di = Delta[i + 1]
        lo = -((s + c) // Di)
        hi = (s - c) // Di
        if parity is not None:
            lo += (parity[i] - lo) % 2
            step = 2
        else:
            step = 1
        if lo > hi:
            return
        if i == 0:
            C = sum(G[j][k] * y[j] * y[k] for j in range(1, 4) for k in range(1, 4))
            nodes += (hi - lo) // step + 1
            if nodes > max_nodes:
                raise EnumerationLimit(f"integer Fincke-Pohst exceeded {max_nodes} nodes")
            inner(lo, hi, step, c, C, y)
            return
        for yi in range(lo, hi + 1, step):
            u = Di * yi + c
            y[i] = yi
            level(i - 1, Rhat - u * u * scale[i])

    level(3, bound * W)
    return nodes


def class_histogram(G, bound, max_nodes=50_000_000):
    """hist[k][q] = number of y in Z^4 with y^t G y = q and class k = y mod 2."""
    hist = [[0] * (bound + 1) for _ in range(16)]
    G00 = G[0][0]

    def inner(lo, hi, step, c, C, y):
        base = ((y[1] & 1) << 1) | ((y[2] & 1) << 2) | ((y[3] & 1) << 3)
        two_c = 2 * c
        for par in (0, 1):
            row = hist[base | par]
            start = lo + ((par - lo) % 2)
            for y0 in range(start, hi + 1, 2):
                row[(G00 * y0 + two_c) * y0 + C] += 1

    nodes = _enumerate(G, bound, None, inner, max_nodes)
    return hist, nodes


def signed_shells(G, a, b, bound, max_nodes=50_000_000):
    """Second backend for exact_theta.complete_shells: identical output format."""
    vector(a, "a", 4)
    vector(b, "b", 4)
    require(all(x in (0, 1) for x in a + b), "binary characteristic required")
    counts = [0] * (bound + 1)
    signed = [0] * (bound + 1)
    G00 = G[0][0]
    a0, b0 = a[0], b[0]

    def inner(lo, hi, step, c, C, y):
        sigma = sum(b[i] * ((y[i] - a[i]) >> 1) for i in range(1, 4)) & 1
        two_c = 2 * c
        if b0 == 0:
            sgn = -1 if sigma else 1
            for y0 in range(lo, hi + 1, 2):
                q = (G00 * y0 + two_c) * y0 + C
                counts[q] += 1
                signed[q] += sgn
        else:
            sgn = -1 if (sigma + ((lo - a0) >> 1)) & 1 else 1
            for y0 in range(lo, hi + 1, 2):
                q = (G00 * y0 + two_c) * y0 + C
                counts[q] += 1
                signed[q] += sgn
                sgn = -sgn

    nodes = _enumerate(G, bound, a, inner, max_nodes)
    return {"algorithm": "integer-fincke-pohst-v1", "complete": True, "norm4_bound": bound,
            "visited_nodes": nodes,
            "shells": [{"norm4": q, "norm": str(Fraction(q, 4)), "vector_count": counts[q],
                        "signed_coefficient": signed[q]} for q in range(bound + 1) if counts[q]]}


def theta_prefix(G, bound, max_nodes=50_000_000):
    hist, _ = class_histogram(G, bound, max_nodes)
    return [sum(hist[k][q] for k in range(16)) for q in range(bound + 1)]


# ----------------------------------------------------------------------------
# index-two sublattices and overlattices
# ----------------------------------------------------------------------------

def kernel_basis(b):
    """Columns spanning {x : x.b even} for a nonzero binary b (determinant 2)."""
    j = next(i for i in range(4) if b[i])
    cols = []
    for i in range(4):
        if i == j:
            cols.append([2 if r == j else 0 for r in range(4)])
        else:
            cols.append([1 if r == i else (-b[i] if r == j else 0) for r in range(4)])
    return transpose(cols)


def overlattice_basis(v):
    """Columns spanning 2Z^4 + Zv for a nonzero binary class v (determinant 8)."""
    j = next(i for i in range(4) if v[i])
    cols = []
    for i in range(4):
        cols.append(list(v) if i == j else [2 if r == i else 0 for r in range(4)])
    return transpose(cols)


def sublattice_theta(hist, b):
    """Theta prefix of the kernel of b from the class histogram of the parent."""
    return [sum(hist[k][q] for k in range(16) if dot(class_vector(k), b) % 2 == 0)
            for q in range(len(hist[0]))]


def overlattice_theta(hist, v):
    """Theta prefix of 2L + Zv from the class histogram of L."""
    k = class_index(v)
    return [hist[0][q] + hist[k][q] for q in range(len(hist[0]))]


def coset_theta(hist, v):
    """Theta prefix of the coset v + 2L (so of v/2 + L after scaling by 1/4)."""
    return list(hist[class_index(v)])


def reduced_lattice(G_parent, B):
    """LLL-reduced Gram of the lattice with columns B inside (Z^4, G_parent)."""
    U, R = lll_reduce(gram_of(G_parent, B))
    return {"basis": multiply(B, U), "gram": R, "determinant": determinant(R)}


# ----------------------------------------------------------------------------
# isometries and automorphisms
# ----------------------------------------------------------------------------

def isometry_search(G, H, max_nodes=1_000_000):
    """Every unimodular U with U^t H U = G (LDL sphere backend), sorted."""
    return ldl_isometries(G, H, max_nodes)


def mod2_action(T):
    """Permutation of the 16 classes of Z^4/2Z^4 induced by an integer matrix T."""
    return [class_index(matvec(T, class_vector(k))) for k in range(16)]


def orbits_on_classes(group):
    """Orbits of a group of integer matrices on the fifteen nonzero classes."""
    perms = [mod2_action(T) for T in group]
    seen, orbits = set(), []
    for k in range(1, 16):
        if k in seen:
            continue
        orbit, frontier = {k}, [k]
        while frontier:
            x = frontier.pop()
            for p in perms:
                if p[x] not in orbit:
                    orbit.add(p[x])
                    frontier.append(p[x])
        seen |= orbit
        orbits.append(sorted(orbit))
    return orbits, perms


# ----------------------------------------------------------------------------
# known pairs (transcribed in sources/check_cs_constructions.py and checked)
# ----------------------------------------------------------------------------

def cs_gram(a, b, c, d, sigma):
    """Conway-Sloane Theorem 1: Gram = V^t diag(a,b,c,d) V / 12 (rational)."""
    s = sigma
    rows = [[3 * s, -1, -1, -1], [1, 3 * s, 1, -1], [1, -1, 3 * s, 1], [1, 1, -1, 3 * s]]
    D = [Fraction(a, 12), Fraction(b, 12), Fraction(c, 12), Fraction(d, 12)]
    return [[sum(D[k] * rows[i][k] * rows[j][k] for k in range(4)) for j in range(4)] for i in range(4)]


def cs_pair(a, b, c, d):
    plus, minus = cs_gram(a, b, c, d, 1), cs_gram(a, b, c, d, -1)
    require(is_integral_matrix(plus) and is_integral_matrix(minus), "non-integral Conway-Sloane member")
    return int_matrix(plus), int_matrix(minus)


def quaternary_gram(coeffs):
    """Even Gram 2G of the form with coefficients (a11,a22,a33,a44,a12,a13,a23,a14,a24,a34)."""
    a11, a22, a33, a44, a12, a13, a23, a14, a24, a34 = coeffs
    return [[2 * a11, a12, a13, a14], [a12, 2 * a22, a23, a24],
            [a13, a23, 2 * a33, a34], [a14, a24, a34, 2 * a44]]


def known_pairs():
    pairs = []
    plus, minus = cs_pair(1, 7, 13, 19)
    pairs.append({"name": "A1 Conway-Sloane (1,7,13,19) = Schiemann det 1729", "G1": plus, "G2": minus,
                  "source": "CS Theorem 1, remark (iii); Earnest-Nipp p.36; LITERATURE-RANK4.md s.3,4"})
    P = [[17, -5, 5, 2], [-5, 13, 6, 5], [5, 6, 9, 5], [2, 5, 5, 5]]
    Q = [[17, -2, 1, 7], [-2, 13, 9, 1], [1, 9, 9, 2], [7, 1, 2, 5]]
    pairs.append({"name": "A2 Conway-Sloane remark (x) (a,b)=(11,1) det 961", "G1": P, "G2": Q,
                  "source": "CS remark (x) transcription in sources/check_cs_constructions.py"})
    g1 = quaternary_gram((6, 8, 10, 15, -4, 0, 5, 1, 4, 8))
    g2 = quaternary_gram((6, 8, 10, 15, -1, 6, 1, 3, -7, 2))
    pairs.append({"name": "A3 Earnest-Nipp det 84721", "G1": g1, "G2": g2,
                  "source": "Earnest-Nipp 1991 p.36 (sources/1311991000.pdf)"})
    for params in ((1, 7, 13, 43), (1, 7, 25, 31), (1, 7, 13, 31), (1, 7, 19, 25), (2, 8, 14, 20), (7, 13, 19, 49)):
        plus, minus = cs_pair(*params)
        pairs.append({"name": f"A4 Conway-Sloane {params} det {prod(params)}", "G1": plus, "G2": minus,
                      "source": "CS Theorem 1 with the transcription verified at (1,7,13,19)"})
    return pairs


# ----------------------------------------------------------------------------
# Task A: index-two neighbours of the known isospectral pairs
# ----------------------------------------------------------------------------

def index_two_family(G, hist, kind):
    """The fifteen index-two sublattices (kind 'sub') or scaled overlattices ('over')."""
    rows = []
    for v in NONZERO:
        B = kernel_basis(v) if kind == "sub" else overlattice_basis(v)
        item = reduced_lattice(G, B)
        theta = sublattice_theta(hist, v) if kind == "sub" else overlattice_theta(hist, v)
        rows.append({"label": v, "basis": item["basis"], "gram": item["gram"],
                     "determinant": item["determinant"], "theta": theta})
    return rows


def first_difference(x, y):
    return next((i for i, (p, q) in enumerate(zip(x, y)) if p != q), None)


def cross_table(fam1, fam2, filter_bound, max_nodes):
    """15 x 15 comparison; complete isometry search on every full-prefix survivor."""
    cells, isometric, new_pairs, searches = [], [], [], 0
    for i, A in enumerate(fam1):
        for j, B in enumerate(fam2):
            diff = first_difference(A["theta"], B["theta"])
            cell = {"i": i, "j": j, "first_theta_difference": diff,
                    "agree_to_filter_bound": diff is None or diff > filter_bound,
                    "agree_to_full_prefix": diff is None}
            if diff is None:
                maps = isometry_search(A["gram"], B["gram"], max_nodes)
                searches += 1
                cell["isometry_count"] = len(maps)
                if maps:
                    cell["isometry"] = maps[0]
                    isometric.append((i, j, maps[0]))
                else:
                    new_pairs.append((i, j))
            cells.append(cell)
    return {"cells": cells, "isometric_pairs": isometric, "isospectral_nonisometric_pairs": new_pairs,
            "isometry_searches": searches}


def klein_from_sublattice_match(G1, S1, G2, S2, U):
    """Klein configuration from an isometry U: (Z^4, Gram S1) -> (Z^4, Gram S2).

    S1 = {basis B1 inside L1}, S2 = {basis B2 inside L2}, U^t Gram(S2) U = Gram(S1).
    Everything is expressed in the coordinates of L2.  Returns L = third
    intermediate lattice with its even characteristic (a, b).
    """
    B1, B2 = S1["basis"], S2["basis"]
    require(multiply(multiply(transpose(U), S2["gram"]), U) == S1["gram"], "not an isometry")
    j1 = next(i for i in range(4) if S1["label"][i])
    j2 = next(i for i in range(4) if S2["label"][i])
    t1 = [int(i == j1) for i in range(4)]                   # in L1, outside S1
    t2 = [int(i == j2) for i in range(4)]                   # in L2, outside S2
    x = solve_coordinates(B1, t1)                            # S1-coordinates of t1 (half-integral)
    phi_t1 = matvec(multiply(B2, U), x)                      # image in L2-coordinates (rational)
    two_phi_t1 = [2 * v for v in phi_t1]
    require(all(v.denominator == 1 for v in two_phi_t1), "image of 2 t1 must be integral")
    two_phi_t1 = [int(v) for v in two_phi_t1]
    cols2 = transpose(B2)
    doubled = {
        "L0": hnf_basis([[2 * v for v in c] for c in cols2]),
        "O1": hnf_basis([[2 * v for v in c] for c in cols2] + [two_phi_t1]),
        "O2": hnf_basis([[2 * v for v in c] for c in cols2] + [[2 * v for v in t2]]),
        "L": hnf_basis([[2 * v for v in c] for c in cols2] + [[p + 2 * q for p, q in zip(two_phi_t1, t2)]]),
        "Lprime": hnf_basis([[2 * v for v in c] for c in cols2] + [two_phi_t1, [2 * v for v in t2]]),
    }
    grams = {}
    for name, C in doubled.items():
        g = gram_of(G2, C)
        quarter = [[Fraction(v, 4) for v in row] for row in g]
        grams[name] = quarter
    for name in ("L0", "O1", "O2", "L"):
        require(is_integral_matrix(grams[name]), f"{name} is not integral")
    GL = int_matrix(grams["L"])
    half_basis = [[Fraction(v, 2) for v in row] for row in doubled["L"]]        # basis of L in L2-coordinates
    b = [0 if in_lattice(B2, [row[i] for row in half_basis]) else 1 for i in range(4)]
    a = solve_coordinates(half_basis, two_phi_t1)
    require(all(v.denominator == 1 for v in a), "a must be integral in the basis of L")
    a = [int(v) % 2 for v in a]
    require(any(a) and any(b) and dot(a, b) % 2 == 0, "characteristic must be even and nontrivial")
    return {"coordinates": "L2 basis", "doubled_bases": doubled,
            "grams": {k: int_matrix(v) if is_integral_matrix(v) else [[str(x) for x in r] for r in v]
                      for k, v in grams.items()},
            "G": GL, "a": a, "b": b, "determinant": determinant(GL),
            "O1_isometric_to_L1_by_construction": True}


def analyse_pair(pair, prefix=2000, filter_bound=200, max_nodes=1_000_000):
    start = time.perf_counter()
    G1, G2 = pair["G1"], pair["G2"]
    d1, d2 = determinant(G1), determinant(G2)
    h1, n1 = class_histogram(G1, prefix)
    h2, n2 = class_histogram(G2, prefix)
    theta1 = [sum(h1[k][q] for k in range(16)) for q in range(prefix + 1)]
    theta2 = [sum(h2[k][q] for k in range(16)) for q in range(prefix + 1)]
    diff = first_difference(theta1, theta2)
    R1, R2 = reduced_lattice(G1, identity(4)), reduced_lattice(G2, identity(4))
    direct = isometry_search(R1["gram"], R2["gram"], max_nodes)
    aut1 = isometry_search(R1["gram"], R1["gram"], max_nodes)
    aut2 = isometry_search(R2["gram"], R2["gram"], max_nodes)
    report = {"name": pair["name"], "source": pair["source"], "G1": G1, "G2": G2,
              "determinants": [d1, d2], "even": [all(G1[i][i] % 2 == 0 for i in range(4)),
                                                 all(G2[i][i] % 2 == 0 for i in range(4))],
              "theta_prefix_bound": prefix, "enumeration_nodes": [n1, n2],
              "vectors_counted": [sum(theta1), sum(theta2)],
              "first_theta_difference": diff, "isospectral_to_prefix": diff is None,
              "theta_first_terms": theta1[:41],
              "reduced": {"G1": R1["gram"], "G2": R2["gram"]},
              "direct_isometries_L1_L2": len(direct), "aut_orders": [len(aut1), len(aut2)],
              "isometry_sphere_bound": "max diagonal of the reduced source Gram, LDL sphere"}
    fams = {("sub", 1): index_two_family(G1, h1, "sub"), ("sub", 2): index_two_family(G2, h2, "sub"),
            ("over", 1): index_two_family(G1, h1, "over"), ("over", 2): index_two_family(G2, h2, "over")}
    for kind in ("sub", "over"):
        table = cross_table(fams[(kind, 1)], fams[(kind, 2)], filter_bound, max_nodes)
        report[f"{kind}_table"] = table
        report[f"{kind}_family_L1"] = [{"label": f["label"], "basis": f["basis"], "gram": f["gram"],
                                       "determinant": f["determinant"], "theta_first_terms": f["theta"][:25]}
                                      for f in fams[(kind, 1)]]
        report[f"{kind}_family_L2"] = [{"label": f["label"], "basis": f["basis"], "gram": f["gram"],
                                       "determinant": f["determinant"], "theta_first_terms": f["theta"][:25]}
                                      for f in fams[(kind, 2)]]
    report["filter_bound"] = filter_bound
    report["klein_configurations"] = []
    for (i, j, U) in report["sub_table"]["isometric_pairs"]:
        report["klein_configurations"].append(
            klein_from_sublattice_match(G1, fams[("sub", 1)][i], G2, fams[("sub", 2)][j], U))
    report["sub_matches"] = len(report["sub_table"]["isometric_pairs"])
    report["over_matches"] = len(report["over_table"]["isometric_pairs"])
    report["verdict"] = ("klein_configuration_found" if report["klein_configurations"]
                         else "no_index_two_neighbour_relation")
    report["seconds"] = time.perf_counter() - start
    return report


def task_a(prefix=2000, filter_bound=200, max_nodes=1_000_000, names=None):
    reports = []
    for pair in known_pairs():
        if names and not any(pair["name"].startswith(n) for n in names):
            continue
        reports.append(analyse_pair(pair, prefix, filter_bound, max_nodes))
    return {"schema": SCHEMA + "-task-a", "theta_prefix_bound": prefix, "filter_bound": filter_bound,
            "max_nodes": max_nodes, "pairs": reports,
            "summary": [{"name": r["name"], "determinants": r["determinants"], "isospectral_to_prefix": r["isospectral_to_prefix"],
                         "direct_isometries": r["direct_isometries_L1_L2"], "aut_orders": r["aut_orders"],
                         "sub_matches": r["sub_matches"], "over_matches": r["over_matches"],
                         "sub_isospectral_nonisometric": len(r["sub_table"]["isospectral_nonisometric_pairs"]),
                         "over_isospectral_nonisometric": len(r["over_table"]["isospectral_nonisometric_pairs"]),
                         "verdict": r["verdict"], "seconds": round(r["seconds"], 2)} for r in reports]}


# ----------------------------------------------------------------------------
# Task C: certification
# ----------------------------------------------------------------------------

def cs_vectors(sigma):
    """Columns v_inf, v_0, v_1, v_2 of Conway-Sloane Theorem 1 in e-coordinates."""
    s = sigma
    return transpose([[3 * s, -1, -1, -1], [1, 3 * s, 1, -1], [1, -1, 3 * s, 1], [1, 1, -1, 3 * s]])


def cs_bijection_certificate():
    """Exact piecewise-isometry proof that L^+ and L^- are isospectral.

    L^+ = V^+ Z^4 lies between 3M and M, where M = 2D4 + Z(1,1,1,1); likewise
    L^- between 3M' and M' with M' = 2D4 + Z(1,1,1,-1).  Every nonzero coset of
    3M in L^+ has a representative with exactly one coordinate divisible by 3;
    changing the sign of that coordinate is an isometry of every diagonal form
    and maps the coset onto a coset of 3M' in L^-.  The nine images are distinct,
    so the pieces assemble to a norm-preserving bijection L^+ -> L^- for every
    diagonal form, in particular for diag(a,b,c,d)/12.
    """
    Vp, Vm = cs_vectors(1), cs_vectors(-1)
    D4 = [[2, -2, 0, 0], [0, 2, -2, 0], [0, 0, 2, -2], [0, 0, 2, 2]]
    M, Mp = hnf_basis(D4 + [[1, 1, 1, 1]]), hnf_basis(D4 + [[1, 1, 1, -1]])
    cols = lambda A: [list(c) for c in zip(*A)]
    require(abs(determinant(Vp)) == 144 and abs(determinant(Vm)) == 144, "wrong Conway-Sloane index")
    require(abs(determinant(M)) == 16 and abs(determinant(Mp)) == 16, "wrong ambient index")
    require(all(in_lattice(M, c) for c in cols(Vp)) and all(in_lattice(Mp, c) for c in cols(Vm)),
            "Conway-Sloane lattices not inside the ambient lattices")
    threeM = [[3 * x for x in r] for r in M]
    threeMp = [[3 * x for x in r] for r in Mp]
    require(all(in_lattice(Vp, c) for c in cols(threeM)) and all(in_lattice(Vm, c) for c in cols(threeMp)),
            "3M not inside L^+ or 3M' not inside L^-")
    reps = []
    for u in product(range(3), repeat=4):
        x = matvec(Vp, list(u))
        if not any(in_lattice(threeM, [p - q for p, q in zip(x, r)]) for r in reps):
            reps.append(x)
    require(len(reps) == 9, "L^+/3M must have order nine")
    cosets = []
    for c in reps:
        zeros = [i for i in range(4) if c[i] % 3 == 0]
        require(len(zeros) in (1, 4), "codeword weight must be three")
        i = 0 if len(zeros) == 4 else zeros[0]
        image = [(-x if k == i else x) for k, x in enumerate(c)]
        require(in_lattice(Vm, image), "flipped representative outside L^-")
        cosets.append({"representative": c, "mod3": [x % 3 for x in c], "flip": i, "image": image})
    for i in range(9):
        for j in range(i):
            require(not in_lattice(threeMp, [p - q for p, q in zip(cosets[i]["image"], cosets[j]["image"])]),
                    "two pieces map onto the same coset of 3M'")
    return {"schema": SCHEMA + "-cs-bijection", "V_plus": Vp, "V_minus": Vm, "M": M, "M_prime": Mp,
            "cosets": cosets, "pieces": 9,
            "statement": "theta of L^+(a,b,c,d) equals theta of L^-(a,b,c,d) for every diagonal form, at every norm",
            "proof": "piecewise coordinate sign change; each piece is an isometry of diag(a,b,c,d)/12 "
                     "mapping a coset of 3M in L^+ onto a coset of 3M' in L^-; images are distinct"}


def cs_parameters(name):
    import re
    m = re.search(r"\((\d+),\s*(\d+),\s*(\d+),\s*(\d+)\)", name)
    return tuple(int(v) for v in m.groups()) if m else None


def transport_characteristic(U, a, b):
    """Characteristic in the basis x = U x': a' = U^-1 a (primal), b' = U^t b (dual)."""
    a_new = matvec(inverse(U), a)
    require(all(v.denominator == 1 for v in a_new), "non-integral transported a")
    return [int(v) % 2 for v in a_new], [v % 2 for v in matvec(transpose(U), b)]


def sublattice_invariant(G, bound=40):
    """Sorted multiset of theta prefixes of the fifteen index-two sublattices."""
    hist, _ = class_histogram(G, bound)
    return sorted(sublattice_theta(hist, b) for b in NONZERO)


def certify_configuration(pair, klein, cutoff_vector_limit=30_000_000, max_nodes=1_000_000,
                          enumeration_nodes=2_000_000_000):
    start = time.perf_counter()
    U, G = lll_reduce(klein["G"])
    a, b = transport_characteristic(U, klein["a"], klein["b"])
    require(any(a) and any(b) and dot(a, b) % 2 == 0, "even nontrivial characteristic required")
    _, pivots = ldl(G)
    det = determinant(G)
    K = kernel_basis(b)
    G0 = gram_of(G, K)
    require(determinant(G0) == 4 * det, "kernel determinant")
    j = next(i for i in range(4) if b[i])
    t = [int(i == j) for i in range(4)]
    cols = [list(c) for c in zip(*K)]
    C1 = hnf_basis([[2 * v for v in c] for c in cols] + [a])
    C2 = hnf_basis([[2 * v for v in c] for c in cols] + [[p + 2 * q for p, q in zip(a, t)]])
    O = []
    for C in (C1, C2):
        g = gram_of(G, C)
        require(all(v % 4 == 0 for row in g for v in row), "coset lattice not integral")
        O.append([[v // 4 for v in row] for row in g])
    require(determinant(O[0]) == det and determinant(O[1]) == det, "coset lattice determinants")
    G1, G2 = pair["G1"], pair["G2"]
    R1, R2 = reduced_lattice(G1, identity(4)), reduced_lattice(G2, identity(4))
    O1r, O2r = reduced_lattice(O[0], identity(4)), reduced_lattice(O[1], identity(4))
    # identify the coset lattices with the known pair (either order)
    links = {}
    for oname, Or in (("O_v1", O1r), ("O_v2", O2r)):
        for pname, Rr, Gp in (("L1", R1, G1), ("L2", R2, G2)):
            maps = isometry_search(Or["gram"], Rr["gram"], max_nodes)
            if maps:
                links[oname] = {"isometric_to": pname, "reduced_isometry": maps[0], "isometry_count": len(maps)}
                break
        require(oname in links, f"{oname} is not isometric to either member")
    require({links["O_v1"]["isometric_to"], links["O_v2"]["isometric_to"]} == {"L1", "L2"},
            "the two coset lattices must realise the two members")
    # non-isometry of the pair: complete LDL sphere search plus an exact invariant
    cross = isometry_search(R1["gram"], R2["gram"], max_nodes)
    require(not cross, "members are isometric: not level 4")
    from lattice_census import sphere
    sphere_bound = max(R1["gram"][i][i] for i in range(4))
    vectors = sphere(R2["gram"], sphere_bound, max_nodes)
    counts = Counter(dot(v, matvec(R2["gram"], v)) for v in vectors)
    inv1, inv2 = sublattice_invariant(R1["gram"]), sublattice_invariant(R2["gram"])
    non_isometry = {"ldl_sphere": {"source_reduced_gram": R1["gram"], "target_reduced_gram": R2["gram"],
                                   "sphere_bound": sphere_bound, "vectors_in_sphere": len(vectors),
                                   "basis_image_candidates": {str(R1["gram"][i][i]): counts[R1["gram"][i][i]] for i in range(4)},
                                   "isometries_found": 0},
                    "index_two_sublattice_theta_invariant": {"bound": 40, "multisets_differ": inv1 != inv2,
                                                             "shared_prefixes": sum(1 for x in inv1 if x in inv2)}}
    require(inv1 != inv2, "sublattice invariant does not separate the pair")
    # full automorphism group of L and the characteristic stabilizer
    from exact_theta import symmetry_data
    group = isometry_search(G, G, max_nodes)
    stabilizer = []
    for T in group:
        if any((x - y) % 2 for x, y in zip(matvec(T, a), a)):
            continue
        dual = matvec(transpose(inverse(T)), b)
        if any((x - y) % 2 for x, y in zip(dual, b)):
            continue
        stabilizer.append({"T": T, "epsilon": symmetry_data(G, a, b, T)["epsilon"]})
    require(stabilizer and all(s["epsilon"] == 0 for s in stabilizer), "stabilizer sign inconsistent with level 4")
    # automorphisms of the sign kernel and their action on the two classes
    U0, G0r = lll_reduce(G0)
    group0 = isometry_search(G0r, G0r, max_nodes)
    v1 = solve_coordinates(K, a)
    v2 = solve_coordinates(K, [p + 2 * q for p, q in zip(a, t)])
    require(all(v.denominator == 1 for v in v1 + v2), "classes must lie in the kernel")
    v1, v2 = [int(v) % 2 for v in v1], [int(v) % 2 for v in v2]
    U0inv = inverse(U0)
    moved = 0
    for T in group0:
        T0 = multiply(multiply(U0, T), U0inv)
        require(is_integral_matrix(T0), "kernel automorphism transport")
        T0 = int_matrix(T0)
        if [x % 2 for x in matvec(T0, v1)] == v2:
            moved += 1
    require(moved == 0, "an automorphism of L0 maps O_v1 onto O_v2: not level 4")
    orbits0, _ = orbits_on_classes([int_matrix(multiply(multiply(U0, T), U0inv)) for T in group0])
    # identity proofs: (a) the Conway-Sloane bijection, directly or through isometries
    params = cs_parameters(pair["name"])
    bijection = {"applies": params is not None, "parameters": params}
    if params is None:
        for candidate in ((7, 13, 19, 49),):
            Pp, Pm = cs_pair(*candidate)
            Rp, Rm = reduced_lattice(Pp, identity(4)), reduced_lattice(Pm, identity(4))
            found = {}
            for mname, Rr in (("L1", R1), ("L2", R2)):
                for sname, Sr in (("L_plus", Rp), ("L_minus", Rm)):
                    maps = isometry_search(Rr["gram"], Sr["gram"], max_nodes)
                    if maps:
                        found[mname] = {"conway_sloane_member": sname, "reduced_isometry": maps[0]}
                        break
            if len(found) == 2 and found["L1"]["conway_sloane_member"] != found["L2"]["conway_sloane_member"]:
                bijection = {"applies": True, "parameters": candidate, "via_isometries": found,
                             "conway_sloane_reduced": {"L_plus": Rp, "L_minus": Rm}}
                break
    # (b) the modular cutoff with two enumeration backends
    item = modular_data(G, a, b)
    B = item["cutoff"]
    estimate = (3.141592653589793 ** 2 / 2) * B * B / (16 * det ** 0.5)
    cutoff = {"cutoff": B, "modular": {k: v for k, v in item.items()}, "estimated_vectors": estimate}
    if estimate <= cutoff_vector_limit:
        from rank4_feasibility import cutoff_build
        t0 = time.perf_counter()
        ldl_doc = cutoff_build(G, a, b, max_nodes=enumeration_nodes)
        t1 = time.perf_counter()
        mine = signed_shells(G, a, b, B, enumeration_nodes)
        t2 = time.perf_counter()
        require(ldl_doc["evidence"]["enumeration_complete"] is True, "LDL cutoff enumeration interrupted")
        require(ldl_doc["evidence"]["enumeration"]["shells"] == mine["shells"], "backends disagree")
        require(ldl_doc["evidence"]["verdict"] == "proved_zero_by_modular_cutoff", "cutoff did not prove the zero")
        cutoff.update({"attempted": True, "ldl_backend": ldl_doc,
                       "integer_backend": {"algorithm": mine["algorithm"], "visited_nodes": mine["visited_nodes"],
                                           "shells_sha256": hashlib.sha256(json.dumps(mine["shells"]).encode()).hexdigest(),
                                           "vectors": sum(s["vector_count"] for s in mine["shells"]),
                                           "all_signed_coefficients_zero": all(s["signed_coefficient"] == 0 for s in mine["shells"])},
                       "seconds": {"ldl": t1 - t0, "integer": t2 - t1}})
    else:
        cutoff.update({"attempted": False, "reason": f"estimated {estimate:.3g} vectors exceeds limit {cutoff_vector_limit}"})
    certificate = {"schema": SCHEMA + "-certificate", "pair": pair["name"], "source": pair["source"],
                   "parameters": params, "G": G, "a": a, "b": b, "determinant": det,
                   "ldl_pivots": [str(p) for p in pivots], "lll_transform_from_construction": U,
                   "construction_G": klein["G"], "construction_a": klein["a"], "construction_b": klein["b"],
                   "kernel_basis": K, "G0": G0, "d0": determinant(G0), "complement_t": t,
                   "coset_lattices": {"O_v1": {"doubled_basis": C1, "gram": O[0], "reduced": O1r},
                                      "O_v2": {"doubled_basis": C2, "gram": O[1], "reduced": O2r}},
                   "members": {"L1": G1, "L2": G2, "L1_reduced": R1, "L2_reduced": R2},
                   "links": links, "non_isometry": non_isometry,
                   "identity": {"bijection": bijection, "modular_cutoff": cutoff,
                                "proved_at_all_levels": bool(bijection["applies"] or cutoff["attempted"])},
                   "automorphisms": group, "stabilizer": stabilizer,
                   "kernel_reduced": {"U0": U0, "gram": G0r, "automorphisms": group0},
                   "kernel_automorphisms": {"order": len(group0), "orbits_on_classes": orbits0,
                                            "v1": v1, "v2": v2, "elements_mapping_v1_to_v2": moved},
                   "level": 4, "seconds": time.perf_counter() - start}
    return certificate


def task_c(task_a_doc, cutoff_vector_limit=30_000_000, max_nodes=1_000_000, names=None):
    certificates = []
    for report in task_a_doc["pairs"]:
        if names and not any(report["name"].startswith(n) for n in names):
            continue
        pair = {"name": report["name"], "source": report["source"], "G1": report["G1"], "G2": report["G2"]}
        for klein in report["klein_configurations"]:
            certificates.append(certify_configuration(pair, klein, cutoff_vector_limit, max_nodes))
    return certificates


# ----------------------------------------------------------------------------
# Task D: bounded exhaustive scan over the accepted census domain
# ----------------------------------------------------------------------------

def census_representatives(path):
    """The 169 isometry classes of integral rank-four lattices with det <= 24.

    Taken from the accepted census pack (schema work7-rank4-census-v1), whose
    coverage and uniqueness were replayed by verify_rank4_census.py.
    """
    raw = Path(path).read_bytes()
    doc = json.loads(raw)
    require(doc["schema"] == "work7-rank4-census-v1" and doc["domain"] == {"rank": 4, "determinant_bound": 24},
            "wrong census pack")
    rows = doc["lattices"]
    require(doc["summary"]["lattice_classes"] == len(rows), "census summary mismatch")
    for row in rows:
        matrix(row["G"], "census Gram", 4)
        require(determinant(row["G"]) == row["determinant"], "census determinant")
    return [(row["G"], row["determinant"]) for row in rows], hashlib.sha256(raw).hexdigest()


def analyse_kernel(G0, prefix=60, max_nodes=1_000_000):
    """Classify all 105 Klein configurations over the sign kernel L0 = (Z^4, G0).

    For classes v1 != v2 the configuration is a zero iff the cosets v1/2 + L0
    and v2/2 + L0 are isospectral; equality through the prefix bound is
    necessary, a difference is an exact proof of non-vanishing.  Levels 1-3 are
    exact (an isometry proves the identity); level 4 is a candidate only.
    """
    d0 = determinant(G0)
    hist, nodes = class_histogram(G0, prefix, max_nodes)
    group = isometry_search(G0, G0, max_nodes)
    orbits, perms = orbits_on_classes(group)
    pairs, tally = [], Counter()
    reduced = {}
    for k1 in range(1, 16):
        for k2 in range(k1 + 1, 16):
            if hist[k1] != hist[k2]:
                tally["nonzero"] += 1
                continue
            k3 = k1 ^ k2
            swap = any(p[k1] == k2 and p[k2] == k1 for p in perms)
            maps = any(p[k1] == k2 for p in perms)
            if swap:
                level = 1
            elif maps:
                level = 2
            else:
                for k in (k1, k2):
                    if k not in reduced:
                        reduced[k] = reduced_lattice(G0, overlattice_basis(class_vector(k)))
                iso = isometry_search(reduced[k1]["gram"], reduced[k2]["gram"], max_nodes)
                level = 3 if iso else 4
            C = overlattice_basis(class_vector(k3))
            g = gram_of(G0, C)
            m = 1 if all(x % 4 == 0 for row in g for x in row) else 2 if all(x % 2 == 0 for row in g for x in row) else 4
            scaled_det = Fraction(m ** 4 * d0, 4)
            require(scaled_det.denominator == 1, "scaled determinant")
            pairs.append({"v1": class_vector(k1), "v2": class_vector(k2), "v3": class_vector(k3), "level": level,
                          "same_orbit": maps, "swap_exists": swap, "L_gram_times_4": g,
                          "smallest_integral_scaling": m, "det_scaled_L": int(scaled_det)})
            tally[f"level_{level}"] += 1
    return {"G0": G0, "determinant": d0, "prefix": prefix, "enumeration_nodes": nodes,
            "aut_order": len(group), "orbits_on_classes": orbits, "zero_pairs": pairs, "tally": dict(tally)}


def task_d(census_path="examples/census-rank4/det24-bound32.json", prefix=60, max_nodes=1_000_000):
    start = time.perf_counter()
    reps, digest = census_representatives(census_path)
    rows, tally = [], Counter()
    level2, level3, level4 = [], [], []
    for G0, d0 in reps:
        row = analyse_kernel(G0, prefix, max_nodes)
        rows.append(row)
        tally.update(row["tally"])
        for p in row["zero_pairs"]:
            entry = {"G0": G0, "det_L0": d0, "det_L": str(Fraction(d0, 4)), **p}
            if p["level"] == 2:
                level2.append(entry)
            elif p["level"] == 3:
                level3.append(entry)
            elif p["level"] == 4:
                level4.append(entry)
    return {"schema": SCHEMA + "-task-d", "universe": {"description": "all isometry classes of positive definite "
            "integral rank-four lattices L0 with det(L0) <= 24, as sign kernels", "classes": len(reps),
            "source": census_path, "source_sha256": digest, "dedup": "exact isometry (accepted census)",
            "isospectrality_filter": f"coset theta prefixes through norm {prefix}",
            "isometry_search": "complete LDL sphere basis-image search (rank4_feasibility.isometries)"},
            "configurations_per_class": 105, "tally": dict(tally), "lattices": rows,
            "level_2_configurations": level2, "level_3_configurations": level3, "level_4_candidates": level4,
            "seconds": time.perf_counter() - start}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", choices=("A", "C", "D", "bijection"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--task-a", type=Path, help="Task A result file for task C")
    parser.add_argument("--census", type=Path, default=Path("examples/census-rank4/det24-bound32.json"))
    parser.add_argument("--prefix", type=int, default=2000)
    parser.add_argument("--filter-bound", type=int, default=200)
    parser.add_argument("--max-nodes", type=int, default=1_000_000)
    parser.add_argument("--cutoff-limit", type=float, default=30_000_000)
    parser.add_argument("--names", nargs="*")
    args = parser.parse_args()
    require(not args.output.exists(), "choose a fresh output path")
    if args.task == "A":
        doc = task_a(args.prefix, args.filter_bound, args.max_nodes, args.names)
        summary = doc["summary"]
    elif args.task == "bijection":
        doc = cs_bijection_certificate()
        summary = {"pieces": doc["pieces"], "statement": doc["statement"]}
    elif args.task == "D":
        doc = task_d(str(args.census), 60 if args.prefix == 2000 else args.prefix, args.max_nodes)
        summary = {"classes": doc["universe"]["classes"], "tally": doc["tally"],
                   "level_2": len(doc["level_2_configurations"]), "level_3": len(doc["level_3_configurations"]),
                   "level_4_candidates": len(doc["level_4_candidates"]), "seconds": round(doc["seconds"], 1)}
    else:
        task_a_doc = json.loads(args.task_a.read_text(encoding="utf-8"))
        certs = task_c(task_a_doc, args.cutoff_limit, args.max_nodes, args.names)
        doc = {"schema": SCHEMA + "-certificates", "certificates": certs}
        summary = [{"pair": c["pair"], "determinant": c["determinant"], "G": c["G"], "a": c["a"], "b": c["b"],
                    "aut_order": len(c["automorphisms"]), "aut_L0_order": c["kernel_automorphisms"]["order"],
                    "cutoff": c["identity"]["modular_cutoff"]["cutoff"],
                    "cutoff_attempted": c["identity"]["modular_cutoff"]["attempted"],
                    "cutoff_verdict": c["identity"]["modular_cutoff"].get("ldl_backend", {}).get("evidence", {}).get("verdict"),
                    "seconds": round(c["seconds"], 1)} for c in certs]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
