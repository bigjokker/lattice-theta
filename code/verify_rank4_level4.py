"""Independent replay of the rank-four level-4 search (second backends only).

Replays results/rank4-level4/task-a*.json (the index-two neighbour tables of
the known isospectral pairs) and examples/rank4-level4/*.json (the Conway-
Sloane bijection certificate and the level-4 configuration certificates).

Independence from rank4_level4_search.py:
  * vectors are enumerated by rectangular boxes with the exact Cauchy-Schwarz
    radii isqrt(bound * (G^-1)_ii) in the outer three coordinates and an exact
    integer quadratic in the innermost one; no LDL recursion is used;
  * isometries and automorphism groups come from verify_rank4_feasibility.
    box_isometries (integer boxes, pairwise joins);
  * every basis, Gram matrix, kernel, coset lattice, transported characteristic
    and bijection piece is recomputed from the stored integers.
A node limit raises EnumerationLimit and the replay reports "unresolved"; a
partial enumeration is never used.  Any inconsistency raises InvalidCertificate.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
from itertools import product
import json
from math import isqrt
from pathlib import Path
import time

from exact_theta import (EnumerationLimit, InvalidCertificate, dot, identity, integer, inverse,
                         ldl, matrix, matvec, multiply, require, symmetry_data, transpose, vector)
from verify_rank4_feasibility import box_isometries, det4, modular_check

CLASSES = [list(v) for v in product((0, 1), repeat=4)]
NONZERO = [c for c in CLASSES if any(c)]


def class_index(v):
    return sum((x & 1) << i for i, x in enumerate(v))


def class_vector(k):
    return [(k >> i) & 1 for i in range(4)]


def gram_of(G, B):
    return multiply(multiply(transpose(B), G), B)


def integral(A):
    return all(Fraction(x).denominator == 1 for row in A for x in row)


def coordinates(B, w):
    return matvec(inverse(B), w)


def member(B, w):
    return all(x.denominator == 1 for x in coordinates(B, w))


# ----------------------------------------------------------------------------
# box enumeration (independent of the builder's integer Fincke-Pohst)
# ----------------------------------------------------------------------------

def _box(G, bound, parity, inner, max_nodes):
    matrix(G, "G", 4)
    integer(bound, "bound")
    require(bound >= 0, "negative bound")
    integer(max_nodes, "max_nodes")
    require(max_nodes > 0, "positive node limit required")
    ldl(G)
    inv = inverse(G)
    radii = [isqrt((bound * inv[i][i]).numerator // (bound * inv[i][i]).denominator) for i in range(4)]
    ranges = []
    for i in (1, 2, 3):
        r = radii[i]
        if parity is None:
            ranges.append(range(-r, r + 1))
        else:
            ranges.append(range(-r + ((parity[i] + r) % 2), r + 1, 2))
    G00 = G[0][0]
    nodes = 0
    for y3 in ranges[2]:
        for y2 in ranges[1]:
            for y1 in ranges[0]:
                nodes += 1
                if nodes > max_nodes:
                    raise EnumerationLimit(f"box enumeration exceeded {max_nodes} nodes")
                y = [0, y1, y2, y3]
                c = G[0][1] * y1 + G[0][2] * y2 + G[0][3] * y3
                C = sum(G[j][k] * y[j] * y[k] for j in range(1, 4) for k in range(1, 4))
                T = c * c - G00 * (C - bound)
                if T < 0:
                    continue
                s = isqrt(T)
                lo = -((s + c) // G00)
                hi = (s - c) // G00
                if parity is not None:
                    lo += (parity[0] - lo) % 2
                    step = 2
                else:
                    step = 1
                if lo > hi:
                    continue
                nodes += (hi - lo) // step + 1
                if nodes > max_nodes:
                    raise EnumerationLimit(f"box enumeration exceeded {max_nodes} nodes")
                inner(lo, hi, step, c, C, y)
    return nodes


def box_class_histogram(G, bound, max_nodes=100_000_000):
    hist = [[0] * (bound + 1) for _ in range(16)]
    G00 = G[0][0]

    def inner(lo, hi, step, c, C, y):
        base = ((y[1] & 1) << 1) | ((y[2] & 1) << 2) | ((y[3] & 1) << 3)
        for y0 in range(lo, hi + 1):
            q = (G00 * y0 + 2 * c) * y0 + C
            require(0 <= q <= bound, "box inner range error")
            hist[base | (y0 & 1)][q] += 1

    nodes = _box(G, bound, None, inner, max_nodes)
    return hist, nodes


def box_signed_shells(G, a, b, bound, max_nodes=100_000_000):
    vector(a, "a", 4)
    vector(b, "b", 4)
    counts, signed = [0] * (bound + 1), [0] * (bound + 1)
    G00 = G[0][0]

    def inner(lo, hi, step, c, C, y):
        for y0 in range(lo, hi + 1, 2):
            q = (G00 * y0 + 2 * c) * y0 + C
            require(0 <= q <= bound, "box inner range error")
            x = [(y0 - a[0]) // 2] + [(y[i] - a[i]) // 2 for i in range(1, 4)]
            counts[q] += 1
            signed[q] += -1 if dot(x, b) % 2 else 1

    nodes = _box(G, bound, a, inner, max_nodes)
    return [{"norm4": q, "norm": str(Fraction(q, 4)), "vector_count": counts[q], "signed_coefficient": signed[q]}
            for q in range(bound + 1) if counts[q]], nodes


def sub_theta(hist, b):
    return [sum(hist[k][q] for k in range(16) if dot(class_vector(k), b) % 2 == 0) for q in range(len(hist[0]))]


def over_theta(hist, v):
    k = class_index(v)
    return [hist[0][q] + hist[k][q] for q in range(len(hist[0]))]


def first_difference(x, y):
    return next((i for i, (p, q) in enumerate(zip(x, y)) if p != q), None)


# ----------------------------------------------------------------------------
# independent constructions of the known pairs
# ----------------------------------------------------------------------------

def cs_gram(a, b, c, d, sigma):
    s = sigma
    rows = [[3 * s, -1, -1, -1], [1, 3 * s, 1, -1], [1, -1, 3 * s, 1], [1, 1, -1, 3 * s]]
    D = [Fraction(a, 12), Fraction(b, 12), Fraction(c, 12), Fraction(d, 12)]
    M = [[sum(D[k] * rows[i][k] * rows[j][k] for k in range(4)) for j in range(4)] for i in range(4)]
    require(integral(M), "non-integral Conway-Sloane Gram")
    return [[int(x) for x in row] for row in M]


def quaternary(coeffs):
    a11, a22, a33, a44, a12, a13, a23, a14, a24, a34 = coeffs
    return [[2 * a11, a12, a13, a14], [a12, 2 * a22, a23, a24], [a13, a23, 2 * a33, a34], [a14, a24, a34, 2 * a44]]


def expected_pair(name):
    import re
    if name.startswith("A2"):
        return ([[17, -5, 5, 2], [-5, 13, 6, 5], [5, 6, 9, 5], [2, 5, 5, 5]],
                [[17, -2, 1, 7], [-2, 13, 9, 1], [1, 9, 9, 2], [7, 1, 2, 5]])
    if name.startswith("A3"):
        return quaternary((6, 8, 10, 15, -4, 0, 5, 1, 4, 8)), quaternary((6, 8, 10, 15, -1, 6, 1, 3, -7, 2))
    m = re.search(r"\((\d+),\s*(\d+),\s*(\d+),\s*(\d+)\)", name)
    require(m is not None, "unknown pair name")
    p = tuple(int(v) for v in m.groups())
    return cs_gram(*p, 1), cs_gram(*p, -1)


def check_family_entry(G, entry, kind):
    """The stored reduced basis spans the kernel of b (sub) or 2Z^4 + Zv (over)."""
    label = vector(entry["label"], "label", 4)
    B = matrix(entry["basis"], "basis", 4)
    require(entry["gram"] == gram_of(G, B), "family Gram mismatch")
    d = abs(det4(B))
    cols = [list(c) for c in zip(*B)]
    if kind == "sub":
        require(d == 2 and all(dot(c, label) % 2 == 0 for c in cols), "basis does not span the kernel")
    else:
        require(d == 8 and all(all(x % 2 == 0 for x in c) or all((x - y) % 2 == 0 for x, y in zip(c, label))
                               for c in cols), "basis does not span 2L + Zv")
    require(entry["determinant"] == det4(entry["gram"]), "family determinant")
    ldl(entry["gram"])


def replay_task_a(doc, max_nodes, enumeration_nodes):
    require(doc["schema"] == "work7-rank4-level4-search-v1-task-a", "wrong Task A schema")
    prefix = integer(doc["theta_prefix_bound"], "prefix")
    filter_bound = integer(doc["filter_bound"], "filter bound")
    reports = []
    for r in doc["pairs"]:
        G1, G2 = matrix(r["G1"], "G1", 4), matrix(r["G2"], "G2", 4)
        E1, E2 = expected_pair(r["name"])
        require(G1 == E1 and G2 == E2, "stored pair differs from the independent construction")
        require(r["determinants"] == [det4(G1), det4(G2)], "determinants")
        h1, n1 = box_class_histogram(G1, prefix, enumeration_nodes)
        h2, n2 = box_class_histogram(G2, prefix, enumeration_nodes)
        t1 = [sum(h1[k][q] for k in range(16)) for q in range(prefix + 1)]
        t2 = [sum(h2[k][q] for k in range(16)) for q in range(prefix + 1)]
        require(r["first_theta_difference"] == first_difference(t1, t2), "theta comparison")
        require(r["theta_first_terms"] == t1[:41] and r["vectors_counted"] == [sum(t1), sum(t2)], "theta prefix")
        for key, G in (("G1", G1), ("G2", G2)):
            R = r["reduced"][key]
            require(ldl(R) and det4(R) == det4(G), "reduced Gram determinant")
            require(box_isometries(R, G, max_nodes), "reduced Gram is not isometric to the pair member")
        require(len(box_isometries(r["reduced"]["G1"], r["reduced"]["G2"], max_nodes)) == r["direct_isometries_L1_L2"],
                "direct isometry count")
        require([len(box_isometries(r["reduced"]["G1"], r["reduced"]["G1"], max_nodes)),
                 len(box_isometries(r["reduced"]["G2"], r["reduced"]["G2"], max_nodes))] == r["aut_orders"], "aut orders")
        cells_checked, iso_checked, none_checked = 0, 0, 0
        for kind in ("sub", "over"):
            fam1, fam2 = r[f"{kind}_family_L1"], r[f"{kind}_family_L2"]
            require(len(fam1) == 15 and len(fam2) == 15 and [f["label"] for f in fam1] == NONZERO
                    and [f["label"] for f in fam2] == NONZERO, "family coverage")
            thetas1, thetas2 = [], []
            for f in fam1:
                check_family_entry(G1, f, kind)
                thetas1.append(sub_theta(h1, f["label"]) if kind == "sub" else over_theta(h1, f["label"]))
            for f in fam2:
                check_family_entry(G2, f, kind)
                thetas2.append(sub_theta(h2, f["label"]) if kind == "sub" else over_theta(h2, f["label"]))
            for f, th in zip(fam1 + fam2, thetas1 + thetas2):
                require(f["theta_first_terms"] == th[:25], "family theta prefix")
            table = r[f"{kind}_table"]
            require(len(table["cells"]) == 225, "table size")
            iso_pairs = {(i, j): U for i, j, U in table["isometric_pairs"]}
            none_pairs = {tuple(p) for p in table["isospectral_nonisometric_pairs"]}
            iso_here, none_here = 0, 0
            for cell in table["cells"]:
                i, j = cell["i"], cell["j"]
                diff = first_difference(thetas1[i], thetas2[j])
                require(cell["first_theta_difference"] == diff, "cross-table theta difference")
                require(cell["agree_to_full_prefix"] == (diff is None)
                        and cell["agree_to_filter_bound"] == (diff is None or diff > filter_bound), "cross-table flags")
                cells_checked += 1
                if diff is None:
                    A, B = fam1[i]["gram"], fam2[j]["gram"]
                    found = box_isometries(A, B, max_nodes)
                    require(cell["isometry_count"] == len(found), "isometry count")
                    if found:
                        U = matrix(cell["isometry"], "isometry", 4)
                        require(abs(det4(U)) == 1 and gram_of(B, U) == A, "stored isometry")
                        require((i, j) in iso_pairs and iso_pairs[(i, j)] == U, "isometric pair list")
                        iso_here += 1
                    else:
                        require((i, j) in none_pairs, "isospectral non-isometric list")
                        none_here += 1
                else:
                    require("isometry_count" not in cell and (i, j) not in iso_pairs and (i, j) not in none_pairs,
                            "isometry data on a non-isospectral cell")
            require(len(iso_pairs) == iso_here and len(none_pairs) == none_here, "table lists")
            iso_checked += iso_here
            none_checked += none_here
        require(r["sub_matches"] == len(r["sub_table"]["isometric_pairs"])
                and r["over_matches"] == len(r["over_table"]["isometric_pairs"]), "match counts")
        require((r["sub_matches"] > 0) == (r["over_matches"] > 0), "sub/over detection must agree")
        require(r["verdict"] == ("klein_configuration_found" if r["sub_matches"] else "no_index_two_neighbour_relation"),
                "verdict")
        require(len(r["klein_configurations"]) == r["sub_matches"], "configuration count")
        for k in r["klein_configurations"]:
            check_klein_construction(G1, G2, k, max_nodes)
        reports.append({"name": r["name"], "determinants": r["determinants"], "box_nodes": [n1, n2],
                        "cells_checked": cells_checked, "isometric_cells": iso_checked,
                        "isospectral_nonisometric_cells": none_checked, "verdict": r["verdict"]})
    return reports


def check_klein_construction(G1, G2, k, max_nodes):
    """L0, O1, O2, L, L' recomputed in L2 coordinates from the stored doubled bases."""
    D = {name: matrix(B, name, 4) for name, B in k["doubled_bases"].items()}
    L0, O1, O2, L, Lp = D["L0"], D["O1"], D["O2"], D["L"], D["Lprime"]
    for name, B, d in (("L0", L0, 32), ("O1", O1, 16), ("O2", O2, 16), ("L", L, 16), ("Lprime", Lp, 8)):
        require(abs(det4(B)) == d, f"doubled {name} has the wrong index")
    cols = lambda B: [list(c) for c in zip(*B)]
    for name in ("O1", "O2", "L", "Lprime"):
        require(all(member(D[name], c) for c in cols(L0)), f"L0 not inside {name}")
    for name in ("O1", "O2", "L"):
        require(all(member(Lp, c) for c in cols(D[name])), f"{name} not inside L'")
    # the three intermediate lattices are distinct, so they are all the intermediate lattices
    require(not all(member(O1, c) for c in cols(O2)) and not all(member(O1, c) for c in cols(L))
            and not all(member(O2, c) for c in cols(L)), "intermediate lattices coincide")
    # O2 is 2 L2 itself (L2 = Z^4 in these coordinates) and L0 = 2 ker b2
    require(all(member(O2, [2 * x for x in c]) for c in identity(4)) and all(member(identity(4), [Fraction(x, 2) for x in c]) for c in cols(O2)),
            "O2 is not L2")
    G = matrix(k["G"], "G", 4)
    require(gram_of(G2, L) == [[4 * x for x in row] for row in G], "Gram of L")
    require(det4(G) == det4(G2) == k["determinant"], "determinant of L")
    g1 = gram_of(G2, O1)
    require(all(x % 4 == 0 for row in g1 for x in row), "O1 integral")
    require(box_isometries([[x // 4 for x in row] for row in g1], G1, max_nodes), "O1 is not isometric to L1")
    a, b = vector(k["a"], "a", 4), vector(k["b"], "b", 4)
    half = [[Fraction(x, 2) for x in row] for row in L]
    require(b == [0 if member(L0, [2 * half[r][i] for r in range(4)]) else 1 for i in range(4)], "b")
    require(any(a) and dot(a, b) % 2 == 0, "characteristic parity")
    # the two cosets a/2 + L0 and a/2 + t + L0 (t in L outside L0) are O1 - L0 and O2 - L0 in some order
    xa = matvec(half, a)
    t = next(i for i in range(4) if b[i])
    xt = [half[r][t] for r in range(4)]
    half_a = [x / 2 for x in xa]
    other = [x + y for x, y in zip(half_a, xt)]
    require(not member(L0, [2 * x for x in half_a]), "a lies in 2L0: trivial characteristic")
    in_O1 = member(O1, [2 * x for x in half_a]) and member(O2, [2 * x for x in other])
    in_O2 = member(O2, [2 * x for x in half_a]) and member(O1, [2 * x for x in other])
    require(in_O1 or in_O2, "the characteristic cosets are not the two coset lattices")


# ----------------------------------------------------------------------------
# certificates
# ----------------------------------------------------------------------------

def replay_bijection(doc):
    require(doc["schema"] == "work7-rank4-level4-search-v1-cs-bijection", "wrong bijection schema")
    Vp, Vm = matrix(doc["V_plus"], "V+", 4), matrix(doc["V_minus"], "V-", 4)
    for sigma, V in ((1, Vp), (-1, Vm)):
        s = sigma
        expected = transpose([[3 * s, -1, -1, -1], [1, 3 * s, 1, -1], [1, -1, 3 * s, 1], [1, 1, -1, 3 * s]])
        require(V == expected, "Conway-Sloane generators")
    M, Mp = matrix(doc["M"], "M", 4), matrix(doc["M_prime"], "M'", 4)
    cols = lambda A: [list(c) for c in zip(*A)]
    # M = 2D4 + Z(1,1,1,1): every column has coordinates of equal parity and sum 0 mod 4; index 16
    for A, residue in ((M, 0), (Mp, 2)):
        require(abs(det4(A)) == 16, "ambient index")
        for c in cols(A):
            require(len({x % 2 for x in c}) == 1, "ambient parity")
            require(sum(c) % 4 == (residue if c[0] % 2 else 0), "ambient sum condition")
    require(all(member(M, c) for c in cols(Vp)) and all(member(Mp, c) for c in cols(Vm)), "L inside M")
    threeM = [[3 * x for x in r] for r in M]
    threeMp = [[3 * x for x in r] for r in Mp]
    require(all(member(Vp, c) for c in cols(threeM)) and all(member(Vm, c) for c in cols(threeMp)), "3M inside L")
    require(abs(det4(Vp)) == 144 and abs(det4(Vm)) == 144, "index nine")
    cosets = doc["cosets"]
    require(len(cosets) == 9 and doc["pieces"] == 9, "nine pieces")
    reps = []
    for item in cosets:
        c = vector(item["representative"], "representative", 4)
        require(member(Vp, c), "representative outside L^+")
        require(all(not member(threeM, [p - q for p, q in zip(c, r)]) for r in reps), "repeated coset")
        reps.append(c)
        zeros = [i for i in range(4) if c[i] % 3 == 0]
        i = integer(item["flip"], "flip")
        require(item["mod3"] == [x % 3 for x in c], "mod 3 pattern")
        require((len(zeros) == 4 and i == 0) or zeros == [i], "flip index")
        image = [(-x if k == i else x) for k, x in enumerate(c)]
        require(item["image"] == image and member(Vm, image), "image outside L^-")
    for i in range(9):
        for j in range(i):
            require(not member(threeMp, [p - q for p, q in zip(cosets[i]["image"], cosets[j]["image"])]),
                    "images not distinct")
    return {"verdict": "verified", "pieces": 9,
            "statement": "L^+ and L^- are isospectral for every diagonal form (piecewise isometry)"}


def replay_certificate(cert, max_nodes, enumeration_nodes):
    require(cert["schema"] == "work7-rank4-level4-search-v1-certificate", "wrong certificate schema")
    G = matrix(cert["G"], "G", 4)
    _, pivots = ldl(G)
    require(cert["ldl_pivots"] == [str(p) for p in pivots], "pivots")
    det = det4(G)
    require(cert["determinant"] == det, "determinant")
    a, b = vector(cert["a"], "a", 4), vector(cert["b"], "b", 4)
    require(all(x in (0, 1) for x in a + b) and any(a) and any(b) and dot(a, b) % 2 == 0, "characteristic")
    U = matrix(cert["lll_transform_from_construction"], "U", 4)
    Gc = matrix(cert["construction_G"], "construction G", 4)
    require(abs(det4(U)) == 1 and gram_of(Gc, U) == G, "construction transport")
    ac = matvec(inverse(U), cert["construction_a"])
    require(all(x.denominator == 1 for x in ac) and [int(x) % 2 for x in ac] == a
            and [x % 2 for x in matvec(transpose(U), cert["construction_b"])] == b, "transported characteristic")
    K = matrix(cert["kernel_basis"], "K", 4)
    cols = lambda A: [list(c) for c in zip(*A)]
    require(abs(det4(K)) == 2 and all(dot(c, b) % 2 == 0 for c in cols(K)), "kernel basis")
    G0 = gram_of(G, K)
    require(cert["G0"] == G0 and cert["d0"] == det4(G0) == 4 * det, "kernel Gram")
    t = vector(cert["complement_t"], "t", 4)
    require(dot(t, b) % 2 == 1, "complement")
    for name, v in (("O_v1", a), ("O_v2", [p + 2 * q for p, q in zip(a, t)])):
        C = matrix(cert["coset_lattices"][name]["doubled_basis"], name, 4)
        require(abs(det4(C)) == 16, "coset lattice index")
        xa = coordinates(K, v)
        require(all(x.denominator == 1 for x in xa), "class outside the kernel")
        for c in cols(C):
            xc = coordinates(K, c)
            require(all(x.denominator == 1 for x in xc), "coset basis outside the kernel")
            require(all(int(x) % 2 == 0 for x in xc) or all((int(x) - int(y)) % 2 == 0 for x, y in zip(xc, xa)),
                    "coset basis outside 2L0 + Zv")
        g = gram_of(G, C)
        require(all(x % 4 == 0 for row in g for x in row), "coset lattice not integral")
        O = [[x // 4 for x in row] for row in g]
        require(cert["coset_lattices"][name]["gram"] == O and det4(O) == det, "coset Gram")
        red = cert["coset_lattices"][name]["reduced"]
        Ur = matrix(red["basis"], "reduced basis", 4)
        require(abs(det4(Ur)) == 1 and gram_of(O, Ur) == red["gram"] == matrix(red["gram"], "reduced", 4), "coset reduction")
    members = cert["members"]
    L1, L2 = matrix(members["L1"], "L1", 4), matrix(members["L2"], "L2", 4)
    E1, E2 = expected_pair(cert["pair"])
    require(L1 == E1 and L2 == E2, "members differ from the independent construction")
    for key, Gm in (("L1", L1), ("L2", L2)):
        red = members[key + "_reduced"]
        Ur = matrix(red["basis"], "member reduction", 4)
        require(abs(det4(Ur)) == 1 and gram_of(Gm, Ur) == red["gram"], "member reduction")
    links = cert["links"]
    require({links["O_v1"]["isometric_to"], links["O_v2"]["isometric_to"]} == {"L1", "L2"}, "links")
    for oname in ("O_v1", "O_v2"):
        Or = cert["coset_lattices"][oname]["reduced"]["gram"]
        Rr = members[links[oname]["isometric_to"] + "_reduced"]["gram"]
        W = matrix(links[oname]["reduced_isometry"], "link", 4)
        require(abs(det4(W)) == 1 and gram_of(Rr, W) == Or, "link isometry")
    # non-isometry: independent complete box search and the sublattice invariant
    R1, R2 = members["L1_reduced"]["gram"], members["L2_reduced"]["gram"]
    ni = cert["non_isometry"]
    require(ni["ldl_sphere"]["source_reduced_gram"] == R1 and ni["ldl_sphere"]["target_reduced_gram"] == R2
            and ni["ldl_sphere"]["isometries_found"] == 0, "non-isometry record")
    require(box_isometries(R1, R2, max_nodes) == [], "members are isometric")
    bound = integer(ni["index_two_sublattice_theta_invariant"]["bound"], "invariant bound")
    h1, _ = box_class_histogram(R1, bound, enumeration_nodes)
    h2, _ = box_class_histogram(R2, bound, enumeration_nodes)
    inv1 = sorted(sub_theta(h1, v) for v in NONZERO)
    inv2 = sorted(sub_theta(h2, v) for v in NONZERO)
    require(inv1 != inv2 and ni["index_two_sublattice_theta_invariant"]["multisets_differ"] is True, "invariant")
    require(ni["index_two_sublattice_theta_invariant"]["shared_prefixes"] == sum(1 for x in inv1 if x in inv2), "shared prefixes")
    # full automorphism group of L and stabilizer signs
    group = box_isometries(G, G, max_nodes)
    for T in cert["automorphisms"]:
        matrix(T, "automorphism", 4)
    require(cert["automorphisms"] == group, "full automorphism group")
    stabilizer = []
    for T in group:
        if any((x - y) % 2 for x, y in zip(matvec(T, a), a)):
            continue
        dual = matvec(transpose(inverse(T)), b)
        if any((x - y) % 2 for x, y in zip(dual, b)):
            continue
        stabilizer.append({"T": T, "epsilon": symmetry_data(G, a, b, T)["epsilon"]})
    require(cert["stabilizer"] == stabilizer and stabilizer and all(s["epsilon"] == 0 for s in stabilizer),
            "stabilizer")
    # automorphisms of the kernel and their action on the two classes
    kr = cert["kernel_reduced"]
    U0 = matrix(kr["U0"], "U0", 4)
    require(abs(det4(U0)) == 1 and gram_of(G0, U0) == kr["gram"], "kernel reduction")
    group0 = box_isometries(kr["gram"], kr["gram"], max_nodes)
    require(kr["automorphisms"] == group0 and cert["kernel_automorphisms"]["order"] == len(group0), "kernel group")
    v1 = coordinates(K, a)
    v2 = coordinates(K, [p + 2 * q for p, q in zip(a, t)])
    v1, v2 = [int(x) % 2 for x in v1], [int(x) % 2 for x in v2]
    require(cert["kernel_automorphisms"]["v1"] == v1 and cert["kernel_automorphisms"]["v2"] == v2, "classes")
    U0inv = inverse(U0)
    moved = 0
    for T in group0:
        T0 = multiply(multiply(U0, T), U0inv)
        require(integral(T0), "kernel transport")
        if [int(x) % 2 for x in matvec(T0, v1)] == v2:
            moved += 1
    require(moved == 0 and cert["kernel_automorphisms"]["elements_mapping_v1_to_v2"] == 0, "level 3 or lower")
    # identity: bijection route
    ident = cert["identity"]
    bij = ident["bijection"]
    routes = []
    if bij["applies"]:
        p = tuple(bij["parameters"])
        Pp, Pm = cs_gram(*p, 1), cs_gram(*p, -1)
        if "via_isometries" in bij:
            via = bij["via_isometries"]
            csr = bij["conway_sloane_reduced"]
            for sname, Gp in (("L_plus", Pp), ("L_minus", Pm)):
                Ur = matrix(csr[sname]["basis"], "cs reduction", 4)
                require(abs(det4(Ur)) == 1 and gram_of(Gp, Ur) == csr[sname]["gram"], "cs reduction")
            require({via["L1"]["conway_sloane_member"], via["L2"]["conway_sloane_member"]} == {"L_plus", "L_minus"},
                    "both members must be reached")
            for mname in ("L1", "L2"):
                W = matrix(via[mname]["reduced_isometry"], "via", 4)
                require(abs(det4(W)) == 1 and gram_of(csr[via[mname]["conway_sloane_member"]]["gram"], W)
                        == members[mname + "_reduced"]["gram"], "via isometry")
        else:
            require(L1 == Pp and L2 == Pm, "members are not the Conway-Sloane pair of the stated parameters")
        routes.append("conway-sloane bijection (cs-bijection certificate) + explicit isometries")
    # identity: modular cutoff route
    mc = ident["modular_cutoff"]
    B = modular_check(G, a, b, mc["modular"])
    require(mc["cutoff"] == B, "cutoff")
    cutoff_result = {"attempted": mc["attempted"], "cutoff": B}
    if mc["attempted"]:
        ldl_doc = mc["ldl_backend"]
        require(ldl_doc["schema"] == "work7-rank4-cutoff-v1" and ldl_doc["G"] == G and ldl_doc["a"] == a
                and ldl_doc["b"] == b and ldl_doc["run"]["norm4_bound"] == B, "cutoff certificate header")
        ev = ldl_doc["evidence"]
        require(ev["enumeration_complete"] is True and ev["enumeration"]["complete"] is True
                and ev["enumeration"]["norm4_bound"] == B, "cutoff enumeration completeness")
        shells, nodes = box_signed_shells(G, a, b, B, enumeration_nodes)
        require(ev["enumeration"]["shells"] == shells, "cutoff shells disagree with the box backend")
        require(all(s["signed_coefficient"] == 0 for s in shells), "a signed coefficient is nonzero")
        require(ev["verdict"] == "proved_zero_by_modular_cutoff", "cutoff verdict")
        digest = hashlib.sha256(json.dumps(shells).encode()).hexdigest()
        require(mc["integer_backend"]["shells_sha256"] == digest
                and mc["integer_backend"]["vectors"] == sum(s["vector_count"] for s in shells), "integer backend record")
        cutoff_result.update({"verdict": "proved_zero_by_modular_cutoff", "box_nodes": nodes,
                              "vectors": sum(s["vector_count"] for s in shells), "shells": len(shells)})
        routes.append("modular cutoff B=%d, complete signed shells by LDL, integer and box backends" % B)
    require(ident["proved_at_all_levels"] == bool(routes), "identity status")
    require(cert["level"] == 4, "level")
    return {"verdict": "verified_level_4" if routes else "verified_candidate_without_identity_proof",
            "pair": cert["pair"], "determinant": det, "G": G, "a": a, "b": b,
            "full_group_order": len(group), "kernel_group_order": len(group0),
            "identity_routes": routes, "cutoff": cutoff_result, "members_isometric": False}


def over_basis(v):
    j = next(i for i in range(4) if v[i])
    return transpose([list(v) if i == j else [2 if r == i else 0 for r in range(4)] for i in range(4)])


def replay_task_d(doc, max_nodes, enumeration_nodes, census_path=None):
    require(doc["schema"] == "work7-rank4-level4-search-v1-task-d", "wrong Task D schema")
    rows = doc["lattices"]
    require(doc["universe"]["classes"] == len(rows) and doc["configurations_per_class"] == 105, "universe size")
    if census_path is not None:
        raw = Path(census_path).read_bytes()
        require(hashlib.sha256(raw).hexdigest() == doc["universe"]["source_sha256"], "census pack hash")
        census = json.loads(raw)
        require(census["schema"] == "work7-rank4-census-v1" and census["domain"] == {"rank": 4, "determinant_bound": 24},
                "census domain")
        require([r["G0"] for r in rows] == [c["G"] for c in census["lattices"]], "universe differs from the census classes")
    tally = Counter()
    listed = {2: [], 3: [], 4: []}
    for row in rows:
        G0 = matrix(row["G0"], "G0", 4)
        d0 = det4(G0)
        require(row["determinant"] == d0 and 1 <= d0 <= 24, "class determinant")
        prefix = integer(row["prefix"], "prefix")
        hist, _ = box_class_histogram(G0, prefix, enumeration_nodes)
        group = box_isometries(G0, G0, max_nodes)
        require(row["aut_order"] == len(group), "automorphism order")
        perms = [[class_index(matvec(T, class_vector(k))) for k in range(16)] for T in group]
        seen = {(p["v1"] and class_index(p["v1"]), class_index(p["v2"])): p for p in row["zero_pairs"]}
        local = Counter()
        for k1 in range(1, 16):
            for k2 in range(k1 + 1, 16):
                zero = hist[k1] == hist[k2]
                if not zero:
                    require((k1, k2) not in seen, "non-isospectral pair listed as a zero")
                    local["nonzero"] += 1
                    continue
                require((k1, k2) in seen, "isospectral pair missing")
                p = seen[(k1, k2)]
                require(p["v3"] == class_vector(k1 ^ k2), "third class")
                swap = any(q[k1] == k2 and q[k2] == k1 for q in perms)
                maps = any(q[k1] == k2 for q in perms)
                if swap:
                    level = 1
                elif maps:
                    level = 2
                else:
                    A = gram_of(G0, over_basis(class_vector(k1)))
                    B = gram_of(G0, over_basis(class_vector(k2)))
                    level = 3 if box_isometries(A, B, max_nodes) else 4
                require(p["level"] == level and p["same_orbit"] == maps and p["swap_exists"] == swap, "level")
                g = gram_of(G0, over_basis(class_vector(k1 ^ k2)))
                require(p["L_gram_times_4"] == g, "Gram of L")
                m = 1 if all(x % 4 == 0 for r in g for x in r) else 2 if all(x % 2 == 0 for r in g for x in r) else 4
                require(p["smallest_integral_scaling"] == m and 4 * p["det_scaled_L"] == m ** 4 * d0, "scaling")
                local[f"level_{level}"] += 1
                if level >= 2:
                    listed[level].append((row["G0"], k1, k2))
        require(row["tally"] == dict(local), "per-class tally")
        tally.update(local)
    require(doc["tally"] == dict(tally), "global tally")
    for level, key in ((2, "level_2_configurations"), (3, "level_3_configurations"), (4, "level_4_candidates")):
        require(len(doc[key]) == len(listed[level]), f"{key} count")
        require([(e["G0"], class_index(e["v1"]), class_index(e["v2"])) for e in doc[key]] == listed[level], key)
    return {"verdict": "verified", "classes": len(rows), "tally": dict(tally),
            "level_3_configurations": len(listed[3]), "level_4_candidates": len(listed[4])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-d", type=Path)
    parser.add_argument("--census", type=Path)
    parser.add_argument("--task-a", type=Path, nargs="*", default=[])
    parser.add_argument("--bijection", type=Path)
    parser.add_argument("--certificates", type=Path, nargs="*", default=[])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-nodes", type=int, default=1_000_000)
    parser.add_argument("--enumeration-nodes", type=int, default=100_000_000)
    args = parser.parse_args()
    if args.output:
        require(not args.output.exists(), "choose a fresh replay path")
    start = time.perf_counter()
    result = {"verdict": "verified", "inputs": {}, "task_a": [], "bijection": None, "certificates": [], "task_d": None}
    try:
        if args.task_d:
            raw = args.task_d.read_bytes()
            result["inputs"][str(args.task_d)] = hashlib.sha256(raw).hexdigest()
            result["task_d"] = replay_task_d(json.loads(raw), args.max_nodes, args.enumeration_nodes, args.census)
        for path in args.task_a:
            raw = path.read_bytes()
            result["inputs"][str(path)] = hashlib.sha256(raw).hexdigest()
            result["task_a"].extend(replay_task_a(json.loads(raw), args.max_nodes, args.enumeration_nodes))
        if args.bijection:
            raw = args.bijection.read_bytes()
            result["inputs"][str(args.bijection)] = hashlib.sha256(raw).hexdigest()
            result["bijection"] = replay_bijection(json.loads(raw))
        for path in args.certificates:
            raw = path.read_bytes()
            result["inputs"][str(path)] = hashlib.sha256(raw).hexdigest()
            doc = json.loads(raw)
            for cert in doc["certificates"]:
                if cert["identity"]["bijection"]["applies"]:
                    require(args.bijection is not None and result["bijection"]["verdict"] == "verified",
                            "certificate relies on the bijection certificate; pass --bijection")
                result["certificates"].append(replay_certificate(cert, args.max_nodes, args.enumeration_nodes))
    except EnumerationLimit as error:
        result = {"verdict": "unresolved", "reason": str(error)}
    except (KeyError, TypeError, IndexError, ValueError, ZeroDivisionError) as error:
        raise InvalidCertificate("malformed level-4 material: " + repr(error)) from error
    result["seconds"] = time.perf_counter() - start
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k != "task_a"} | {"task_a_pairs": len(result.get("task_a", []))}, indent=1))
    return 0 if result["verdict"] == "verified" else 2


if __name__ == "__main__":
    raise SystemExit(main())
