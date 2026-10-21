"""Smaller rank-four counterexamples from the Eisenstein construction (item 3).

Family F (brief `Fable-rank4-eisenstein-search.md`, section 3): binary
Hermitian Z[w]-lattices with integral trace form, parameters (p, q, delta),
delta = m + n w, trace-form Gram

    G0 = [[2p, -p,   n,  -m ],
          [-p, 2p,  m-n,  n ],
          [ n, m-n, 2q,  -q ],
          [-m,  n,  -q,  2q ]]        on the Z-basis (e1, w e1, e2, w e2),

with det G0 = (3pq - N(delta))^2 and the order-three isometry R (multiplication
by w).  The line Z[w]e1 modulo 2 gives three index-two overlattices
O_v1, O_v2, O_v3 (v1 = e1, v2 = w e1, v3 = e1 + w e1) permuted cyclically by R,
hence isospectral; they are integral iff p is even and delta is in 2Z[w].
Then L = O_v3 with the even characteristic whose sign kernel is L0 and whose
cosets are v1/2 + L0, v2/2 + L0 has an identically vanishing theta function,
det L = D'^2 with D' = 3aq - 2N(beta) (p = 2a, delta = 2 beta).  The zero is a
counterexample to the rank-four symmetry converse exactly when no element of
Aut(L0) swaps v1 and v2 modulo 2L0 (level 2 of the hierarchy of
LITERATURE-RANK4.md section 1).

Tasks (subcommands): A formulas and control, B complete enumeration of F for
5 <= D' <= 41, C certification, D extension to D' <= 154, E the mixed-type
lemma.  Standard library only; exact arithmetic for every claim.  Reuses the
accepted modules rank4_level4_search, rank4_feasibility, exact_theta and
rank4_counterexample without modifying them.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path
import random
import time

from exact_theta import (EnumerationLimit, InvalidCertificate, complete_shells, dot, identity,
                         integer, inverse, ldl, matrix, matvec, multiply, require, symmetry_data,
                         transpose, vector)
from rank4_feasibility import cutoff_build, determinant, isometries, modular_data
from rank4_level4_search import (analyse_kernel, class_vector, class_index, gram_of, hnf_basis,
                                 kernel_basis, lll_reduce, orbits_on_classes, overlattice_basis,
                                 reduced_lattice, signed_shells, transport_characteristic)
import rank4_counterexample

SCHEMA = "work7-rank4-eisenstein-v1"
R = [[0, -1, 0, 0], [1, -1, 0, 0], [0, 0, 0, -1], [0, 0, 1, -1]]
V1, V2, V3 = [1, 0, 0, 0], [0, 1, 0, 0], [1, 1, 0, 0]
OMEGA_POWERS = [(1, 0), (0, 1), (-1, -1)]            # 1, w, w^2 as (m, n)
LINE_BASIS_CHANGES = [((1, 0), (0, 0), (0, 0), (1, 0)),      # line (1,0): e1' = f1
                      ((1, 0), (1, 0), (0, 0), (1, 0)),      # line (1,1): e1' = f1 + f2
                      ((1, 0), (0, 1), (0, 0), (1, 0)),      # line (1,w): e1' = f1 + w f2
                      ((1, 0), (-1, -1), (0, 0), (1, 0)),    # line (1,w^2): e1' = f1 + w^2 f2
                      ((0, 0), (1, 0), (1, 0), (0, 0))]      # line (0,1): e1' = f2, e2' = f1
LINE_NAMES = ["(1,0)", "(1,1)", "(1,w)", "(1,w^2)", "(0,1)"]


# ----------------------------------------------------------------------------
# Eisenstein arithmetic (elements as (m, n) = m + n w, w^2 = -1 - w)
# ----------------------------------------------------------------------------

def emul(x, y):
    (m1, n1), (m2, n2) = x, y
    return (m1 * m2 - n1 * n2, m1 * n2 + n1 * m2 - n1 * n2)


def econj(x):
    m, n = x
    return (m - n, -n)


def enorm(x):
    m, n = x
    return m * m - m * n + n * n


def mult_matrix(x):
    """Integer matrix of multiplication by x on the Z-basis (e, w e)."""
    m, n = x
    return [[m, -n], [n, m - n]]


def trace_gram(p, q, m, n):
    return [[2 * p, -p, n, -m], [-p, 2 * p, m - n, n], [n, m - n, 2 * q, -q], [-m, n, -q, 2 * q]]


def params_from_gram(G):
    """(p, q, m, n) of a trace-form Gram; raises if G does not have the shape."""
    matrix(G, "G0", 4)
    require(G[0][0] % 2 == 0 and G[2][2] % 2 == 0, "trace form has odd diagonal")
    p, q, n, m = G[0][0] // 2, G[2][2] // 2, G[0][2], -G[0][3]
    require(trace_gram(p, q, m, n) == G, "not a trace-form Gram matrix")
    return p, q, m, n


def basis_change(x, y, z, w):
    """Z-matrix of the Z[w]-linear change e1' = x f1 + y f2, e2' = z f1 + w f2."""
    X, Y, Z, W = mult_matrix(x), mult_matrix(y), mult_matrix(z), mult_matrix(w)
    return [X[0] + Z[0], X[1] + Z[1], Y[0] + W[0], Y[1] + W[1]]


def is_unit(x):
    return enorm(x) == 1


def line_parameters(p, q, m, n):
    """Whether the line Z[w]e1 has integral overlattices, and then (a, q, beta, D')."""
    if p % 2 or m % 2 or n % 2:
        return None
    a, beta = p // 2, (m // 2, n // 2)
    Dp = 3 * a * q - 2 * enorm(beta)
    return {"a_param": a, "q": q, "beta": list(beta), "D_prime": Dp}


# ----------------------------------------------------------------------------
# Task A: formulas and control
# ----------------------------------------------------------------------------

def overlattice_gram(G0, v):
    """Gram of O_v = L0 + Z(v/2); raises if it is not integral."""
    g = gram_of(G0, overlattice_basis(v))
    require(all(x % 4 == 0 for row in g for x in row), "overlattice not integral")
    return [[x // 4 for x in row] for row in g]


def task_a(seed=20261020, trials=60):
    rng = random.Random(seed)
    checks = []
    for _ in range(trials):
        while True:
            p, q = rng.randint(1, 9), rng.randint(1, 9)
            m, n = rng.randint(-6, 6), rng.randint(-6, 6)
            if 3 * p * q > enorm((m, n)):
                break
        G0 = trace_gram(p, q, m, n)
        ldl(G0)
        require(multiply(multiply(transpose(R), G0), R) == G0, "R is not an isometry")
        require(determinant(G0) == (3 * p * q - enorm((m, n))) ** 2, "determinant formula")
        require(params_from_gram(G0) == (p, q, m, n), "parameter read-back")
        # a random Z[w]-basis change preserves the shape and the invariant 3pq - N(delta)
        while True:
            x, y, z, w = [(rng.randint(-2, 2), rng.randint(-2, 2)) for _ in range(4)]
            det = (emul(x, w)[0] - emul(y, z)[0], emul(x, w)[1] - emul(y, z)[1])
            if is_unit(det):
                break
        T = basis_change(x, y, z, w)
        require(abs(determinant(T)) == 1 and multiply(T, R) == multiply(R, T), "basis change")
        p2, q2, m2, n2 = params_from_gram(gram_of(G0, T))
        require(3 * p2 * q2 - enorm((m2, n2)) == 3 * p * q - enorm((m, n)), "invariant")
        row = {"p": p, "q": q, "m": m, "n": n, "det": determinant(G0), "basis_change": [x, y, z, w]}
        if p % 2 == 0 and m % 2 == 0 and n % 2 == 0:
            Dp = (3 * p * q - enorm((m, n))) // 2
            for v in (V1, V2, V3):
                require(determinant(overlattice_gram(G0, v)) == Dp * Dp, "overlattice determinant")
            row["line_integral"] = True
            row["D_prime"] = Dp
        else:
            # at least one of the three overlattices is not integral
            bad = 0
            for v in (V1, V2, V3):
                try:
                    overlattice_gram(G0, v)
                except InvalidCertificate:
                    bad += 1
            require(bad > 0, "parity condition is not necessary")
            row["line_integral"] = False
        checks.append(row)
    # the control of section 4
    control = control_check()
    return {"schema": SCHEMA + "-task-a", "random_checks": len(checks), "seed": seed,
            "all_passed": True, "control": control, "sample": checks[:5]}


def control_check():
    p, q, m, n = 8, 14, -4, 2
    G0 = trace_gram(p, q, m, n)
    require(G0 == [[16, -8, 2, 4], [-8, 16, -6, 2], [2, -6, 28, -14], [4, 2, -14, 28]], "control G0")
    require(multiply(multiply(transpose(R), G0), R) == G0, "control R")
    require(determinant(G0) == 94864 == (3 * p * q - enorm((m, n))) ** 2, "control determinant")
    O1, O2, O3 = (overlattice_gram(G0, v) for v in (V1, V2, V3))
    require(determinant(O1) == determinant(O2) == determinant(O3) == 23716, "control coset determinants")
    G = rank4_counterexample.G
    L1 = [[48, 12, -8, -2], [12, 4, -3, 1], [-8, -3, 28, -14], [-2, 1, -14, 28]]
    L2 = [[48, 36, -8, -2], [36, 28, -7, 0], [-8, -7, 28, -14], [-2, 0, -14, 28]]
    links = {"O3_to_G": len(isometries(O3, G, 5_000_000)), "O1_to_L1": len(isometries(O1, L1, 5_000_000)),
             "O2_to_L2": len(isometries(O2, L2, 5_000_000)), "Aut_O3": len(isometries(O3, O3, 5_000_000))}
    require(links == {"O3_to_G": 2, "O1_to_L1": 2, "O2_to_L2": 2, "Aut_O3": 2}, "control isometries")
    res = analyse_kernel(G0, prefix=60, max_nodes=5_000_000)
    require(res["aut_order"] == 6 and res["tally"] == {"level_2": 15, "nonzero": 90}, "control analyse_kernel")
    require(all(p["smallest_integral_scaling"] == 1 and p["det_scaled_L"] == 23716 for p in res["zero_pairs"]),
            "control scalings")
    counts = {"5..41": len(reduced_triples(5, 41)), "5..154": len(reduced_triples(5, 154))}
    require(counts == {"5..41": 498, "5..154": 5539}, "control triple counts")
    return {"G0": G0, "D_prime": 154, "links": links, "aut_L0_order": 6, "tally": res["tally"],
            "triple_counts": counts}


# ----------------------------------------------------------------------------
# Task B: complete enumeration of F
# ----------------------------------------------------------------------------

def reduced_triples(Dmin, Dmax):
    """All (p0, delta0, q0) with 3 p0 q0 - N(delta0) = 2D', p0 <= sqrt(D'), N(delta0) <= p0^2, q0 >= p0.

    Covering: every binary Hermitian Z[w]-lattice with integral trace form has a
    basis (f1, f2) with p0 = h(f1,f1) minimal, q0 >= p0 and delta0 reduced modulo
    p0 sqrt(-3) Z[w] (covering radius p0), so N(delta0) <= p0^2.  Duplicates are
    harmless; the list is a superset of the reduced forms.
    """
    out = []
    for Dp in range(Dmin, Dmax + 1):
        E = 2 * Dp
        for p0 in range(1, isqrt(Dp) + 1):
            for m0 in range(-2 * p0, 2 * p0 + 1):
                for n0 in range(-2 * p0, 2 * p0 + 1):
                    N = m0 * m0 - m0 * n0 + n0 * n0
                    if N > p0 * p0 or (E + N) % (3 * p0):
                        continue
                    q0 = (E + N) // (3 * p0)
                    if q0 >= p0:
                        out.append((Dp, p0, (m0, n0), q0))
    return out


def primitive_part(G):
    c = 0
    for row in G:
        for x in row:
            c = gcd(c, abs(x))
    return c, [[x // c for x in row] for row in G]


def analyse_member(Dp, p0, delta0, q0, line_index, prefix=60, max_nodes=2_000_000):
    """Transform the line to Z[w]e1', test the parity condition, classify."""
    x, y, z, w = LINE_BASIS_CHANGES[line_index]
    T = basis_change(x, y, z, w)
    G0_red = trace_gram(p0, q0, delta0[0], delta0[1])
    G0 = gram_of(G0_red, T)
    p, q, m, n = params_from_gram(G0)
    require(3 * p * q - enorm((m, n)) == 2 * Dp, "invariant lost in the basis change")
    params = line_parameters(p, q, m, n)
    record = {"D_prime": Dp, "reduced_triple": {"p0": p0, "delta0": list(delta0), "q0": q0},
              "line": LINE_NAMES[line_index], "basis_change": T, "p": p, "q": q, "delta": [m, n],
              "G0": G0, "member": params is not None}
    if params is None:
        record["excluded"] = "line fails the parity condition (p even, delta in 2Z[w])"
        return record
    require(params["D_prime"] == Dp, "D' formula")
    record.update(params)
    res = analyse_kernel(G0, prefix, max_nodes)
    group = isometries(G0, G0, max_nodes)
    require(len(group) == res["aut_order"], "group order mismatch between the two searches")
    perms = [[class_index(matvec(T2, class_vector(k))) for k in range(16)] for T2 in group]
    k1, k2, k3 = class_index(V1), class_index(V2), class_index(V3)
    line_pairs = {}
    for ka, kb in ((k1, k2), (k1, k3), (k2, k3)):
        pr = next(pz for pz in res["zero_pairs"] if class_index(pz["v1"]) == ka and class_index(pz["v2"]) == kb)
        line_pairs[f"{ka},{kb}"] = pr["level"]
    levels = set(line_pairs.values())
    require(len(levels) == 1 and levels <= {1, 2}, "line pairs must share a level in {1,2}")
    level = levels.pop()
    pr = next(pz for pz in res["zero_pairs"] if class_index(pz["v1"]) == k1 and class_index(pz["v2"]) == k2)
    require(pr["smallest_integral_scaling"] == 1, "member must have integral L")
    GL4 = pr["L_gram_times_4"]
    GL = [[x // 4 for x in row] for row in GL4]
    require(determinant(GL) == Dp * Dp, "det L")
    c, GLp = primitive_part(GL)
    swap = None
    kinds = Counter()
    for T2, perm in zip(group, perms):
        if perm[k1] == k2 and perm[k2] == k1:
            TR, RT = multiply(T2, R), multiply(R, T2)
            kind = ("conjugate-linear" if TR == multiply(multiply(R, R), T2) else
                    "linear" if TR == RT else "does not normalise R")
            kinds[kind] += 1
            if swap is None:
                swap = {"T": T2, "kind": kind}
    require((level == 1) == (swap is not None), "level/swap consistency")
    # the characteristic on L = O_v3 in the basis of L (columns of C/2, C = overlattice basis of v3)
    C = overlattice_basis(V3)
    half = [[Fraction(v, 2) for v in row] for row in C]
    b = [0 if all(x.denominator == 1 for x in [half[r][i] for r in range(4)]) else 1 for i in range(4)]
    a = matvec(inverse(half), V1)
    require(all(x.denominator == 1 for x in a), "a integral")
    a = [int(x) % 2 for x in a]
    require(any(a) and any(b) and dot(a, b) % 2 == 0, "characteristic")
    U, GLr = lll_reduce(GLp)
    ar, br = transport_characteristic(U, a, b)
    off_line = [pz for pz in res["zero_pairs"] if not line_of_same_orbit(pz, perms)]
    record.update({"level": level, "aut_L0_order": len(group), "orbits_on_classes": res["orbits_on_classes"],
                   "tally_105": res["tally"], "line_pair_levels": line_pairs, "swap": swap,
                   "swap_kinds": dict(kinds), "L_gram": GL, "scaling_c": c, "L_primitive_gram": GLp,
                   "primitive_determinant": determinant(GLp), "L_reduced": GLr, "L_reduced_transform": U,
                   "a": ar, "b": br, "a_unreduced": a, "b_unreduced": b,
                   "accidental_zero_pairs": [{"v1": pz["v1"], "v2": pz["v2"], "level": pz["level"]} for pz in off_line],
                   "level_3_or_4_pairs": [pz for pz in res["zero_pairs"] if pz["level"] >= 3]})
    return record


def line_of_same_orbit(pz, perms):
    """True when v2 is R v1 or R^2 v1, i.e. the pair lies on one of the five F4-lines."""
    k1, k2 = class_index(pz["v1"]), class_index(pz["v2"])
    rp = [class_index(matvec(R, class_vector(k))) for k in range(16)]
    return rp[k1] == k2 or rp[rp[k1]] == k2


def same_example(rec1, rec2, max_nodes=2_000_000):
    """Isometry of the primitive L's carrying (a,b) of rec1 to (a,b) of rec2, or None."""
    if rec1["primitive_determinant"] != rec2["primitive_determinant"]:
        return None
    maps = isometries(rec1["L_reduced"], rec2["L_reduced"], max_nodes)
    for U in maps:                       # U^t G2 U = G1: x1 -> U x1 is the isometry L1 -> L2
        a2 = [v % 2 for v in matvec(U, rec1["a"])]
        b2 = matvec(transpose(inverse(U)), rec1["b"])
        require(all(v.denominator == 1 for v in b2), "dual transport")
        b2 = [int(v) % 2 for v in b2]
        if a2 == rec2["a"] and b2 == rec2["b"]:
            return U
    return None


def task_b(Dmin=5, Dmax=41, prefix=60, max_nodes=2_000_000, progress=None):
    start = time.perf_counter()
    triples = reduced_triples(Dmin, Dmax)
    members, excluded = [], Counter()
    for i, (Dp, p0, delta0, q0) in enumerate(triples):
        for line in range(5):
            rec = analyse_member(Dp, p0, delta0, q0, line, prefix, max_nodes)
            if rec["member"]:
                members.append(rec)
            else:
                excluded["parity"] += 1
        if progress and i % progress == 0:
            print(f"  triples {i}/{len(triples)} members {len(members)} {time.perf_counter() - start:.0f}s", flush=True)
    # deduplicate by exact isometry of (primitive L, characteristic)
    representatives = []
    for idx, rec in enumerate(members):
        match = None
        for j, rep in enumerate(representatives):
            U = same_example(members[rep["member_index"]], rec, max_nodes)
            if U is not None:
                match = (j, U)
                break
        if match is None:
            rec["example_index"] = len(representatives)
            rec["isometry_to_representative"] = identity(4)
            representatives.append({"member_index": idx, "D_prime": rec["D_prime"], "level": rec["level"],
                                    "primitive_determinant": rec["primitive_determinant"], "scaling_c": rec["scaling_c"]})
        else:
            rec["example_index"] = match[0]
            rec["isometry_to_representative"] = match[1]
    per_D = {}
    for Dp in range(Dmin, Dmax + 1):
        trip = sum(1 for t in triples if t[0] == Dp)
        mem = [m for m in members if m["D_prime"] == Dp]
        reps = [r for r in representatives if r["D_prime"] == Dp]
        l2 = [r for r in reps if r["level"] == 2]
        per_D[Dp] = {"reduced_triples": trip, "line_candidates": 5 * trip, "members": len(mem),
                     "distinct_examples": len(reps), "level_1": sum(r["level"] == 1 for r in reps),
                     "level_2": len(l2), "non_primitive_members": sum(m["scaling_c"] > 1 for m in mem),
                     "smallest_level_2_primitive_determinant": min((r["primitive_determinant"] for r in l2), default=None)}
    level2 = [r for r in representatives if r["level"] == 2]
    doc = {"schema": SCHEMA + "-task-b", "universe": {
               "description": "family F: binary Hermitian Z[w]-lattices with integral trace form, all five "
                              "F4-lines, members = lines with integral overlattices (p even, delta in 2Z[w])",
               "D_prime_range": [Dmin, Dmax], "reduced_triples": len(triples), "line_candidates": 5 * len(triples),
               "reduction": "p0 = min h <= sqrt(D'), q0 >= p0, N(delta0) <= p0^2 (box |m0|,|n0| <= 2 p0)",
               "dedup": "exact isometry of the primitive L carrying (a,b) mod 2 to (a',b') mod 2 "
                        "(rank4_feasibility.isometries, complete)",
               "isospectrality_filter": f"coset theta prefixes through norm {prefix} (analyse_kernel)",
               "level_test": "swap in the complete Aut(L0) = level 1, else level 2 (R maps v1 to v2)"},
           "excluded": {"lines_failing_parity": excluded["parity"]},
           "members": members, "representatives": representatives, "per_D_prime": per_D,
           "level_2_examples": level2,
           "accidental_pairs": [{"member": i, "pairs": m["accidental_zero_pairs"]} for i, m in enumerate(members) if m["accidental_zero_pairs"]],
           "level_3_or_4_candidates": [{"member": i, "pairs": m["level_3_or_4_pairs"]} for i, m in enumerate(members) if m["level_3_or_4_pairs"]],
           "seconds": time.perf_counter() - start}
    doc["summary"] = {"reduced_triples": len(triples), "members": len(members), "distinct_examples": len(representatives),
                      "level_1_examples": sum(r["level"] == 1 for r in representatives), "level_2_examples": len(level2),
                      "smallest_level_2_primitive_determinant": min((r["primitive_determinant"] for r in level2), default=None),
                      "accidental_zero_pairs": sum(len(m["accidental_zero_pairs"]) for m in members),
                      "level_3_or_4_candidates": sum(len(m["level_3_or_4_pairs"]) for m in members),
                      "seconds": round(doc["seconds"], 1)}
    return doc


# ----------------------------------------------------------------------------
# Task C: certification of a level-2 member
# ----------------------------------------------------------------------------

def decomposition_test(G, group):
    """L is decomposable iff some involution T != +-I in Aut(L) has (I+T)/2 integral.

    If L = A (+) B then T = (-1 on A, +1 on B) is such an involution; conversely
    (I+T)/2 integral is the orthogonal projection onto the fixed sublattice and
    splits L = L^T (+) L^{-T}.  Returns the splitting involution or None.
    """
    I = identity(4)
    for T in group:
        if T == I or T == [[-x for x in row] for row in I]:
            continue
        if multiply(T, T) != I:
            continue
        if all((x + y) % 2 == 0 for rT, rI in zip(T, I) for x, y in zip(rT, rI)):
            return T
    return None


def manuscript_schema_certificate(G, a, b, K, G0, t, group, max_nodes=2_000_000):
    """Route (iii): the ambient-coset certificate of rank4_counterexample, when G0/2 is integral."""
    if any(x % 2 for row in G0 for x in row):
        return {"available": False, "reason": "G0/2 is not integral (q odd); the manuscript schema needs an even ambient lattice"}
    M = [[x // 2 for x in row] for row in G0]
    Kinv = inverse(K)
    V = [[2 * x for x in row] for row in Kinv]
    require(all(Fraction(x).denominator == 1 for row in V for x in row), "embedding")
    V = [[int(x) for x in row] for row in V]
    plus = [x % 4 for x in matvec(V, a)]
    minus = [x % 4 for x in matvec(V, [p + 2 * q for p, q in zip(a, t)])]
    amb = None
    for cand in (R, multiply(R, R)):
        if all((x - y) % 4 == 0 for x, y in zip(matvec(cand, minus), plus)):
            amb = cand
            break
    require(amb is not None, "neither R nor R^2 carries the negative coset onto the positive one")
    doc = {"schema": "work7-rank4-ambient-coset-counterexample-v1", "G": G, "a": a, "b": b,
           "ambient_metric": M, "embedding": V, "ambient_isometry": amb, "index_two_basis": K,
           "common_coset_basis_change": identity(4), "complementary_vector": t,
           "positive_coset_representative": plus, "negative_coset_representative": minus,
           "automorphisms": group, "corroboration": complete_shells(G, a, b, 64)}
    replay = rank4_counterexample.replay(doc, max_nodes)
    require(replay["verdict"] == "verified_rank4_symmetry_converse_counterexample", "manuscript replay failed")
    return {"available": True, "document": doc, "replay_verdict": replay["verdict"],
            "replay_full_group_order": replay["full_group_order"]}


def certify_member(rec, cutoff_vector_limit=30_000_000, max_nodes=2_000_000, enumeration_nodes=2_000_000_000):
    """Section-5 certificate for a level-2 member record of task_b."""
    start = time.perf_counter()
    require(rec["member"] and rec["level"] == 2, "only level-2 members are certified")
    G0_std = rec["G0"]
    c = rec["scaling_c"]
    require(all(x % c == 0 for row in G0_std for x in row), "scaling must divide the kernel Gram")
    G0 = [[x // c for x in row] for row in G0_std]          # kernel of the primitive L
    # L = O_v3 in L0-coordinates: basis C/2 with C = overlattice basis of v3; Gram GL = C^t G0_std C / 4
    C = overlattice_basis(V3)
    GL = gram_of(G0_std, C)
    require(all(x % (4 * c) == 0 for row in GL for x in row), "L Gram divisibility")
    GL = [[x // (4 * c) for x in row] for row in GL]
    require(GL == rec["L_primitive_gram"], "primitive L mismatch")
    half = [[Fraction(v, 2) for v in row] for row in C]
    Kfull = inverse(half)                                     # L0-standard basis vectors in L-coordinates
    require(all(x.denominator == 1 for row in Kfull for x in row), "kernel basis integral")
    K_un = [[int(x) for x in row] for row in Kfull]
    require(gram_of(GL, K_un) == G0, "kernel Gram")
    U, G = lll_reduce(GL)
    require(G == rec["L_reduced"], "reduced L mismatch")
    a, b = transport_characteristic(U, rec["a_unreduced"], rec["b_unreduced"])
    require([a, b] == [rec["a"], rec["b"]], "characteristic mismatch")
    Uinv = inverse(U)
    K = multiply(Uinv, K_un)
    require(all(Fraction(x).denominator == 1 for row in K for x in row), "kernel transport")
    K = [[int(x) for x in row] for row in K]
    require(gram_of(G, K) == G0 and abs(determinant(K)) == 2, "kernel basis in reduced coordinates")
    require(all(dot(col, b) % 2 == 0 for col in transpose(K)), "kernel basis in ker b")
    det = determinant(G)
    require(determinant(G0) == 4 * det, "d0")
    _, pivots = ldl(G)
    j = next(i for i in range(4) if b[i])
    t = [int(i == j) for i in range(4)]
    cols = transpose(K)
    C1 = hnf_basis([[2 * v for v in col] for col in cols] + [a])
    C2 = hnf_basis([[2 * v for v in col] for col in cols] + [[p + 2 * q for p, q in zip(a, t)]])
    O = []
    for Cx in (C1, C2):
        g = gram_of(G, Cx)
        require(all(x % 4 == 0 for row in g for x in row), "coset lattice not integral")
        O.append([[x // 4 for x in row] for row in g])
        require(determinant(O[-1]) == det, "coset determinant")
    v1 = matvec(inverse(K), a)
    v2 = matvec(inverse(K), [p + 2 * q for p, q in zip(a, t)])
    require(all(x.denominator == 1 for x in v1 + v2), "classes in the kernel")
    v1, v2 = [int(x) % 2 for x in v1], [int(x) % 2 for x in v2]
    # route (i): the kernel rotation
    R0 = R
    require(multiply(multiply(transpose(R0), G0), R0) == G0 and multiply(R0, multiply(R0, R0)) == identity(4), "R0")
    require([x % 2 for x in matvec(R0, v1)] == v2 or [x % 2 for x in matvec(R0, v2)] == v1, "R0 does not relate the classes")
    rotation = {"R0": R0, "maps_v1_to_v2": [x % 2 for x in matvec(R0, v1)] == v2,
                "maps_v2_to_v1": [x % 2 for x in matvec(R0, v2)] == v1}
    # route (ii): modular cutoff
    item = modular_data(G, a, b)
    B = item["cutoff"]
    estimate = (3.141592653589793 ** 2 / 2) * B * B / (16 * det ** 0.5)
    cutoff = {"cutoff": B, "modular": item, "estimated_vectors": estimate, "attempted": False}
    if estimate <= cutoff_vector_limit:
        t0 = time.perf_counter()
        ldl_doc = cutoff_build(G, a, b, max_nodes=enumeration_nodes)
        t1 = time.perf_counter()
        mine = signed_shells(G, a, b, B, enumeration_nodes)
        t2 = time.perf_counter()
        require(ldl_doc["evidence"]["enumeration_complete"] is True, "LDL cutoff interrupted")
        require(ldl_doc["evidence"]["enumeration"]["shells"] == mine["shells"], "backends disagree")
        require(ldl_doc["evidence"]["verdict"] == "proved_zero_by_modular_cutoff", "cutoff did not prove the zero")
        cutoff.update({"attempted": True, "ldl_backend": ldl_doc,
                       "integer_backend": {"algorithm": mine["algorithm"], "visited_nodes": mine["visited_nodes"],
                                           "shells_sha256": hashlib.sha256(json.dumps(mine["shells"]).encode()).hexdigest(),
                                           "vectors": sum(s["vector_count"] for s in mine["shells"])},
                       "seconds": {"ldl": t1 - t0, "integer": t2 - t1}})
    else:
        cutoff["reason"] = f"estimated {estimate:.3g} vectors exceeds the limit {cutoff_vector_limit}"
    # full group of L, stabilizer signs, indecomposability
    group = isometries(G, G, max_nodes)
    stabilizer = []
    for T in group:
        if any((x - y) % 2 for x, y in zip(matvec(T, a), a)):
            continue
        dual = matvec(transpose(inverse(T)), b)
        if any((x - y) % 2 for x, y in zip(dual, b)):
            continue
        stabilizer.append({"T": T, "epsilon": symmetry_data(G, a, b, T)["epsilon"]})
    require(stabilizer and all(s["epsilon"] == 0 for s in stabilizer), "a symmetry witness exists")
    split = decomposition_test(G, group)
    require(split is None, "L is decomposable")
    # route (iii)
    manuscript = manuscript_schema_certificate(G, a, b, K, G0, t, group, max_nodes)
    # Aut(L0) and its action
    group0 = isometries(G0, G0, max_nodes)
    perms = [[class_index(matvec(T, class_vector(k))) for k in range(16)] for T in group0]
    k1, k2 = class_index(v1), class_index(v2)
    mapping = sum(1 for p in perms if p[k1] == k2)
    swapping = sum(1 for p in perms if p[k1] == k2 and p[k2] == k1)
    require(mapping >= 1 and swapping == 0, "level is not 2")
    orbits, _ = orbits_on_classes(group0)
    return {"schema": SCHEMA + "-level2-certificate", "member": {"a": rec["a_param"], "q": rec["q"], "beta": rec["beta"],
            "p": rec["p"], "delta": rec["delta"], "D_prime": rec["D_prime"], "line": rec["line"],
            "reduced_triple": rec["reduced_triple"]}, "G0_standard": G0_std, "scaling_c": c,
            "G": G, "determinant": det, "primitive": True, "ldl_pivots": [str(x) for x in pivots],
            "lll_transform": U, "L_gram_in_kernel_coordinates": GL, "a": a, "b": b,
            "kernel_basis": K, "G0": G0, "d0": determinant(G0), "complement_t": t,
            "coset_lattices": {"O_v1": {"doubled_basis": C1, "gram": O[0]}, "O_v2": {"doubled_basis": C2, "gram": O[1]}},
            "classes": {"v1": v1, "v2": v2},
            "identity": {"kernel_rotation": rotation, "modular_cutoff": cutoff, "manuscript_schema": manuscript},
            "automorphisms": group, "stabilizer": stabilizer, "indecomposable": True,
            "indecomposability_argument": ("Aut(L) = {I,-I}: a decomposable lattice has at least four automorphisms"
                                           if len(group) == 2 else "no involution of Aut(L) other than +-I has (I+T)/2 integral"),
            "kernel_automorphisms": {"order": len(group0), "orbits_on_classes": orbits,
                                     "elements_mapping_v1_to_v2": mapping, "elements_swapping_v1_v2": swapping},
            "level": 2, "seconds": time.perf_counter() - start}


def task_c(task_b_doc, limit_examples=5, cutoff_vector_limit=30_000_000, max_nodes=2_000_000, det_bound=1729):
    reps = sorted((r for r in task_b_doc["representatives"] if r["level"] == 2),
                  key=lambda r: (r["primitive_determinant"], r["D_prime"]))
    members = task_b_doc["members"]
    chosen, seen_dets = [], []
    for r in reps:
        d = r["primitive_determinant"]
        if d >= det_bound:
            continue
        if d not in seen_dets:
            seen_dets.append(d)
        if len(seen_dets) <= limit_examples or d in seen_dets[:limit_examples]:
            chosen.append(r)
    if not chosen and reps:
        chosen = [reps[0]]
    certificates = [certify_member(members[r["member_index"]], cutoff_vector_limit, max_nodes) for r in chosen]
    return {"schema": SCHEMA + "-certificates", "det_bound": det_bound, "certificates": certificates,
            "summary_table": [{"D_prime": r["D_prime"], "primitive_determinant": r["primitive_determinant"],
                               "member_index": r["member_index"]} for r in reps]}


# ----------------------------------------------------------------------------
# Task D: extension and structural tally
# ----------------------------------------------------------------------------

def task_d(Dmin=5, Dmax=154, prefix=60, max_nodes=2_000_000, progress=None, task_b_doc=None):
    doc = task_b(Dmin, Dmax, prefix, max_nodes, progress)
    doc["schema"] = SCHEMA + "-task-d"
    kinds = Counter()
    for m in doc["members"]:
        if m["level"] == 1:
            kinds.update(m["swap_kinds"])
    doc["structural_tally"] = {"swapping_elements_by_kind_over_all_level_1_members": dict(kinds),
                               "level_1_members_with_conjugate_linear_swap": sum(1 for m in doc["members"] if m["level"] == 1 and m["swap_kinds"].get("conjugate-linear")),
                               "level_1_members_with_only_non_normalising_swaps": sum(1 for m in doc["members"] if m["level"] == 1 and set(m["swap_kinds"]) == {"does not normalise R"})}
    # recover the det-23,716 example
    target_G, target_a, target_b = rank4_counterexample.G, rank4_counterexample.A, rank4_counterexample.B
    Ut, Gt = lll_reduce(target_G)
    at, bt = transport_characteristic(Ut, target_a, target_b)
    target = {"L_reduced": Gt, "a": at, "b": bt, "primitive_determinant": 23716}
    hits = []
    for i, m in enumerate(doc["members"]):
        if m["D_prime"] == 154 and m["primitive_determinant"] == 23716:
            U = same_example(m, target, max_nodes)
            if U is not None:
                hits.append({"member": i, "a_param": m["a_param"], "q": m["q"], "beta": m["beta"], "isometry": U})
    doc["det_23716_recovered"] = {"members_matching": len(hits), "hits": hits[:3]}
    return doc


# ----------------------------------------------------------------------------
# Task E: the mixed-type lemma
# ----------------------------------------------------------------------------

def mixed_lattice(N, c, glue_N):
    """Glue N (rank two, R-fixed) with c A2 along a group of order three.

    Returns the Gram of L0 in an integral basis, the integer matrix of R, the
    three classes of the unique F4-line, and the candidate automorphisms
    sigma_k = (-1 on N) (+) (w^k kappa on c A2) in the same basis.
    """
    A2 = [[2 * c, -c], [-c, 2 * c]]
    Gdiag = [N[0] + [0, 0], N[1] + [0, 0], [0, 0] + A2[0], [0, 0] + A2[1]]
    glue = [Fraction(x) for x in glue_N] + [Fraction(1, 3), Fraction(2, 3)]     # (1/sqrt(-3)) r1 = -(r1 + 2 w r1)/3
    generators = [[3 * int(i == j) for i in range(4)] for j in range(4)] + [[int(3 * x) for x in glue]]
    B3 = hnf_basis(generators)                                   # basis of 3 L0 in the Q-basis
    G9 = gram_of(Gdiag, B3)
    require(all(x % 9 == 0 for row in G9 for x in row), "glue is not integral")
    G = [[x // 9 for x in row] for row in G9]
    rho = mult_matrix((0, 1))
    RQ = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0] + rho[0], [0, 0] + rho[1]]
    B3inv = inverse(B3)
    Rnew = multiply(multiply(B3inv, RQ), B3)
    require(all(x.denominator == 1 for row in Rnew for x in row), "glue is not R-stable")
    Rnew = [[int(x) for x in row] for row in Rnew]
    require(multiply(multiply(transpose(Rnew), G), Rnew) == G, "R is not an isometry of the glued lattice")
    coords = lambda v: [int(3 * x) for x in matvec(B3inv, v)]   # coordinates in the basis B3/3
    line = [[x % 2 for x in coords([0, 0, 1, 0])], [x % 2 for x in coords([0, 0, 0, 1])], [x % 2 for x in coords([0, 0, 1, 1])]]
    kappa = [[1, -1], [0, -1]]
    sigmas = []
    for k in range(3):
        rk = identity(2)
        for _ in range(k):
            rk = multiply(rk, rho)
        blk = multiply(rk, kappa)
        SQ = [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0] + blk[0], [0, 0] + blk[1]]
        S = multiply(multiply(B3inv, SQ), B3)
        if all(x.denominator == 1 for row in S for x in row):
            S = [[int(x) for x in row] for row in S]
            sigmas.append({"k": k, "T": S, "is_isometry": multiply(multiply(transpose(S), G), S) == G})
        else:
            sigmas.append({"k": k, "T": None, "is_isometry": False, "reason": "not integral on L0"})
    return {"N": N, "c": c, "glue_N": [str(x) for x in glue_N], "G0": G, "determinant": determinant(G),
            "R": Rnew, "line_classes": line, "sigmas": sigmas}


def task_e(prefix=60, max_nodes=2_000_000):
    cases = [([[1, 0], [0, 3]], 1, [0, Fraction(1, 3)]), ([[2, 1], [1, 2]], 2, [Fraction(1, 3), Fraction(1, 3)]),
             ([[4, -2], [-2, 4]], 1, [Fraction(1, 3), Fraction(2, 3)]), ([[3, 0], [0, 3]], 1, [Fraction(1, 3), 0]),
             ([[1, 0], [0, 3]], 4, [0, Fraction(1, 3)]), ([[2, 1], [1, 2]], 5, [Fraction(1, 3), Fraction(1, 3)]),
             ([[3, 0], [0, 3]], 4, [Fraction(1, 3), 0]), ([[1, 0], [0, 3]], 7, [0, Fraction(1, 3)])]
    rows = []
    for N, c, glue in cases:
        mixed = mixed_lattice(N, c, glue)
        G0 = mixed["G0"]
        res = analyse_kernel(G0, prefix, max_nodes)
        group = isometries(G0, G0, max_nodes)
        perms = [[class_index(matvec(T, class_vector(k))) for k in range(16)] for T in group]
        k1, k2, k3 = (class_index(v) for v in mixed["line_classes"])
        require(k1 ^ k2 == k3 and len({k1, k2, k3}) == 3, "line classes")
        pr = next(p for p in res["zero_pairs"] if {class_index(p["v1"]), class_index(p["v2"])} == {k1, k2})
        swap = next((T for T, p in zip(group, perms) if p[k1] == k2 and p[k2] == k1), None)
        sigma_swaps = []
        for s in mixed["sigmas"]:
            if s["T"] is not None and s["is_isometry"]:
                p = [class_index(matvec(s["T"], class_vector(k))) for k in range(16)]
                sigma_swaps.append({"k": s["k"], "swaps_v1_v2": p[k1] == k2 and p[k2] == k1, "fixes_v3": p[k3] == k3})
        rows.append({"N": N, "c": c, "glue_N": mixed["glue_N"], "G0": G0, "determinant": mixed["determinant"],
                     "R": mixed["R"], "line_classes": mixed["line_classes"], "aut_L0_order": len(group),
                     "line_level": pr["level"], "swapping_element": swap, "sigma_candidates": mixed["sigmas"],
                     "sigma_swaps": sigma_swaps, "tally_105": res["tally"]})
        require(pr["level"] == 1 and any(s["swaps_v1_v2"] for s in sigma_swaps), "lemma fails for this mixed lattice")
    return {"schema": SCHEMA + "-task-e", "cases": rows, "lemma_holds_in_all_cases": True,
            "statement": "every mixed-type order-three Klein configuration tested is level 1; the swapping "
                         "element is (-1 on N) (+) (w^k kappa on c A2) for some k"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", choices=("A", "B", "C", "D", "E"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--task-b", type=Path, help="Task B result file (for C)")
    parser.add_argument("--dmin", type=int, default=5)
    parser.add_argument("--dmax", type=int, default=41)
    parser.add_argument("--prefix", type=int, default=60)
    parser.add_argument("--max-nodes", type=int, default=2_000_000)
    parser.add_argument("--cutoff-limit", type=float, default=30_000_000)
    parser.add_argument("--limit-examples", type=int, default=5)
    parser.add_argument("--det-bound", type=int, default=1729)
    parser.add_argument("--progress", type=int, default=0)
    args = parser.parse_args()
    require(not args.output.exists(), "choose a fresh output path")
    if args.task == "A":
        doc = task_a()
        summary = {"random_checks": doc["random_checks"], "control": {k: v for k, v in doc["control"].items() if k != "G0"}}
    elif args.task == "B":
        doc = task_b(args.dmin, args.dmax, args.prefix, args.max_nodes, args.progress or None)
        summary = doc["summary"]
    elif args.task == "C":
        doc = task_c(json.loads(args.task_b.read_text(encoding="utf-8")), args.limit_examples, args.cutoff_limit,
                     args.max_nodes, args.det_bound)
        summary = [{"D_prime": c["member"]["D_prime"], "determinant": c["determinant"], "G": c["G"], "a": c["a"], "b": c["b"],
                    "aut_order": len(c["automorphisms"]), "aut_L0_order": c["kernel_automorphisms"]["order"],
                    "cutoff": c["identity"]["modular_cutoff"]["cutoff"], "cutoff_attempted": c["identity"]["modular_cutoff"]["attempted"],
                    "manuscript_schema": c["identity"]["manuscript_schema"]["available"], "seconds": round(c["seconds"], 1)}
                   for c in doc["certificates"]]
    elif args.task == "D":
        doc = task_d(args.dmin, max(args.dmax, 154) if args.dmax == 41 else args.dmax, args.prefix, args.max_nodes,
                     args.progress or None)
        summary = {**doc["summary"], "structural_tally": doc["structural_tally"], "det_23716_recovered": doc["det_23716_recovered"]["members_matching"]}
    else:
        doc = task_e(args.prefix, args.max_nodes)
        summary = [{"N": r["N"], "c": r["c"], "det": r["determinant"], "aut": r["aut_L0_order"], "line_level": r["line_level"],
                    "sigma_swaps": r["sigma_swaps"]} for r in doc["cases"]]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
