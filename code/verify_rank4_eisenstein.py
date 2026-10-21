"""Independent replay of the Eisenstein (item 3) search: box backends only.

Replays results/rank4-eisenstein/task-b*.json and task-d*.json (the family
enumeration), examples/rank4-eisenstein/*.json (level-2 certificates) and
results/rank4-eisenstein/task-e.json (the mixed-type lemma).

Independence from rank4_eisenstein_search.py: the Eisenstein arithmetic,
trace form, basis changes, kernel and coset lattices are re-derived here from
the stored integers; vectors are enumerated with the box backends of
verify_rank4_level4 (exact Cauchy-Schwarz radii, no LDL recursion);
isometries and groups come from verify_rank4_feasibility.box_isometries; the
cutoff metadata from verify_rank4_feasibility.modular_check; the manuscript-
schema route from rank4_counterexample.replay unchanged.  Nothing stored is
trusted: every field is recomputed and compared.  Node limits raise
EnumerationLimit (verdict "unresolved"); any mismatch raises InvalidCertificate.
"""
import argparse
from collections import Counter
from fractions import Fraction
import hashlib
import json
from math import gcd, isqrt
from pathlib import Path
import time

from exact_theta import (EnumerationLimit, InvalidCertificate, dot, identity, integer, inverse, ldl,
                         matrix, matvec, multiply, require, symmetry_data, transpose, vector)
from verify_rank4_feasibility import box_isometries, det4, modular_check
from verify_rank4_level4 import box_class_histogram, box_signed_shells
import rank4_counterexample

R = [[0, -1, 0, 0], [1, -1, 0, 0], [0, 0, 0, -1], [0, 0, 1, -1]]
V1, V2, V3 = [1, 0, 0, 0], [0, 1, 0, 0], [1, 1, 0, 0]
LINES = [((1, 0), (0, 0), (0, 0), (1, 0)), ((1, 0), (1, 0), (0, 0), (1, 0)), ((1, 0), (0, 1), (0, 0), (1, 0)),
         ((1, 0), (-1, -1), (0, 0), (1, 0)), ((0, 0), (1, 0), (1, 0), (0, 0))]
LINE_NAMES = ["(1,0)", "(1,1)", "(1,w)", "(1,w^2)", "(0,1)"]


def class_index(v):
    return sum((x & 1) << i for i, x in enumerate(v))


def class_vector(k):
    return [(k >> i) & 1 for i in range(4)]


def gram_of(G, B):
    return multiply(multiply(transpose(B), G), B)


def enorm(x):
    m, n = x
    return m * m - m * n + n * n


def mult_matrix(x):
    m, n = x
    return [[m, -n], [n, m - n]]


def trace_gram(p, q, m, n):
    return [[2 * p, -p, n, -m], [-p, 2 * p, m - n, n], [n, m - n, 2 * q, -q], [-m, n, -q, 2 * q]]


def basis_change(x, y, z, w):
    X, Y, Z, W = mult_matrix(x), mult_matrix(y), mult_matrix(z), mult_matrix(w)
    return [X[0] + Z[0], X[1] + Z[1], Y[0] + W[0], Y[1] + W[1]]


def params(G):
    matrix(G, "G0", 4)
    require(G[0][0] % 2 == 0 and G[2][2] % 2 == 0, "odd diagonal")
    p, q, n, m = G[0][0] // 2, G[2][2] // 2, G[0][2], -G[0][3]
    require(trace_gram(p, q, m, n) == G, "not a trace form")
    return p, q, m, n


def over_basis(v):
    j = next(i for i in range(4) if v[i])
    return transpose([list(v) if i == j else [2 if r == i else 0 for r in range(4)] for i in range(4)])


def reduced_triples(Dmin, Dmax):
    out = []
    for Dp in range(Dmin, Dmax + 1):
        for p0 in range(1, isqrt(Dp) + 1):
            for m0 in range(-2 * p0, 2 * p0 + 1):
                for n0 in range(-2 * p0, 2 * p0 + 1):
                    N = m0 * m0 - m0 * n0 + n0 * n0
                    if N > p0 * p0 or (2 * Dp + N) % (3 * p0):
                        continue
                    q0 = (2 * Dp + N) // (3 * p0)
                    if q0 >= p0:
                        out.append((Dp, p0, (m0, n0), q0))
    return out


def gcd_of(G):
    c = 0
    for row in G:
        for x in row:
            c = gcd(c, abs(x))
    return c


def unimodular(U):
    matrix(U, "U", 4)
    require(abs(det4(U)) == 1, "not unimodular")
    return U


def transport(U, a, b):
    """Characteristic in the basis x = U x'."""
    a2 = matvec(inverse(U), a)
    require(all(v.denominator == 1 for v in a2), "transport")
    return [int(v) % 2 for v in a2], [v % 2 for v in matvec(transpose(U), b)]


def push(U, a, b):
    """Characteristic pushed forward along x -> U x."""
    b2 = matvec(transpose(inverse(U)), b)
    require(all(v.denominator == 1 for v in b2), "dual push")
    return [v % 2 for v in matvec(U, a)], [int(v) % 2 for v in b2]


def line_characteristic(G0):
    """(a, b) of L = O_v3 in the basis C/2, C = over_basis(V3), plus the Gram times four."""
    C = over_basis(V3)
    half = [[Fraction(v, 2) for v in row] for row in C]
    b = [0 if all(half[r][i].denominator == 1 for r in range(4)) else 1 for i in range(4)]
    a = matvec(inverse(half), V1)
    require(all(x.denominator == 1 for x in a), "a")
    return [int(x) % 2 for x in a], b, gram_of(G0, C)


def classify_pairs(G0, prefix, max_nodes, enumeration_nodes):
    """Levels of all 105 pairs over L0 by box backends; returns (levels dict, group, perms, tally)."""
    hist, _ = box_class_histogram(G0, prefix, enumeration_nodes)
    group = box_isometries(G0, G0, max_nodes)
    perms = [[class_index(matvec(T, class_vector(k))) for k in range(16)] for T in group]
    levels, tally = {}, Counter()
    for k1 in range(1, 16):
        for k2 in range(k1 + 1, 16):
            if hist[k1] != hist[k2]:
                tally["nonzero"] += 1
                continue
            swap = any(p[k1] == k2 and p[k2] == k1 for p in perms)
            maps = any(p[k1] == k2 for p in perms)
            if swap:
                level = 1
            elif maps:
                level = 2
            else:
                A = gram_of(G0, over_basis(class_vector(k1)))
                B = gram_of(G0, over_basis(class_vector(k2)))
                level = 3 if box_isometries(A, B, max_nodes) else 4
            levels[(k1, k2)] = level
            tally[f"level_{level}"] += 1
    return levels, group, perms, dict(tally)


def on_r_line(k1, k2):
    rp = [class_index(matvec(R, class_vector(k))) for k in range(16)]
    return rp[k1] == k2 or rp[rp[k1]] == k2


def replay_member(m, prefix, max_nodes, enumeration_nodes):
    Dp = integer(m["D_prime"], "D'")
    t = m["reduced_triple"]
    p0, q0, d0 = integer(t["p0"], "p0"), integer(t["q0"], "q0"), (integer(t["delta0"][0], "m0"), integer(t["delta0"][1], "n0"))
    require(3 * p0 * q0 - enorm(d0) == 2 * Dp and p0 <= isqrt(Dp) and q0 >= p0 and enorm(d0) <= p0 * p0, "triple outside the universe")
    li = LINE_NAMES.index(m["line"])
    T = basis_change(*LINES[li])
    require(m["basis_change"] == T and abs(det4(T)) == 1 and multiply(T, R) == multiply(R, T), "basis change")
    G0 = gram_of(trace_gram(p0, q0, d0[0], d0[1]), T)
    require(m["G0"] == G0, "G0")
    p, q, mm, nn = params(G0)
    require([m["p"], m["q"], m["delta"]] == [p, q, [mm, nn]], "parameters")
    member = p % 2 == 0 and mm % 2 == 0 and nn % 2 == 0
    require(m["member"] == member, "parity condition")
    if not member:
        return None
    require(m["a_param"] == p // 2 and m["beta"] == [mm // 2, nn // 2] and 3 * (p // 2) * q - 2 * enorm((mm // 2, nn // 2)) == Dp, "line parameters")
    levels, group, perms, tally = classify_pairs(G0, prefix, max_nodes, enumeration_nodes)
    require(m["aut_L0_order"] == len(group) and m["tally_105"] == tally, "group order or tally")
    k1, k2, k3 = class_index(V1), class_index(V2), class_index(V3)
    for ka, kb in ((k1, k2), (k1, k3), (k2, k3)):
        require(m["line_pair_levels"][f"{ka},{kb}"] == levels[(ka, kb)] == m["level"], "line level")
    swaps = [(T2, pm) for T2, pm in zip(group, perms) if pm[k1] == k2 and pm[k2] == k1]
    require((m["level"] == 1) == bool(swaps), "level/swap")
    kinds = Counter()
    for T2, _ in swaps:
        TR = multiply(T2, R)
        kinds["conjugate-linear" if TR == multiply(multiply(R, R), T2) else "linear" if TR == multiply(R, T2) else "does not normalise R"] += 1
    require(m["swap_kinds"] == dict(kinds), "swap kinds")
    if swaps:
        S = matrix(m["swap"]["T"], "swap", 4)
        require(gram_of(G0, S) == G0 and any(S == T2 for T2, _ in swaps), "stored swapping element")
    else:
        require(m["swap"] is None, "swap stored at level 2")
    a_un, b_un, GL4 = line_characteristic(G0)
    require(all(x % 4 == 0 for row in GL4 for x in row), "L integral")
    GL = [[x // 4 for x in row] for row in GL4]
    require(m["L_gram"] == GL and det4(GL) == Dp * Dp, "Gram of L")
    c = gcd_of(GL)
    GLp = [[x // c for x in row] for row in GL]
    require(m["scaling_c"] == c and m["L_primitive_gram"] == GLp and m["primitive_determinant"] == det4(GLp), "primitive part")
    U = unimodular(m["L_reduced_transform"])
    require(gram_of(GLp, U) == m["L_reduced"], "reduced L")
    require([m["a_unreduced"], m["b_unreduced"]] == [a_un, b_un], "unreduced characteristic")
    require(list(transport(U, a_un, b_un)) == [m["a"], m["b"]], "transported characteristic")
    require(any(m["a"]) and any(m["b"]) and dot(m["a"], m["b"]) % 2 == 0, "characteristic parity")
    acc = [{"v1": class_vector(k1_), "v2": class_vector(k2_), "level": lv} for (k1_, k2_), lv in sorted(levels.items())
           if not on_r_line(k1_, k2_)]
    require(m["accidental_zero_pairs"] == acc, "accidental pairs")
    require([(p_["v1"], p_["v2"], p_["level"]) for p_ in m["level_3_or_4_pairs"]]
            == [(class_vector(k1_), class_vector(k2_), lv) for (k1_, k2_), lv in sorted(levels.items()) if lv >= 3], "level 3/4 list")
    return {"level": m["level"], "c": c}


def replay_enumeration(doc, max_nodes, enumeration_nodes):
    require(doc["schema"] in ("work7-rank4-eisenstein-v1-task-b", "work7-rank4-eisenstein-v1-task-d"), "schema")
    Dmin, Dmax = (integer(x, "range") for x in doc["universe"]["D_prime_range"])
    triples = reduced_triples(Dmin, Dmax)
    require(doc["universe"]["reduced_triples"] == len(triples) and doc["universe"]["line_candidates"] == 5 * len(triples), "universe size")
    prefix = int(doc["universe"]["isospectrality_filter"].split("norm")[1].split()[0])
    members = doc["members"]
    # every (triple, line) passing the parity test must be a member, in order
    expected = []
    for (Dp, p0, d0, q0) in triples:
        for li in range(5):
            G0 = gram_of(trace_gram(p0, q0, d0[0], d0[1]), basis_change(*LINES[li]))
            p, q, mm, nn = params(G0)
            if p % 2 == 0 and mm % 2 == 0 and nn % 2 == 0:
                expected.append((Dp, p0, list(d0), q0, LINE_NAMES[li]))
    got = [(m["D_prime"], m["reduced_triple"]["p0"], m["reduced_triple"]["delta0"], m["reduced_triple"]["q0"], m["line"]) for m in members]
    require(got == expected, "member list differs from the independent parity scan")
    require(doc["excluded"] == {"lines_failing_parity": 5 * len(triples) - len(members)}, "excluded count")
    results = []
    for m in members:
        require(m["member"] is True, "non-member stored")
        results.append(replay_member(m, prefix, max_nodes, enumeration_nodes))
    # deduplication
    reps = doc["representatives"]
    for i, rep in enumerate(reps):
        m = members[rep["member_index"]]
        require(m["example_index"] == i and m["isometry_to_representative"] == identity(4), "representative record")
        require([rep["D_prime"], rep["level"], rep["primitive_determinant"], rep["scaling_c"]]
                == [m["D_prime"], m["level"], m["primitive_determinant"], m["scaling_c"]], "representative data")
    for m in members:
        rep = members[reps[m["example_index"]]["member_index"]]
        U = unimodular(m["isometry_to_representative"])
        require(gram_of(m["L_reduced"], U) == rep["L_reduced"], "isometry to representative")
        require(list(push(U, rep["a"], rep["b"])) == [m["a"], m["b"]], "characteristic not carried to the representative")
        require(m["level"] == rep["level"], "level differs inside an example class")
    for i, r1 in enumerate(reps):
        for r2 in reps[:i]:
            if r1["primitive_determinant"] != r2["primitive_determinant"]:
                continue
            m1, m2 = members[r1["member_index"]], members[r2["member_index"]]
            for U in box_isometries(m1["L_reduced"], m2["L_reduced"], max_nodes):
                require(list(push(U, m1["a"], m1["b"])) != [m2["a"], m2["b"]], "two representatives are the same example")
    # tables
    per_D = {}
    for Dp in range(Dmin, Dmax + 1):
        trip = sum(1 for t in triples if t[0] == Dp)
        mem = [m for m in members if m["D_prime"] == Dp]
        rr = [r for r in reps if r["D_prime"] == Dp]
        l2 = [r for r in rr if r["level"] == 2]
        per_D[str(Dp)] = {"reduced_triples": trip, "line_candidates": 5 * trip, "members": len(mem),
                          "distinct_examples": len(rr), "level_1": sum(r["level"] == 1 for r in rr), "level_2": len(l2),
                          "non_primitive_members": sum(m["scaling_c"] > 1 for m in mem),
                          "smallest_level_2_primitive_determinant": min((r["primitive_determinant"] for r in l2), default=None)}
    require({str(k): v for k, v in doc["per_D_prime"].items()} == per_D, "per-D' table")
    l2 = [r for r in reps if r["level"] == 2]
    require(doc["level_2_examples"] == l2, "level-2 list")
    s = doc["summary"]
    require([s["reduced_triples"], s["members"], s["distinct_examples"], s["level_1_examples"], s["level_2_examples"],
             s["smallest_level_2_primitive_determinant"], s["accidental_zero_pairs"], s["level_3_or_4_candidates"]]
            == [len(triples), len(members), len(reps), sum(r["level"] == 1 for r in reps), len(l2),
                min((r["primitive_determinant"] for r in l2), default=None),
                sum(len(m["accidental_zero_pairs"]) for m in members), sum(len(m["level_3_or_4_pairs"]) for m in members)], "summary")
    return {"verdict": "verified", "D_prime_range": [Dmin, Dmax], "reduced_triples": len(triples), "members": len(members),
            "distinct_examples": len(reps), "level_2_examples": len(l2),
            "smallest_level_2_primitive_determinant": s["smallest_level_2_primitive_determinant"]}


def decomposable(G, group):
    I = identity(4)
    for T in group:
        if T == I or T == [[-x for x in row] for row in I] or multiply(T, T) != I:
            continue
        if all((x + y) % 2 == 0 for rT, rI in zip(T, I) for x, y in zip(rT, rI)):
            return T
    return None


def replay_certificate(cert, max_nodes, enumeration_nodes):
    require(cert["schema"] == "work7-rank4-eisenstein-v1-level2-certificate", "certificate schema")
    mem = cert["member"]
    p, q, (mm, nn) = integer(mem["p"], "p"), integer(mem["q"], "q"), mem["delta"]
    G0_std = trace_gram(p, q, mm, nn)
    require(cert["G0_standard"] == G0_std and p == 2 * mem["a"] and [mm, nn] == [2 * x for x in mem["beta"]], "member parameters")
    require(mem["D_prime"] == 3 * mem["a"] * q - 2 * enorm(tuple(mem["beta"])), "D'")
    c = integer(cert["scaling_c"], "c")
    require(c >= 1 and all(x % c == 0 for row in G0_std for x in row), "scaling")
    G0 = [[x // c for x in row] for row in G0_std]
    require(cert["G0"] == G0, "G0")
    a_un, b_un, GL4 = line_characteristic(G0_std)
    require(all(x % (4 * c) == 0 for row in GL4 for x in row), "L divisibility")
    GL = [[x // (4 * c) for x in row] for row in GL4]
    require(cert["L_gram_in_kernel_coordinates"] == GL and gcd_of(GL) == 1 and cert["primitive"] is True, "primitive L")
    U = unimodular(cert["lll_transform"])
    G = matrix(cert["G"], "G", 4)
    require(gram_of(GL, U) == G, "reduced Gram")
    _, pivots = ldl(G)
    det = det4(G)
    require(cert["determinant"] == det == (mem["D_prime"] ** 2) // c ** 4 and cert["ldl_pivots"] == [str(x) for x in pivots], "determinant")
    a, b = vector(cert["a"], "a", 4), vector(cert["b"], "b", 4)
    require(list(transport(U, a_un, b_un)) == [a, b], "characteristic")
    require(all(x in (0, 1) for x in a + b) and any(a) and any(b) and dot(a, b) % 2 == 0, "characteristic parity")
    K = matrix(cert["kernel_basis"], "K", 4)
    require(abs(det4(K)) == 2 and gram_of(G, K) == G0 and all(dot(col, b) % 2 == 0 for col in transpose(K)), "kernel basis")
    require(cert["d0"] == det4(G0) == 4 * det, "d0")
    t = vector(cert["complement_t"], "t", 4)
    require(dot(t, b) % 2 == 1, "complement")
    Kinv = inverse(K)
    classes = {}
    for name, v in (("v1", a), ("v2", [x + 2 * y for x, y in zip(a, t)])):
        xv = matvec(Kinv, v)
        require(all(x.denominator == 1 for x in xv), "class in kernel")
        classes[name] = [int(x) % 2 for x in xv]
        C = matrix(cert["coset_lattices"]["O_" + name]["doubled_basis"], name, 4)
        require(abs(det4(C)) == 16, "coset index")
        for col in transpose(C):
            xc = matvec(Kinv, col)
            require(all(x.denominator == 1 for x in xc), "coset basis outside the kernel")
            require(all(int(x) % 2 == 0 for x in xc) or all((int(x) - int(y)) % 2 == 0 for x, y in zip(xc, xv)), "coset basis")
        g = gram_of(G, C)
        require(all(x % 4 == 0 for row in g for x in row), "coset lattice not integral")
        O = [[x // 4 for x in row] for row in g]
        require(cert["coset_lattices"]["O_" + name]["gram"] == O and det4(O) == det, "coset Gram")
    require(cert["classes"] == classes, "classes")
    v1, v2 = classes["v1"], classes["v2"]
    # route (i)
    rot = cert["identity"]["kernel_rotation"]
    R0 = matrix(rot["R0"], "R0", 4)
    require(R0 == R and gram_of(G0, R0) == G0 and multiply(R0, multiply(R0, R0)) == identity(4), "kernel rotation")
    m12 = [x % 2 for x in matvec(R0, v1)] == v2
    m21 = [x % 2 for x in matvec(R0, v2)] == v1
    require(rot["maps_v1_to_v2"] == m12 and rot["maps_v2_to_v1"] == m21 and (m12 or m21), "rotation does not relate the cosets")
    routes = ["kernel rotation R0 (isometry of L0 carrying O_v1 onto O_v2 or back)"]
    # route (ii)
    mc = cert["identity"]["modular_cutoff"]
    B = modular_check(G, a, b, mc["modular"])
    require(mc["cutoff"] == B, "cutoff")
    cutoff = {"attempted": mc["attempted"], "cutoff": B}
    if mc["attempted"]:
        ldl_doc = mc["ldl_backend"]
        require(ldl_doc["schema"] == "work7-rank4-cutoff-v1" and ldl_doc["G"] == G and ldl_doc["a"] == a and ldl_doc["b"] == b
                and ldl_doc["run"]["norm4_bound"] == B, "cutoff header")
        ev = ldl_doc["evidence"]
        require(ev["enumeration_complete"] is True and ev["enumeration"]["complete"] is True
                and ev["enumeration"]["norm4_bound"] == B, "cutoff completeness")
        shells, nodes = box_signed_shells(G, a, b, B, enumeration_nodes)
        require(ev["enumeration"]["shells"] == shells, "cutoff shells disagree with the box backend")
        require(all(s["signed_coefficient"] == 0 for s in shells) and ev["verdict"] == "proved_zero_by_modular_cutoff", "cutoff verdict")
        require(mc["integer_backend"]["shells_sha256"] == hashlib.sha256(json.dumps(shells).encode()).hexdigest()
                and mc["integer_backend"]["vectors"] == sum(s["vector_count"] for s in shells), "integer backend record")
        cutoff.update({"verdict": "proved_zero_by_modular_cutoff", "box_nodes": nodes, "vectors": sum(s["vector_count"] for s in shells)})
        routes.append("modular cutoff B=%d" % B)
    # groups
    group = box_isometries(G, G, max_nodes)
    for T in cert["automorphisms"]:
        matrix(T, "automorphism", 4)
    require(cert["automorphisms"] == group, "full group")
    stab = []
    for T in group:
        if any((x - y) % 2 for x, y in zip(matvec(T, a), a)):
            continue
        dual = matvec(transpose(inverse(T)), b)
        if any((x - y) % 2 for x, y in zip(dual, b)):
            continue
        stab.append({"T": T, "epsilon": symmetry_data(G, a, b, T)["epsilon"]})
    require(cert["stabilizer"] == stab and stab and all(s["epsilon"] == 0 for s in stab), "stabilizer")
    require(cert["indecomposable"] is True and decomposable(G, group) is None, "indecomposability")
    group0 = box_isometries(G0, G0, max_nodes)
    perms = [[class_index(matvec(T, class_vector(k))) for k in range(16)] for T in group0]
    k1, k2 = class_index(v1), class_index(v2)
    mapping = sum(1 for pm in perms if pm[k1] == k2)
    swapping = sum(1 for pm in perms if pm[k1] == k2 and pm[k2] == k1)
    ka = cert["kernel_automorphisms"]
    require(ka["order"] == len(group0) and ka["elements_mapping_v1_to_v2"] == mapping and ka["elements_swapping_v1_v2"] == swapping
            and mapping >= 1 and swapping == 0, "kernel group action")
    seen, orbits = set(), []
    for k in range(1, 16):
        if k in seen:
            continue
        orb, frontier = {k}, [k]
        while frontier:
            x = frontier.pop()
            for pm in perms:
                if pm[x] not in orb:
                    orb.add(pm[x])
                    frontier.append(pm[x])
        seen |= orb
        orbits.append(sorted(orb))
    require(ka["orbits_on_classes"] == orbits, "orbits")
    # route (iii)
    ms = cert["identity"]["manuscript_schema"]
    if ms["available"]:
        doc = ms["document"]
        require(doc["G"] == G and doc["a"] == a and doc["b"] == b and doc["index_two_basis"] == K
                and doc["ambient_metric"] == [[x // 2 for x in row] for row in G0] and all(x % 2 == 0 for row in G0 for x in row)
                and doc["complementary_vector"] == t and doc["automorphisms"] == group, "manuscript document fields")
        replay = rank4_counterexample.replay(doc, max_nodes)
        require(replay["verdict"] == "verified_rank4_symmetry_converse_counterexample" == ms["replay_verdict"], "manuscript replay")
        routes.append("manuscript ambient-coset schema (rank4_counterexample.replay)")
    else:
        require(any(x % 2 for row in G0 for x in row), "manuscript schema wrongly declared unavailable")
    require(cert["level"] == 2, "level")
    return {"verdict": "verified_level_2", "D_prime": mem["D_prime"], "determinant": det, "G": G, "a": a, "b": b,
            "full_group_order": len(group), "kernel_group_order": len(group0), "identity_routes": routes, "cutoff": cutoff}


def replay_task_e(doc, prefix, max_nodes, enumeration_nodes):
    require(doc["schema"] == "work7-rank4-eisenstein-v1-task-e", "Task E schema")
    rows = []
    for case in doc["cases"]:
        G0 = matrix(case["G0"], "G0", 4)
        ldl(G0)
        Rm = matrix(case["R"], "R", 4)
        require(gram_of(G0, Rm) == G0 and multiply(Rm, multiply(Rm, Rm)) == identity(4) and Rm != identity(4), "order-three isometry")
        # mixed type: the fixed sublattice of R has rank two
        fixed = [[Rm[i][j] - int(i == j) for j in range(4)] for i in range(4)]
        rank = sum(1 for d in ldl_rank(fixed))
        require(rank == 2, "R must fix a rank-two sublattice")
        N = matrix(case["N"], "N", 2)
        c = integer(case["c"], "c")
        require(case["determinant"] == det4(G0) == (N[0][0] * N[1][1] - N[0][1] * N[1][0]) * 3 * c * c // 9, "glue determinant")
        cls = [vector(v, "class", 4) for v in case["line_classes"]]
        k1, k2, k3 = (class_index(v) for v in cls)
        require(k1 ^ k2 == k3 and len({k1, k2, k3}) == 3, "line")
        rp = [class_index(matvec(Rm, class_vector(k))) for k in range(16)]
        require({rp[k1], rp[k2], rp[k3]} == {k1, k2, k3} and rp[k1] != k1, "line is not an R-orbit")
        levels, group, perms, tally = classify_pairs(G0, prefix, max_nodes, enumeration_nodes)
        require(case["aut_L0_order"] == len(group) and case["tally_105"] == tally, "group or tally")
        require(levels[(min(k1, k2), max(k1, k2))] == 1 == case["line_level"], "line level")
        S = matrix(case["swapping_element"], "swap", 4)
        require(gram_of(G0, S) == G0 and S in group, "swapping element")
        ps = [class_index(matvec(S, class_vector(k))) for k in range(16)]
        require(ps[k1] == k2 and ps[k2] == k1, "stored element does not swap")
        found = False
        for s in case["sigma_candidates"]:
            if s["T"] is None:
                continue
            T = matrix(s["T"], "sigma", 4)
            iso = gram_of(G0, T) == G0
            require(s["is_isometry"] == iso, "sigma isometry flag")
            if iso:
                pt = [class_index(matvec(T, class_vector(k))) for k in range(16)]
                sw = next(x for x in case["sigma_swaps"] if x["k"] == s["k"])
                require(sw["swaps_v1_v2"] == (pt[k1] == k2 and pt[k2] == k1) and sw["fixes_v3"] == (pt[k3] == k3), "sigma action")
                # sigma acts as -1 on the R-fixed sublattice and anti-commutes with R on the complement
                require(multiply(T, Rm) == multiply(multiply(Rm, Rm), T), "sigma is not conjugate-linear")
                found = found or sw["swaps_v1_v2"]
        require(found, "no sigma swaps the line")
        rows.append({"N": N, "c": c, "determinant": det4(G0), "line_level": 1, "aut_order": len(group)})
    require(doc["lemma_holds_in_all_cases"] is True, "lemma flag")
    return {"verdict": "verified", "cases": rows}


def ldl_rank(A):
    """Pivots of a rational row reduction (rank helper)."""
    W = [[Fraction(x) for x in row] for row in A]
    pivots = []
    r = 0
    for col in range(4):
        piv = next((i for i in range(r, 4) if W[i][col] != 0), None)
        if piv is None:
            continue
        W[r], W[piv] = W[piv], W[r]
        pivots.append(W[r][col])
        for i in range(4):
            if i != r and W[i][col] != 0:
                f = W[i][col] / W[r][col]
                W[i] = [x - f * y for x, y in zip(W[i], W[r])]
        r += 1
    return pivots


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--enumeration", type=Path, nargs="*", default=[], help="task-b / task-d files")
    parser.add_argument("--certificates", type=Path, nargs="*", default=[])
    parser.add_argument("--task-e", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--max-nodes", type=int, default=2_000_000)
    parser.add_argument("--enumeration-nodes", type=int, default=200_000_000)
    parser.add_argument("--prefix", type=int, default=60)
    args = parser.parse_args()
    if args.output:
        require(not args.output.exists(), "choose a fresh replay path")
    start = time.perf_counter()
    result = {"verdict": "verified", "inputs": {}, "enumerations": [], "certificates": [], "task_e": None}
    try:
        for path in args.enumeration:
            raw = path.read_bytes()
            result["inputs"][str(path)] = hashlib.sha256(raw).hexdigest()
            result["enumerations"].append(replay_enumeration(json.loads(raw), args.max_nodes, args.enumeration_nodes))
        for path in args.certificates:
            raw = path.read_bytes()
            result["inputs"][str(path)] = hashlib.sha256(raw).hexdigest()
            doc = json.loads(raw)
            for cert in doc["certificates"]:
                result["certificates"].append(replay_certificate(cert, args.max_nodes, args.enumeration_nodes))
        if args.task_e:
            raw = args.task_e.read_bytes()
            result["inputs"][str(args.task_e)] = hashlib.sha256(raw).hexdigest()
            result["task_e"] = replay_task_e(json.loads(raw), args.prefix, args.max_nodes, args.enumeration_nodes)
    except EnumerationLimit as error:
        result = {"verdict": "unresolved", "reason": str(error)}
    except (KeyError, TypeError, IndexError, ValueError, ZeroDivisionError) as error:
        raise InvalidCertificate("malformed Eisenstein material: " + repr(error)) from error
    result["seconds"] = time.perf_counter() - start
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=1))
    return 0 if result["verdict"] == "verified" else 2


if __name__ == "__main__":
    raise SystemExit(main())
