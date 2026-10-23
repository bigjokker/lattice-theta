"""Level-4 pairing search in minus type (Astra's plan of 2026-10-22).

Target.  An even, orthogonally indecomposable, positive definite rank-four
lattice L with det L = d = 5 mod 8 and an even characteristic (a, b) with
Theta_L[a, b] = 0.  By RANK4-STAGE2-STRATEGY.md sections 5-6 any such zero has
level 4: the two odd lattices N1, N2 below are isospectral and non-isometric.
The census through d = 997 found no zero (RANK4-LOCAL-OBSTRUCTION-REVIEW.md
section 5), so the first new determinant is 1005.

Pipeline for one determinant d (all steps exact, standard library only):

 1. classes  complete list of even classes of det d: rank4_parity_search.class_list,
             whose box carries the written coverage argument COVERAGE.  Decomposable
             classes are dropped (uniform orthogonal-sum corollary).
 2. planes   V = L/2L, Q(x) = x^t G x / 2, B(x, y) = x^t G y mod 2.  The 90 relevant
             characteristics with a != c (c = G^-1 b mod 2) fall into the 15
             B-isotropic planes; by plane invariance (strategy note section 2) the
             zero/nonzero status depends only on the plane, and automorphisms reduce
             further to plane orbits.  In minus type a plane has exactly one singular
             vector; take a nonsingular (Q(a) = 1) and c singular.
 3. pair     M = {x in L : B(c, x) = 0}, t in L \\ M, O1 = M + Z a/2, O2 = M + Z(a/2 + t),
             A_i = 2 G_{O_i}: odd integral, det 16d, common even sublattice M (Gram
             2 G_M); the third index-two overlattice of M is L (Gram 2 G_L).
 4. reject   one enumeration of L per class, histogram of y mod 4L by norm, gives
             theta(N1) - theta(N2) for every plane at once (the signed characteristic
             series).  One nonzero coefficient eliminates the plane, exactly.
 5. certify  survivors only: levels N_i of 2 A_i, Sturm bound at lcm(N1, N2), direct
             theta series of A1 and A2, complete isometry search A1 -> A2,
             indecomposability of L, and an independent proof of the zero by the
             project's modular cutoff (rank4_feasibility.cutoff_build).

Tasks
  controls   the three controls (det-21 decomposable zero, Shiota q=151, det-1729
             level-4 coset pair) plus the plane-structure checks.
  regress    run steps 2-5 on stored census classes and compare with stored verdicts.
  benchmark  steps 1-5 on an already-censused determinant, compared with its chunk.
  run        the search itself on new determinants (refuses d <= 997 and needs --go).
  replay     independent replay of a run chunk.
"""
import argparse
from collections import Counter
from fractions import Fraction
import gzip
import hashlib
import os
import platform
from itertools import combinations, product
import json
from math import gcd
from multiprocessing import get_context
from pathlib import Path
import sys
import time

from exact_theta import EnumerationLimit, InvalidCertificate, complete_shells, dot, inverse, matvec, require
from rank4_feasibility import cutoff_build, determinant, isometries
from rank4_level4_search import _enumerate, gram_of, hnf_basis, kernel_basis, lll_reduce, theta_prefix
import rank4_parity_search as builder

SCHEMA = "work7-rank4-pairing-v1"
OUT = Path("results/rank4-pairing")
CENSUSED_MAX = 997                       # minus-type census complete through here
BITS = [tuple((k >> i) & 1 for i in range(4)) for k in range(16)]   # k <-> rank4_level4_search.class_index
LADDER = (200, 800, 3200)                # norm4 = y^t G y bounds for the rejection step
I4 = [[int(i == j) for j in range(4)] for i in range(4)]


def lcm(x, y):
    return x * y // gcd(x, y)


# ----------------------------------------------------------------------------
# step 2: the quadratic space L/2L and its isotropic planes
# ----------------------------------------------------------------------------

def norm(G, x):
    return dot(x, matvec(G, x))


def q_value(G, k):
    n = norm(G, list(BITS[k]))
    require(n % 2 == 0, "Q needs an even lattice")
    return (n // 2) % 2


def b_value(G, j, k):
    return dot(list(BITS[j]), matvec(G, list(BITS[k]))) % 2


def isotropic_planes(G):
    """All 2-dimensional B-isotropic subspaces, as sorted triples of class indices."""
    planes = set()
    for u in range(1, 16):
        for v in range(u + 1, 16):
            if b_value(G, u, v) == 0:
                planes.add(tuple(sorted((u, v, u ^ v))))
    return sorted(planes)


def mod2_perm(T):
    return [sum((x % 2) << i for i, x in enumerate(matvec(T, list(BITS[k])))) for k in range(16)]


def plane_orbits(planes, group):
    """Orbit representatives; for every plane the index of an element taking its rep to it."""
    perms = [mod2_perm(T) for T in group]
    index = {p: i for i, p in enumerate(planes)}
    rep_of, carrier = {}, {}
    reps = []
    for p in planes:
        if p in rep_of:
            continue
        reps.append(p)
        for g, perm in enumerate(perms):
            image = tuple(sorted(perm[k] for k in p))
            require(image in index, "automorphism does not preserve the isotropic planes")
            if image not in rep_of:
                rep_of[image], carrier[image] = p, g
    return reps, {p: {"rep": rep_of[p], "element": carrier[p]} for p in planes}


def plane_basis(G, plane):
    """Minus type: (a, c) with Q(a) = 1 and c the unique singular vector of the plane."""
    singular = [k for k in plane if q_value(G, k) == 0]
    require(len(singular) == 1, "a minus-type isotropic plane has exactly one singular vector")
    a = min(k for k in plane if k not in singular)
    return list(BITS[a]), list(BITS[singular[0]])


def characteristic_b(G, c):
    return [x % 2 for x in matvec(G, list(c))]


# ----------------------------------------------------------------------------
# step 3: the odd-odd pair and the third member
# ----------------------------------------------------------------------------

def columns(B):
    return [[B[r][j] for r in range(4)] for j in range(4)]


def coset_pair(G, a, c):
    """M, O1, O2 and L as doubled bases inside L, with the integral Grams A_i = 2 G_{O_i}."""
    a, c = list(a), list(c)
    d = determinant(G)
    require(builder.is_even(G) and d % 2 == 1, "L must be even with odd determinant")
    require(any(a) and any(c) and a != c, "a, c must be distinct nonzero classes")
    b = characteristic_b(G, c)
    require(dot(a, b) % 2 == 0, "B(a, c) must vanish")
    K = kernel_basis(b)                                     # columns span M = {x : b.x even}
    j = next(i for i in range(4) if b[i])
    t = [int(i == j) for i in range(4)]
    twoM = [[2 * x for x in col] for col in columns(K)]
    gens = {"O1": twoM + [a], "O2": twoM + [[p + 2 * q for p, q in zip(a, t)]],
            "L": twoM + [[2 * x for x in t]], "M": twoM}
    out = {"G": G, "determinant": d, "a": a, "c": c, "b": b, "t": t, "kernel_basis": K}
    for name, g in gens.items():
        C = hnf_basis(g)
        raw = gram_of(G, C)                                 # = 4 G_name
        require(all(x % 2 == 0 for row in raw for x in row), "doubled Gram must be even")
        A = [[x // 2 for x in row] for row in raw]          # = 2 G_name
        out[name] = {"doubled_basis": C, "gram": A, "determinant": determinant(A)}
    for name in ("O1", "O2"):
        U, R = lll_reduce(out[name]["gram"])
        out[name]["reduced"] = {"transform": U, "gram": R}
    return out


def check_configuration(pair):
    """The algebraic facts of step 3; raises on any failure."""
    d = pair["determinant"]
    A1, A2, AM, AL = (pair[k]["gram"] for k in ("O1", "O2", "M", "L"))
    require(pair["O1"]["determinant"] == pair["O2"]["determinant"] == 16 * d, "det A_i = 16 d")
    require(pair["L"]["determinant"] == 16 * d and pair["M"]["determinant"] == 64 * d, "det 2G_L, 2G_M")
    require(not builder.is_even(A1) and not builder.is_even(A2), "N1, N2 must be odd")
    require(builder.is_even(AL) and builder.is_even(AM), "L and M must be even in this normalization")
    lattices = {k: pair[k]["doubled_basis"] for k in ("O1", "O2", "L", "M")}
    inv = {k: inverse(B) for k, B in lattices.items()}

    def contains(outer, inner):
        return all(x.denominator == 1 for row in (matvec(inv[outer], col) for col in columns(lattices[inner])) for x in row)
    for k in ("O1", "O2", "L"):
        require(contains(k, "M"), f"M must lie in {k}")
        require(abs(determinant(lattices[k])) * 2 == abs(determinant(lattices["M"])), f"[{k}:M] = 2")
    require(not contains("O1", "O2") and not contains("O1", "L") and not contains("O2", "L"),
            "the three index-two overlattices must be distinct")
    # the even sublattice of N_i is M: its basis vectors have even A_i-norm, and the coset vector is odd
    for k in ("O1", "O2"):
        C = lattices[k]
        odd_cols = [j for j in range(4) if norm(pair["G"], [x for x in (row[j] for row in C)]) // 2 % 2]
        require(len(odd_cols) >= 1, f"{k} has an odd vector")
    return {"det_A": 16 * d, "N_odd": True, "M_index_two": True, "distinct": True}


# ----------------------------------------------------------------------------
# step 4: shared rejection by the mod-4L histogram of L
# ----------------------------------------------------------------------------

def mod4_histogram(G, bound, max_nodes=200_000_000):
    """hist[r][q]: number of y in L with y^t G y = q and y = r mod 4 (r in 0..255)."""
    hist = [[0] * (bound + 1) for _ in range(256)]
    G00 = G[0][0]

    def inner(lo, hi, step, c, C, y):
        base = ((y[1] & 3) << 2) | ((y[2] & 3) << 4) | ((y[3] & 3) << 6)
        two_c = 2 * c
        for y0 in range(lo, hi + 1):
            hist[base | (y0 & 3)][(G00 * y0 + two_c) * y0 + C] += 1

    nodes = _enumerate(G, bound, None, inner, max_nodes)
    if builder.is_even(G):
        require(not any(row[q] for row in hist for q in range(1, bound + 1, 2)),
                "an even lattice produced an odd norm: histogram bug")
    return hist, nodes


def coset_counts(hist, a):
    """Vector counts of the coset a + 2L by norm4, from the mod-4L histogram."""
    out = [0] * len(hist[0])
    for r in range(256):
        if any((((r >> (2 * i)) & 3) - a[i]) % 2 for i in range(4)):
            continue
        for q, n in enumerate(hist[r]):
            if n:
                out[q] += n
    return out


def signed_series(hist, a, b):
    """sum over y = a mod 2 of (-1)^(b.(y-a)/2) q^(y^t G y): theta(N1) - theta(N2) at A-norm q/2."""
    out = [0] * len(hist[0])
    for r in range(256):
        rv = [(r >> (2 * i)) & 3 for i in range(4)]
        if any((x - y) % 2 for x, y in zip(rv, a)):
            continue
        sign = -1 if dot([((x - y) // 2) % 2 for x, y in zip(rv, a)], b) % 2 else 1
        for q, n in enumerate(hist[r]):
            if n:
                out[q] += sign * n
    return out


def first_nonzero(series):
    return next(((q, v) for q, v in enumerate(series) if v), None)


# ----------------------------------------------------------------------------
# step 5: Sturm certificate for a pair of integral quaternary lattices
# ----------------------------------------------------------------------------

def level_of_even(E):
    """Least N with N E^{-1} integral with even diagonal (E even, positive definite)."""
    require(builder.is_even(E), "level is defined here for even Gram matrices")
    inv = inverse(E)
    N = 1
    for i in range(4):
        for j in range(4):
            x = inv[i][j] / 2 if i == j else inv[i][j]
            N = lcm(N, Fraction(x).denominator)
    return N


def factor(n):
    out, p = [], 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            out.append((p, e))
        p += 1
    if n > 1:
        out.append((n, 1))
    return out


def sturm_bound(N, k=2):
    """floor(k/12 [SL2(Z) : Gamma0(N)]); PARI's mfsturm([N, k])."""
    index = 1
    for p, e in factor(N):
        index *= p ** (e - 1) * (p + 1)
    return k * index // 12


def theta_certificate(A1, A2, max_nodes=500_000_000):
    """Exact decision of theta(A1) == theta(A2), theta_A = sum q^(x^t A x).

    theta_A = theta_{2A} in the normalization q^(x^t (2A) x / 2), a weight-2 form on
    Gamma0(level(2A)) with character (det(2A)/.).  Equal determinants give equal
    characters; at N = lcm of the two levels, equality of the coefficients of
    q^0..q^sturm_bound(N) proves equality (Sturm).
    """
    require(determinant(A1) == determinant(A2), "different determinants cannot be isospectral here")
    E1, E2 = [[2 * x for x in r] for r in A1], [[2 * x for x in r] for r in A2]
    N1, N2 = level_of_even(E1), level_of_even(E2)
    N = lcm(N1, N2)
    bound = sturm_bound(N)
    t1, t2 = theta_prefix(A1, bound, max_nodes), theta_prefix(A2, bound, max_nodes)
    diff = next((q for q in range(bound + 1) if t1[q] != t2[q]), None)
    return {"levels": [N1, N2], "N": N, "character_discriminant": 16 * determinant(A1),
            "sturm_bound": bound, "equal": diff is None, "first_difference": diff,
            "coefficients_compared": bound + 1, "theta_head": t1[:min(len(t1), 40)]}


def two_adic_profile(A):
    """(s, e, odd part of det mod 8) of A/2, the similarity fingerprint of the local obstruction review."""
    def v2(x):
        k = 0
        while x % 2 == 0:
            x //= 2
            k += 1
        return k
    det = determinant(A)
    s = min(v2(abs(x)) for r in A for x in r if x) - 1          # of A/2
    e = v2(det) - 4 - 4 * s
    odd = det >> v2(det)
    return {"s": s, "e": e, "odd_det_mod_8": odd % 8}


def target_filter(A1, A2):
    """Necessary conditions on a candidate pair (A_i odd integral, det 16d, d = 5 mod 8)."""
    reasons = []
    d1, d2 = determinant(A1), determinant(A2)
    if d1 != d2:
        reasons.append("determinants differ")
    if builder.is_even(A1) or builder.is_even(A2):
        reasons.append("A_i must be odd")
    if d1 % 16 or (d1 // 16) % 8 != 5:
        reasons.append("det A_i must be 16 d with d = 5 mod 8")
    for A in (A1, A2):
        prof = two_adic_profile(A)
        if prof != {"s": -1, "e": 4, "odd_det_mod_8": 5}:
            reasons.append(f"2-adic profile of A/2 is {prof}, need s=-1, e=4, odd det = 5 mod 8")
            break
    return {"passes": not reasons, "reasons": reasons}


VERDICTS = {
    "nonzero": "nonzero: theta(A1) != theta(A2) below the Sturm bound",
    "hit": "HIT: zero with non-isometric pair in an indecomposable minus-type lattice",
    "discrepancy": "DISCREPANCY: zero with isometric pair in an indecomposable minus-type lattice "
                   "(contradicts RANK4-STAGE2-STRATEGY.md sections 5-6; investigate)",
    "decomposable": "zero in a decomposable lattice (outside the target)",
    "other_type": "zero outside minus type (outside the target)",
    "unresolved": "unresolved: an enumeration hit its resource limit",
}


DEFAULT_LIMITS = {"max_nodes": 5_000_000,        # automorphism groups, isometry searches, decomposability
                  "enum_nodes": 200_000_000,     # rejection histograms (build) and replay enumerations
                  "theta_nodes": 500_000_000,    # survivor theta series and signed coefficient
                  "cutoff_nodes": 400_000_000}   # signed-theta modular cutoff of a zero


def limits_of(limits=None, **overrides):
    out = dict(DEFAULT_LIMITS)
    out.update(limits or {})
    out.update({k: v for k, v in overrides.items() if v is not None})
    require(set(out) == set(DEFAULT_LIMITS) and all(isinstance(v, int) and v > 0 for v in out.values()),
            "limits must be positive integers for max_nodes, enum_nodes, theta_nodes, cutoff_nodes")
    return out


def certify_plane(G, a, c, group=None, max_nodes=None, theta_nodes=None, cutoff_nodes=None, limits=None, sink=None):
    """Certificate for one (a, c); used for survivors of step 4 and for controls.

    Unequal theta series below the Sturm bound: nonzero, with the differing
    coefficient and its signed counterpart (no isometry search or cutoff needed).
    Equal: a zero, proved a second time by the signed-theta modular cutoff, then
    classified (hit / discrepancy / decomposable / other type).  A resource limit at
    any stage gives "unresolved" with the evidence gathered so far, never a verdict;
    a contradictory certificate still raises, after the evidence has been written to the event
    sink (a ProofDisagreement: an error requiring investigation, never a verdict).  Every zero
    (hit, discrepancy, decomposable, other type) and every unresolved survivor is written to the
    sink as soon as it is classified, before the determinant finishes or is replayed.
    """
    lim = limits_of(limits, max_nodes=max_nodes, theta_nodes=theta_nodes, cutoff_nodes=cutoff_nodes)
    start = time.perf_counter()
    pair = coset_pair(G, a, c)
    conf = check_configuration(pair)
    A1, A2 = pair["O1"]["reduced"]["gram"], pair["O2"]["reduced"]["gram"]
    out = {"G": G, "determinant": pair["determinant"], "a": pair["a"], "c": pair["c"], "b": pair["b"],
           "t": pair["t"], "A1": A1, "A2": A2, "configuration": conf, "limits": lim}
    stage = "theta"
    try:
        theta = theta_certificate(A1, A2, lim["theta_nodes"])
        out["theta"] = theta
        if not theta["equal"]:
            stage = "nonzero evidence"
            n = theta["first_difference"]
            t1, t2 = theta_prefix(A1, n, lim["theta_nodes"]), theta_prefix(A2, n, lim["theta_nodes"])
            signed = complete_shells(G, pair["a"], pair["b"], 2 * n, lim["theta_nodes"])
            coeff = {s["norm4"]: s["signed_coefficient"] for s in signed["shells"]}.get(2 * n, 0)
            if not coeff == t1[n] - t2[n] != 0:
                _disagreement(sink, out, "theta difference vs signed coefficient",
                              {"A_norm": n, "theta_A1": t1[n], "theta_A2": t2[n], "signed_coefficient": coeff})
            out["nonzero_evidence"] = {"A_norm": n, "theta_A1": t1[n], "theta_A2": t2[n], "norm4": 2 * n,
                                       "signed_coefficient": coeff}
            return out | {"kind": "nonzero", "verdict": VERDICTS["nonzero"], "seconds": time.perf_counter() - start}
        stage = "automorphism group / decomposability"
        if group is None:
            group = isometries(G, G, lim["max_nodes"])
        dec = builder.decomposability(G, group, lim["max_nodes"])
        out["decomposable"] = dec["decomposable"]
        stage = "signed-theta cutoff"
        cut = cutoff_build(G, pair["a"], pair["b"], max_nodes=lim["cutoff_nodes"])
        ev = cut["evidence"]
        if not ev["enumeration_complete"]:
            raise EnumerationLimit(ev["reason"])
        if ev["verdict"] != "proved_zero_by_modular_cutoff":
            _disagreement(sink, out, "pair theta equal but signed-theta cutoff not zero",
                          {"cutoff_verdict": ev["verdict"], "cutoff_certificate": cut})
        out["signed_theta_cutoff"] = cut
        stage = "isometry A1 -> A2"
        iso = isometries(A1, A2, lim["max_nodes"])
    except EnumerationLimit as error:
        res = out | {"kind": "unresolved", "verdict": VERDICTS["unresolved"], "stage": stage,
                     "reason": f"{stage}: {error}", "seconds": time.perf_counter() - start}
        if sink:
            sink("unresolved_survivor", {"G": G, "a": res["a"], "c": res["c"], "stage": stage,
                                         "reason": res["reason"], "limits": lim, "certificate": res})
        return res
    arf = builder.arf_type(G)
    if dec["decomposable"]:
        kind = "decomposable"
    elif arf["type"] != "minus":
        kind = "other_type"
    else:
        kind = "discrepancy" if iso else "hit"
    res = out | {"kind": kind, "verdict": VERDICTS[kind], "isometric": bool(iso), "isometry": iso[0] if iso else None,
                 "arf": arf, "target_filter": target_filter(A1, A2), "seconds": time.perf_counter() - start}
    if sink:
        sink("candidate", {"candidate_kind": kind, "G": G, "a": res["a"], "c": res["c"], "b": res["b"],
                           "A1": A1, "A2": A2, "certificate": res,
                           "verification": "unverified: written before the independent replay"})
    return res


def _disagreement(sink, partial, what, evidence):
    """Save everything known, then raise: two independent computations of the same fact disagree."""
    if sink:
        sink("disagreement", {"what": what, "G": partial["G"], "a": partial["a"], "c": partial["c"],
                              "b": partial["b"], "A1": partial["A1"], "A2": partial["A2"],
                              "limits": partial["limits"], "theta": partial.get("theta"), "evidence": evidence})
    raise ProofDisagreement(f"{what} (G={partial['G']}, a={partial['a']}, c={partial['c']})")


# ----------------------------------------------------------------------------
# one class, one determinant
# ----------------------------------------------------------------------------

def decide_class(G, group=None, ladder=LADDER, max_nodes=None, enum_nodes=None,
                 keep_decomposable=False, certify=True, cutoff_nodes=None, limits=None, sink=None):
    lim = limits_of(limits, max_nodes=max_nodes, enum_nodes=enum_nodes, cutoff_nodes=cutoff_nodes)
    max_nodes, enum_nodes = lim["max_nodes"], lim["enum_nodes"]
    start = time.perf_counter()
    d = determinant(G)
    require(builder.is_even(G) and d % 8 == 5, "minus-type search needs even L with det = 5 mod 8")
    require(builder.arf_type(G)["type"] == "minus", "det = 5 mod 8 must give minus type")
    if group is None:
        group = isometries(G, G, max_nodes)
    dec = builder.decomposability(G, group, max_nodes)
    row = {"G": G, "determinant": d, "aut_order": len(group), "decomposable": dec["decomposable"]}
    if dec["decomposable"] and not keep_decomposable:
        return row | {"status": "skipped_decomposable", "seconds": time.perf_counter() - start}
    planes = isotropic_planes(G)
    require(len(planes) == 15, "a nondegenerate 4-dim symplectic F2-space has 15 Lagrangian planes")
    reps, carriers = plane_orbits(planes, group)
    pending = {p: plane_basis(G, p) for p in reps}
    decided = {}
    near = None
    for bound in ladder:
        hist, nodes = mod4_histogram(G, bound, enum_nodes)
        if near is None:
            near = near_misses(G, hist, reps, bound)
            if sink:                                 # written now, so a later failure cannot lose them
                for nm in near:
                    sink("near_miss", {"G": G, **nm,
                                       "note": "telemetry only: balanced >= 2 occupied shells in all three shifts"})
        for p in list(pending):
            # Three shifts: each nonzero vector of the plane as a (c any other one; c and c + a give the
            # same vanishing pattern).  A zero must vanish in all of them (plane invariance), so the
            # smallest nonzero witness among the three rejects the plane.  RANK4-PAIRING-SEARCH-DISCOVERY-REVIEW.md s.2
            best = None
            for k in p:
                a = list(BITS[k])
                c = list(BITS[next(x for x in p if x != k)])
                b = characteristic_b(G, c)
                hit = first_nonzero(signed_series(hist, a, b))
                if hit and (best is None or hit[0] < best[0][0]):
                    best = (hit, a, c, b)
            if best:
                hit, a, c, b = best
                require(hit[0] % 2 == 0, "odd witness norm in an even lattice")
                decided[p] = {"plane": list(p), "a": a, "c": c, "b": b, "verdict": "nonzero", "shifts_checked": 3,
                              "norm4": hit[0], "coefficient": hit[1], "A_norm": hit[0] // 2, "ladder_bound": bound}
                del pending[p]
        if not pending:
            break
    survivors = []
    for p, (a, c) in pending.items():
        rec = {"plane": list(p), "a": a, "c": c, "kind": "survived_prefix", "verdict": "survived_prefix",
               "ladder_bound": ladder[-1]}
        if certify:
            rec["certificate"] = certify_plane(G, a, c, group, limits=lim, sink=sink)
            rec["kind"], rec["verdict"] = rec["certificate"]["kind"], rec["certificate"]["verdict"]
        survivors.append(rec)
    kinds = Counter(s["kind"] for s in survivors)
    return row | {"status": "decided", "planes": 15, "plane_orbits": len(reps), "near_misses": near or [],
                  "orbit_map": {",".join(map(str, p)): {"rep": list(v["rep"]), "element": v["element"]}
                                for p, v in carriers.items()},
                  "automorphisms": group, "decided": list(decided.values()), "survivors": survivors,
                  "hits": kinds["hit"], "survivor_kinds": dict(kinds),
                  "max_first_nonzero_norm4": max((v["norm4"] for v in decided.values()), default=None),
                  "seconds": time.perf_counter() - start}


def near_misses(G, hist, planes, bound, min_depth=2):
    """Planes whose leading occupied shells balance in all three shifts, at least min_depth of them.

    Telemetry only (RANK4-PAIRING-SEARCH-PREFLIGHT-RECOMMENDATIONS.md s.2E); never part of a decision.
    Empty norm bins are not shells.  depth_is_lower_bound: some shift balanced every occupied shell
    found within the histogram bound, so its true depth may be larger.
    """
    out = []
    for p in planes:
        shifts = []
        for k in p:
            a = list(BITS[k])
            c = list(BITS[next(x for x in p if x != k)])
            series = signed_series(hist, a, characteristic_b(G, c))
            counts = coset_counts(hist, a)
            occupied = [q for q in range(len(counts)) if counts[q]]
            prefix = next((i for i, q in enumerate(occupied) if series[q]), len(occupied))
            shifts.append({"a": a, "c": c, "balanced_prefix": prefix, "occupied_within_bound": len(occupied),
                           "profile": [(q, counts[q], series[q]) for q in occupied[:prefix + 1]]})
        depth = min(x["balanced_prefix"] for x in shifts)
        if depth >= min_depth:
            out.append({"plane": list(p), "joint_depth": depth, "histogram_bound": bound,
                        "depth_is_lower_bound": any(x["balanced_prefix"] == x["occupied_within_bound"] for x in shifts),
                        "shifts": shifts})
    return out


ZERO_KINDS = ("hit", "discrepancy", "decomposable", "other_type")


def summarize(classes, rows):
    """Chunk summary; replay recomputes it from checked evidence and requires equality."""
    kinds = Counter(s["kind"] for r in rows for s in r.get("survivors", []))
    return {"candidates": classes["candidates"], "classes": classes["classes"],
            "decomposable": sum(r["decomposable"] for r in rows),
            "indecomposable": sum(not r["decomposable"] for r in rows),
            "plane_orbits": sum(r.get("plane_orbits", 0) for r in rows),
            "survivors": sum(len(r.get("survivors", [])) for r in rows),
            "survivors_nonzero": kinds["nonzero"],
            "zeros": sum(kinds[k] for k in ZERO_KINDS),
            "hits": kinds["hit"], "discrepancies": kinds["discrepancy"],
            "unresolved": kinds["unresolved"] + kinds["survived_prefix"],
            "max_first_nonzero_norm4": max((r.get("max_first_nonzero_norm4") or 0 for r in rows), default=0)}


_SPHERES = {}


def _sphere_cached(H, radius, max_nodes):
    """Vectors of H with norm <= radius grouped by norm, and H v; cached per H (the stored
    representatives are compared against many candidates).  A larger radius recomputes."""
    from lattice_census import sphere
    key = tuple(map(tuple, H))
    hit = _SPHERES.get(key)
    if hit is None or hit[0] < radius:
        vectors = sphere(H, radius, max_nodes)
        Hv = {v: matvec(H, v) for v in vectors}
        by_norm = {}
        for v in vectors:
            by_norm.setdefault(dot(v, Hv[v]), []).append(v)
        hit = _SPHERES[key] = (radius, by_norm, Hv)
    return hit[1], hit[2]


def first_isometry(G, H, max_nodes=1_000_000):
    """[U] for the first unimodular U with U^t H U = G found, or [] (rank4_feasibility.isometries, early exit).

    Same candidate images and pruning as the complete search; it stops at the first match,
    which is all the class-list deduplication uses (maps[0]).  The complete search is kept for
    automorphism groups, and verify_rank4_parity.replay_class_list re-checks every transport
    and the pairwise non-isometry of the representatives independently.
    """
    from rank4_feasibility import admit
    admit(G)
    admit(H)
    by_norm, Himages = _sphere_cached(H, max(G[i][i] for i in range(4)), max_nodes)
    images = [by_norm.get(G[i][i], []) for i in range(4)]
    nodes = 0

    def extend(cols):
        nonlocal nodes
        if len(cols) == 4:
            T = [[cols[j][i] for j in range(4)] for i in range(4)]
            return T if abs(determinant(T)) == 1 else None
        j = len(cols)
        for v in images[j]:
            nodes += 1
            if nodes > max_nodes:
                raise EnumerationLimit("rank-four isometry search interrupted")
            if all(dot(u, Himages[v]) == G[i][j] for i, u in enumerate(cols)):
                found = extend(cols + [list(v)])
                if found:
                    return found
        return None

    T = extend([])
    return [T] if T else []


def class_list(d, max_nodes=5_000_000, fast=True):
    """builder.class_list for one determinant; with fast=True the dedup stops at the first isometry.

    builder.class_list calls the module-level name `isometries` only for the pairwise dedup
    searches, and uses only maps[0]; the patch is undone afterwards.
    """
    if not fast:
        return builder.class_list(d, even_only=True, det_filter="odd", dmin=d, max_nodes=max_nodes)
    original = builder.isometries
    builder.isometries = first_isometry
    _SPHERES.clear()
    try:
        doc = builder.class_list(d, even_only=True, det_filter="odd", dmin=d, max_nodes=max_nodes)
    finally:
        builder.isometries = original
        _SPHERES.clear()
    doc["dedup"] += "; this list: first-isometry search (rank4_pairing_search.first_isometry)"
    return doc


# ----------------------------------------------------------------------------
# durable records: atomic writes, provenance, events, checkpoints
# (RANK4-PAIRING-SEARCH-PREFLIGHT-RECOMMENDATIONS.md section 2, items A-C)
# ----------------------------------------------------------------------------

class ProofDisagreement(InvalidCertificate):
    """Two independent computations of the same zero/nonzero fact disagree: investigate, never accept."""


class DeterminantIncomplete(Exception):
    """Some classes of a determinant failed; the completed class list and rows stay in the checkpoint."""


SOURCES = ("rank4_pairing_search.py", "rank4_parity_search.py", "rank4_level4_search.py", "rank4_feasibility.py",
           "exact_theta.py", "verify_rank4_parity.py", "verify_rank4_feasibility.py", "verify_rank4_level4.py",
           "verify_lattice_census.py", "lattice_census.py")
HERE = Path(__file__).resolve().parent


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sha256_json(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


_PROVENANCE = {}


def provenance():
    """Schema, source hashes, Python and platform (computed once per process)."""
    if not _PROVENANCE:
        _PROVENANCE.update(schema=SCHEMA, python=sys.version.split()[0], platform=platform.platform(),
                           sources_sha256={s: sha256_file(HERE / s) for s in SOURCES if (HERE / s).exists()})
    return dict(_PROVENANCE)


def atomic_write_json(path, obj, gz=False, indent=None):
    """Write to a temporary file in the same directory, then os.replace: readers see the old or the
    new file, never a partial one.  The temporary file is removed if anything fails."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f"{path.name}.{os.getpid()}.tmp")
    try:
        if gz:
            with gzip.open(tmp, "wt", encoding="utf-8") as f:
                json.dump(obj, f, default=str)
        else:
            tmp.write_text(json.dumps(obj, indent=indent, default=str) + "\n", encoding="utf-8")
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            tmp.unlink()


def write_chunk_file(path, doc):
    atomic_write_json(path, doc, gz=True)


class EventSink:
    """One file per event, named by determinant, kind and a hash of its content (so repeats are
    idempotent and parallel workers never share a file), written atomically at once."""

    def __init__(self, out_dir, determinant):
        self.dir = Path(out_dir) / "events"
        self.d = determinant
        self.written = []

    def __call__(self, kind, payload):
        body = {"schema": SCHEMA + "-event", "kind": kind, "determinant": self.d, **payload}
        key = sha256_json(body)[:16]
        path = self.dir / f"det{self.d:05d}-{kind}-{key}.json"
        if not path.exists():
            atomic_write_json(path, body | {"event_id": key, "written": time.strftime("%Y-%m-%d %H:%M:%S"),
                                            "provenance": provenance()}, indent=1)
        self.written.append(str(path).replace("\\", "/"))
        return path


def events_of(out_dir, d):
    """All event files of one determinant (from any run), by kind."""
    out = {}
    for f in sorted((Path(out_dir) / "events").glob(f"det{d:05d}-*.json")):
        kind = f.name.split("-")[1]
        out.setdefault(kind, []).append(str(f).replace("\\", "/"))
    return out


CHECKPOINT_SECONDS = 2.0               # at most this much finished class work is lost in a crash


def checkpoint_path(out_dir, d):
    return Path(out_dir) / "work" / f"det{d:05d}.work.json.gz"


def run_determinant(d, out_dir=OUT, max_nodes=None, ladder=LADDER, cutoff_nodes=None,
                    keep_decomposable=False, limits=None):
    """One chunk, checkpointed after the class list and after every class.

    A failure of the class list (the class universe is then not known complete) leaves only a
    failure event; nothing partial is kept as a coverage claim.  A failure of one class (resource
    limit, proof disagreement, anything else) is recorded with its G and stage; the remaining
    classes are still processed, the checkpoint keeps every finished row, and DeterminantIncomplete
    is raised.  Calling again resumes: finished rows are reused and only the failed or missing classes
    are recomputed (with the current limits).  Only a complete document becomes a chunk, which then
    goes to the unchanged full replay.  keep_decomposable is for tests and controls only.
    """
    require(d % 8 == 5, "det must be 5 mod 8")
    lim = limits_of(limits, max_nodes=max_nodes, cutoff_nodes=cutoff_nodes)
    out_dir = Path(out_dir)
    path = out_dir / f"det{d:05d}.json.gz"
    if path.exists():
        raise FileExistsError(path)
    sink = EventSink(out_dir, d)
    ck_path = checkpoint_path(out_dir, d)
    start = time.perf_counter()
    attempt = {"time": time.strftime("%Y-%m-%d %H:%M:%S"), "limits": lim, "provenance": provenance()}
    if ck_path.exists():
        ck = load_doc(ck_path)
        require(ck.get("schema") == SCHEMA + "-work" and ck["determinant"] == d and ck["ladder"] == list(ladder)
                and ck["keep_decomposable"] == keep_decomposable, "checkpoint does not match this run")
        require(sha256_json(ck["classes_doc"]) == ck["classes_doc_sha256"], "checkpoint class list corrupted")
        classes = ck["classes_doc"]
        rows = {int(k): v for k, v in ck["rows"].items()}
        ck["attempts"].append(attempt)
    else:
        try:
            classes = class_list(d, max_nodes=lim["max_nodes"])
        except Exception as error:
            sink("failure", {"stage": "class_list", "error": f"{type(error).__name__}: {error}", "limits": lim,
                             "note": "class universe not known complete; rerun this stage with a larger budget"})
            raise
        rows = {}
        ck = {"schema": SCHEMA + "-work", "determinant": d, "ladder": list(ladder),
              "keep_decomposable": keep_decomposable, "classes_doc": classes,
              "classes_doc_sha256": sha256_json(classes), "rows": {}, "failures": [], "attempts": [attempt]}
    reps = classes["representatives"]

    last = [0.0]

    def save(force=False):
        """Atomic checkpoint; throttled to one write per CHECKPOINT_SECONDS unless forced (failures, end)."""
        if not force and time.perf_counter() - last[0] < CHECKPOINT_SECONDS:
            return
        ck["rows"] = {str(i): r for i, r in sorted(rows.items())}
        ck["pending"] = [i for i in range(len(reps)) if i not in rows]
        atomic_write_json(ck_path, ck, gz=True)
        last[0] = time.perf_counter()
    failures = []
    save(force=True)
    for i, rep in enumerate(reps):
        if i in rows:
            continue
        try:
            rows[i] = decide_class(rep["G"], ladder=ladder, keep_decomposable=keep_decomposable, limits=lim, sink=sink)
        except Exception as error:
            failures.append({"index": i, "G": rep["G"], "error": f"{type(error).__name__}: {error}", "limits": lim,
                             "stage": "proof disagreement" if isinstance(error, ProofDisagreement) else "decide_class"})
            sink("failure", failures[-1])
            ck["failures"] = failures
            save(force=True)
            continue
        ck["failures"] = failures
        save()
    save(force=True)
    if failures:
        raise DeterminantIncomplete(f"det {d}: {len(failures)} of {len(reps)} classes failed "
                                    f"({', '.join(sorted({f['stage'] for f in failures}))}); "
                                    f"completed work kept in {ck_path}")
    ordered = [rows[i] for i in range(len(reps))]
    doc = {"schema": SCHEMA + "-chunk", "determinant": d, "classes_doc": classes, "lattices": ordered,
           "parameters": {"ladder": list(ladder), "limits": lim, "keep_decomposable": keep_decomposable},
           "provenance": provenance(), "build_attempts": ck["attempts"],
           "summary": summarize(classes, ordered), "seconds": time.perf_counter() - start}
    write_chunk_file(path, doc)
    ck_path.unlink()
    return {"determinant": d, "file": str(path).replace("\\", "/"), **doc["summary"],
            "near_misses": sum(len(r.get("near_misses", [])) for r in ordered), "events": sink.written,
            "seconds": doc["seconds"]}


def retry_unresolved(d, out_dir=OUT, limits=None):
    """Re-certify the unresolved survivors of an existing chunk with new limits.

    The original chunk is kept unchanged under out_dir/superseded/ (its hash in the file name); the new
    chunk records the retry (time, previous file hash, limits) and is then replayed like any other.
    """
    lim = limits_of(limits)
    path = Path(out_dir) / f"det{d:05d}.json.gz"
    doc = load_doc(path)
    require(doc["determinant"] == d, "chunk determinant")
    pending = [(row, s) for row in doc["lattices"] for s in row.get("survivors", []) if s["kind"] == "unresolved"]
    require(pending, f"det {d}: no unresolved survivors to retry")
    old_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    sink = EventSink(out_dir, d)
    for row, s in pending:
        cert = certify_plane(row["G"], s["a"], s["c"], row["automorphisms"], limits=lim, sink=sink)
        s["certificate"], s["kind"], s["verdict"] = cert, cert["kind"], cert["verdict"]
    for row, _ in pending:
        kinds = Counter(x["kind"] for x in row["survivors"])
        row["hits"], row["survivor_kinds"] = kinds["hit"], dict(kinds)
    doc["summary"] = summarize(doc["classes_doc"], doc["lattices"])
    doc.setdefault("retries", []).append({"time": time.strftime("%Y-%m-%d %H:%M:%S"), "limits": lim,
                                          "previous_sha256": old_hash, "survivors_retried": len(pending),
                                          "provenance": provenance()})
    keep = Path(out_dir) / "superseded" / f"det{d:05d}.{time.strftime('%Y%m%d-%H%M%S')}.{old_hash[:12]}.json.gz"
    keep.parent.mkdir(parents=True, exist_ok=True)
    keep.write_bytes(path.read_bytes())
    write_chunk_file(path, doc)
    return {"determinant": d, "retried": len(pending), "superseded_copy": str(keep).replace("\\", "/"),
            **doc["summary"]}


# ----------------------------------------------------------------------------
# independent replay of a chunk
# ----------------------------------------------------------------------------
#
# Fixes RANK4-PAIRING-SEARCH-REVIEW.md findings 1-3.  The plane step has a second
# implementation here (tuple vectors, no class-index/XOR helpers); enumerations use
# the verifier's box backends; the pair is rebuilt from G, a, c in a different basis.

def _vbits():
    return [v for v in product((0, 1), repeat=4) if any(v)]


def _vB(G, u, v):
    return sum(G[i][j] * u[i] * v[j] for i in range(4) for j in range(4)) % 2


def _vQ(G, v):
    return (sum(G[i][j] * v[i] * v[j] for i in range(4) for j in range(4)) // 2) % 2


def _vkey(vectors):
    return tuple(sorted(sum(x << i for i, x in enumerate(v)) for v in vectors))


def _vact(T, v):
    return tuple(sum(T[i][j] * v[j] for j in range(4)) % 2 for i in range(4))


def _vplanes(G):
    """All B-isotropic planes as frozensets of tuple vectors (second implementation)."""
    spans = set()
    for u, v in combinations(_vbits(), 2):
        if _vB(G, u, v) == 0:
            spans.add(frozenset((u, v, tuple((x + y) % 2 for x, y in zip(u, v)))))
    require(len(spans) == 15, "replay: 15 isotropic planes")
    return spans


def _int_hnf(gens):
    """Row echelon basis of the Z-span of integer vectors in Z^4 (independent of hnf_basis)."""
    rows = [list(g) for g in gens if any(g)]
    basis = []
    for col in range(4):
        piv = [r for r in rows if r[col]]
        rest = [r for r in rows if not r[col]]
        while len(piv) > 1:
            piv.sort(key=lambda r: abs(r[col]))
            p = piv[0]
            nxt = [p]
            for r in piv[1:]:
                q = r[col] // p[col]
                r = [x - q * y for x, y in zip(r, p)]
                (nxt if r[col] else rest).append(r)
            piv = nxt
        require(len(piv) == 1, "replay: generators do not have full rank")
        basis.append(piv[0])
        rows = [r for r in rest if any(r)]
    require(not rows, "replay: leftover generators")
    return basis


def _independent_pair(G, a, c):
    """2 G_{O1}, 2 G_{O2} from G, a, c, with M spanned by a different basis than kernel_basis."""
    b = [sum(G[i][j] * c[j] for j in range(4)) % 2 for i in range(4)]
    j = max(i for i in range(4) if b[i])                      # kernel_basis uses the first such index
    e = [[int(r == i) for r in range(4)] for i in range(4)]
    M = [e[i] if not b[i] else [x + y for x, y in zip(e[i], e[j])] for i in range(4) if i != j] + [[2 * x for x in e[j]]]
    require(all(sum(x * y for x, y in zip(m, b)) % 2 == 0 for m in M), "replay: M basis")
    t = e[j]
    out = []
    for extra in (list(a), [x + 2 * y for x, y in zip(a, t)]):
        Bs = _int_hnf([[2 * x for x in m] for m in M] + [extra])
        raw = [[sum(u[p] * G[p][q] * v[q] for p in range(4) for q in range(4)) for v in Bs] for u in Bs]
        require(all(x % 2 == 0 for r in raw for x in r), "replay: doubled Gram parity")
        out.append([[x // 2 for x in r] for r in raw])
    return out, b


def _level_even(E):
    """Least N with N E^{-1} integral, even diagonal (second implementation)."""
    inv = inverse(E)
    N = 1
    for i in range(4):
        for j in range(4):
            den = Fraction(inv[i][j] / (2 if i == j else 1)).denominator
            N = N * den // gcd(N, den)
    return N


def _sturm(N):
    index, n, p = N, N, 2
    while p * p <= n:
        if n % p == 0:
            index = index // p * (p + 1)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        index = index // n * (n + 1)
    return index // 6


def _box_theta(A, bound, nodes):
    """theta_A coefficients 0..bound by the verifier's box backend (y = 2x, norm4 = 4 x^t A x)."""
    from verify_rank4_level4 import box_signed_shells
    shells, _ = box_signed_shells(A, [0] * 4, [0] * 4, 4 * bound, nodes)
    out = [0] * (bound + 1)
    for s in shells:
        require(s["norm4"] % 4 == 0, "replay: theta parity")
        out[s["norm4"] // 4] = s["vector_count"]
    return out


def _replay_survivor(G, d, rec, group, decomposable, max_nodes, enum_nodes):
    from verify_rank4_feasibility import box_isometries, cutoff_replay
    from verify_rank4_level4 import box_signed_shells
    cert = rec.get("certificate")
    require(isinstance(cert, dict), "replay: survivor without certificate")
    a, c = list(rec["a"]), list(rec["c"])
    require(cert.get("G") == G and cert.get("a") == a and cert.get("c") == c,
            "replay: certificate does not belong to this plane")
    kind = cert.get("kind")
    require(kind in VERDICTS and rec.get("kind") == kind and rec.get("verdict") == cert.get("verdict") == VERDICTS[kind],
            "replay: survivor verdict")
    if kind == "unresolved":
        return kind
    (P1, P2), b = _independent_pair(G, a, c)
    require(cert["b"] == b, "replay: characteristic b")
    A1, A2 = cert["A1"], cert["A2"]
    for A, P in ((A1, P1), (A2, P2)):
        require(len(A) == 4 and all(len(r) == 4 for r in A), "replay: A_i shape")
        require(determinant(A) == determinant(P) == 16 * d, "replay: det A_i = 16 d")
        require(box_isometries(A, P, max_nodes), "replay: stored A_i is not the Gram of the pair built from G, a, c")
    N1, N2 = _level_even([[2 * x for x in r] for r in A1]), _level_even([[2 * x for x in r] for r in A2])
    N = N1 * N2 // gcd(N1, N2)
    bound = _sturm(N)
    th = cert["theta"]
    require(th["levels"] == [N1, N2] and th["N"] == N and th["sturm_bound"] == bound, "replay: level or Sturm bound")
    if kind == "nonzero":
        ev = cert["nonzero_evidence"]
        n = ev["A_norm"]
        require(isinstance(n, int) and 0 < n <= bound, "replay: difference index")
        t1, t2 = _box_theta(A1, n, enum_nodes), _box_theta(A2, n, enum_nodes)
        shells, _ = box_signed_shells(G, a, b, 2 * n, enum_nodes)
        coeff = {s["norm4"]: s["signed_coefficient"] for s in shells}.get(2 * n, 0)
        require(t1[n] == ev["theta_A1"] and t2[n] == ev["theta_A2"]
                and coeff == ev["signed_coefficient"] == t1[n] - t2[n] != 0, "replay: nonzero evidence")
        return kind
    t1, t2 = _box_theta(A1, bound, enum_nodes), _box_theta(A2, bound, enum_nodes)
    require(t1 == t2, "replay: claimed zero but theta(A1) != theta(A2) below the Sturm bound")
    cut = cert["signed_theta_cutoff"]
    require(cut["G"] == G and cut["a"] == a and cut["b"] == b, "replay: cutoff certificate for another characteristic")
    require(cutoff_replay(cut, enum_nodes)["verdict"] == "proved_zero_by_modular_cutoff",
            "replay: cutoff does not prove the zero")
    iso = bool(box_isometries(A1, A2, max_nodes))
    expected = "decomposable" if decomposable else ("discrepancy" if iso else "hit")
    require(kind == expected and cert["isometric"] == iso and cert["decomposable"] == decomposable,
            "replay: zero classification")
    if kind == "hit":
        require(target_filter(A1, A2)["passes"], "replay: a hit must pass the target filters")
    return kind


def replay_chunk(path, max_nodes=5_000_000, enum_nodes=200_000_000, domain=None):
    """Verdict "verified" only if every claim of the chunk is re-derived; "unresolved" if a survivor is.

    domain: an optional verify_rank4_parity.SharedDomain covering this determinant (one interval scan for a
    batch of replays).  Any malformed field is an InvalidCertificate, never a crash or an acceptance."""
    try:
        return _replay_chunk(path, max_nodes, enum_nodes, domain)
    except InvalidCertificate:
        raise
    except (KeyError, TypeError, IndexError, ValueError, AttributeError, ZeroDivisionError) as error:
        raise InvalidCertificate(f"malformed pairing chunk: {type(error).__name__}: {error}") from error


def _replay_chunk(path, max_nodes, enum_nodes, domain=None):
    import verify_rank4_parity as verifier
    from verify_rank4_feasibility import box_isometries
    from verify_rank4_level4 import box_signed_shells
    with gzip.open(path, "rt", encoding="utf-8") as f:
        doc = json.load(f)
    start = time.perf_counter()
    require(doc.get("schema") == SCHEMA + "-chunk", "wrong schema")
    d = doc["determinant"]
    cd = doc["classes_doc"]
    require(isinstance(d, int) and d % 8 == 5, "replay: chunk determinant")
    require(cd["determinant_minimum"] == cd["determinant_bound"] == d and cd["even_only"] is True
            and cd["odd_lattices_only"] is False and cd["determinant_filter"] == "odd", "replay: class-list domain")
    name = Path(path).name
    if name.startswith("det"):
        require(name == f"det{d:05d}.json.gz", "replay: file name does not match the determinant")
    keep = doc.get("parameters", {}).get("keep_decomposable", False)
    listing = verifier.replay_class_list(cd, max_nodes, domain)
    rows = doc["lattices"]
    require([r["G"] for r in rows] == [r["G"] for r in cd["representatives"]], "replay: rows do not match the class list")
    stats, kinds, max_first = Counter(), Counter(), 0
    for row in rows:
        G = row["G"]
        require(determinant(G) == d == row["determinant"] and builder.is_even(G), "replay: row determinant")
        group = box_isometries(G, G, max_nodes)
        decomposable = verifier.splitting_involution(group) is not None
        require(decomposable == (verifier.indecomposable_components(G, max_nodes) >= 2) == row["decomposable"],
                "replay: decomposability")
        require(row["aut_order"] == len(group), "replay: automorphism group order")
        stats["decomposable" if decomposable else "indecomposable"] += 1
        if row["status"] == "skipped_decomposable":
            require(decomposable, "replay: an indecomposable class was skipped")
            continue
        require(row["status"] == "decided" and (not decomposable or keep), "replay: row status")
        require(row["automorphisms"] == group, "replay: automorphism group")
        require(sum(_vQ(G, v) == 0 for v in _vbits()) == 5, "replay: L must be minus type")
        planes = {_vkey(p): p for p in _vplanes(G)}
        omap = {}
        for k, v in row["orbit_map"].items():
            key = tuple(int(x) for x in k.split(","))
            require(k == ",".join(map(str, key)), "replay: orbit-map key format")
            omap[key] = v
        require(set(omap) == set(planes), "replay: orbit map must list exactly the 15 isotropic planes")
        orbits = {frozenset(_vkey(_vact(T, v) for v in p) for T in group) for p in planes.values()}
        require(row["plane_orbits"] == len(orbits), "replay: orbit count")
        reps = {tuple(v["rep"]) for v in omap.values()}
        require(reps <= set(planes) and len(reps) == len(orbits) and all(len(o & reps) == 1 for o in orbits),
                "replay: one representative per orbit")
        for key, v in omap.items():
            e = v["element"]
            require(isinstance(e, int) and 0 <= e < len(group), "replay: carrier index")
            rep = planes[tuple(v["rep"])]
            require(_vkey(_vact(group[e], x) for x in rep) == key, "replay: carrier does not map rep to plane")
        dec_planes = [tuple(r["plane"]) for r in row["decided"]]
        sur_planes = [tuple(s["plane"]) for s in row["survivors"]]
        require(len(set(dec_planes)) == len(dec_planes) and len(set(sur_planes)) == len(sur_planes)
                and not set(dec_planes) & set(sur_planes) and set(dec_planes) | set(sur_planes) == reps,
                "replay: exactly one decision per orbit representative")
        for r in row["decided"] + row["survivors"]:
            a, c = tuple(r["a"]), tuple(r["c"])
            plane = planes[tuple(r["plane"])]
            require(a != c and a in plane and c in plane, "replay: (a, c) not in the plane")
        row_first = None
        for r in row["decided"]:
            a, c = list(r["a"]), list(r["c"])
            b = [sum(G[i][j] * c[j] for j in range(4)) % 2 for i in range(4)]
            require(r["b"] == b and r["verdict"] == "nonzero", "replay: decision record")
            q = r["norm4"]
            require(isinstance(q, int) and q > 0, "replay: norm4")
            shells, _ = box_signed_shells(G, a, b, q, enum_nodes)
            coeffs = {s["norm4"]: s["signed_coefficient"] for s in shells}
            require(all(coeffs.get(x, 0) == 0 for x in range(q)), "replay: an earlier nonzero coefficient")
            require(coeffs.get(q, 0) == r["coefficient"] != 0, "replay: stored nonzero coefficient")
            row_first = max(row_first or 0, q)
            stats["orbit_decisions"] += 1
        require(row.get("max_first_nonzero_norm4") == row_first, "replay: row max first nonzero norm4")
        max_first = max(max_first, row_first or 0)
        row_kinds = Counter()
        for s in row["survivors"]:
            row_kinds[_replay_survivor(G, d, s, group, decomposable, max_nodes, enum_nodes)] += 1
            stats["survivors"] += 1
        require(row["hits"] == row_kinds["hit"], "replay: row hit count")
        if "survivor_kinds" in row:
            require(row["survivor_kinds"] == dict(row_kinds), "replay: row survivor kinds")
        kinds += row_kinds
        stats["plane_orbits"] += len(orbits)
    recomputed = {"candidates": listing["candidates"], "classes": listing["classes"],
                  "decomposable": stats["decomposable"], "indecomposable": stats["indecomposable"],
                  "plane_orbits": stats["plane_orbits"], "survivors": stats["survivors"],
                  "survivors_nonzero": kinds["nonzero"], "zeros": sum(kinds[k] for k in ZERO_KINDS),
                  "hits": kinds["hit"], "discrepancies": kinds["discrepancy"], "unresolved": kinds["unresolved"],
                  "max_first_nonzero_norm4": max_first}
    summary = doc["summary"]
    for key in ("candidates", "classes", "decomposable", "indecomposable", "plane_orbits", "survivors", "hits"):
        require(key in summary, f"replay: summary lacks {key}")
    for key, value in summary.items():
        require(key in recomputed and recomputed[key] == value, f"replay: summary {key}")
    for key, value in recomputed.items():
        if key not in summary:      # chunks written before these fields existed
            require(value == 0, f"replay: summary lacks nonzero {key}")
    verdict = "unresolved" if kinds["unresolved"] else "verified"
    return {"file": str(path).replace("\\", "/"), "determinant": d, "verdict": verdict, **recomputed,
            "orbit_decisions_replayed": stats["orbit_decisions"], "survivors_replayed": stats["survivors"],
            "survivor_kinds": dict(kinds), "seconds": time.perf_counter() - start}


# ----------------------------------------------------------------------------
# controls and regressions
# ----------------------------------------------------------------------------

DET21 = [[2, 1, 0, 0], [1, 2, 0, 0], [0, 0, 2, 1], [0, 0, 1, 4]]
DET21_AC = ([1, 0, 1, 0], [0, 1, 0, 1])          # strategy note section 3; b = Gc = (1,0,1,0)
SHIOTA_A = [[6, 2, -1, 1], [2, 12, 5, 4], [-1, 5, 16, 6], [1, 4, 6, 28]]
SHIOTA_B = [[6, 0, -2, -3], [0, 12, 3, 4], [-2, 3, 14, 2], [-3, 4, 2, 28]]


def control_det21():
    G = DET21
    a0, c0 = DET21_AC
    require(characteristic_b(G, c0) == [1, 0, 1, 0], "det-21 characteristic")
    plane = tuple(sorted({builder_index(a0), builder_index(c0), builder_index(a0) ^ builder_index(c0)}))
    a, c = plane_basis(G, plane)                 # Q(a) = 1 basis of the same plane
    cert = certify_plane(G, a, c)
    row = decide_class(G, keep_decomposable=True, certify=False)
    # the zero plane, or the orbit representative standing for it, must survive the rejection step
    rep = tuple(row["orbit_map"][",".join(map(str, plane))]["rep"])
    require(rep in [tuple(s["plane"]) for s in row["survivors"]], "det-21 zero plane must survive rejection")
    require(cert["configuration"]["det_A"] == 16 * 21, "configuration")
    require(cert["theta"]["equal"], "det-21 pair must be isospectral (the zero)")
    cutoff_verdict = cert["signed_theta_cutoff"]["evidence"]["verdict"]
    require(cutoff_verdict == "proved_zero_by_modular_cutoff", "cutoff proves the zero")
    require(cert["decomposable"], "det-21 must fail indecomposability")
    require(cert["kind"] == "decomposable", "det-21 must be classified as a zero in a decomposable lattice")
    return {"control": "det-21 decomposable zero", "plane": list(plane), "a": a, "c": c,
            "passes_configuration": True, "theta_equal": True, "levels": cert["theta"]["levels"],
            "sturm_bound": cert["theta"]["sturm_bound"], "cutoff": cutoff_verdict, "kind": cert["kind"],
            "isometric": cert["isometric"], "decomposable": True, "target_filter": cert["target_filter"],
            "rejected_because": "decomposable", "plane_orbits": row["plane_orbits"],
            "surviving_orbits_in_class": [s["plane"] for s in row["survivors"]]}


def builder_index(v):
    return sum((x % 2) << i for i, x in enumerate(v))


def control_shiota():
    cert = theta_certificate(SHIOTA_A, SHIOTA_B)
    iso = isometries(SHIOTA_A, SHIOTA_B, 5_000_000)
    filt = target_filter(SHIOTA_A, SHIOTA_B)
    require(determinant(SHIOTA_A) == determinant(SHIOTA_B) == 151 ** 2, "Shiota determinants")
    require(cert["equal"] and not iso, "Shiota pair: isospectral and non-isometric")
    require(not filt["passes"], "Shiota must fail the target filters")
    return {"control": "Shiota q=151", "theta": cert, "isometric": False, "target_filter": filt}


def control_det1729():
    doc = json.loads(Path("examples/rank4-level4/det1729-level4.json").read_text(encoding="utf-8"))
    c = doc["certificates"][0]
    A1 = c["coset_lattices"]["O_v1"]["reduced"]["gram"]
    A2 = c["coset_lattices"]["O_v2"]["reduced"]["gram"]
    cert = theta_certificate(A1, A2)
    iso = isometries(A1, A2, 5_000_000)
    filt = target_filter(A1, A2)
    require(cert["equal"] and not iso, "det-1729 coset pair: isospectral and non-isometric")
    require(not filt["passes"], "det-1729 (odd L) must fail the target filters")
    return {"control": "det-1729 level-4 coset pair (odd L)", "theta": cert, "isometric": False, "target_filter": filt}


def load_doc(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as f:
        return json.load(f)


def regress_row(r, all_orderings=True):
    """Compare plane-based rejection with the stored verdict of every characteristic of one class."""
    G = r["G"]
    stored = {}
    for dct in r["parity"]["decisions"]:
        a = builder_index(dct["a"])
        cs = builder.solve_mod2(G, dct["b"])
        require(len(cs) == 1, "c")
        stored[(a, builder_index(cs[0]))] = dct
    planes = isotropic_planes(G)
    require(len(planes) == 15, "15 planes")
    typ = builder.arf_type(G)["type"]
    hist, _ = mod4_histogram(G, max(LADDER))
    out = Counter()
    for p in planes:
        verdicts = set()
        for a in p:
            for c in p:
                if a == c or not (all_orderings or (typ == "minus" and [list(BITS[a]), list(BITS[c])] == list(plane_basis(G, p)))):
                    continue
                dct = stored[(a, c)]
                b = characteristic_b(G, BITS[c])
                require(b == dct["b"], "b = G c mod 2")
                fn = first_nonzero(signed_series(hist, list(BITS[a]), b))
                if dct["verdict"] == "nonzero":
                    f = dct["first_nonzero"]
                    require(fn is not None and fn[0] == f["norm4"] and fn[1] == f["coefficient"],
                            f"first nonzero mismatch {fn} vs {f}")
                    verdicts.add("nonzero")
                else:
                    require(dct["verdict"].startswith("zero"), "unresolved stored verdict")
                    require(fn is None, "a stored zero has a nonzero histogram coefficient")
                    verdicts.add("zero")
                out["characteristics"] += 1
        require(len(verdicts) == 1, "plane invariance violated")
        out["planes_" + verdicts.pop()] += 1
    return out


def regress_source(path, limit=None):
    doc = load_doc(path)
    total = Counter()
    rows = [r for r in doc["lattices"] if not r["decomposable"]][:limit]
    for r in rows:
        total += regress_row(r)
        total["classes"] += 1
        total["type_" + builder.arf_type(r["G"])["type"]] += 1
    return dict(total)


def certify_stored_zeros(path):
    """Certify every stored zero plane with the full step-5 chain (Sturm, cutoff, isometry, decomposability)."""
    doc = load_doc(path)
    rows = []
    for r in doc["lattices"]:
        if r["decomposable"]:
            continue
        zeros = [d for d in r["parity"]["decisions"] if d["verdict"].startswith("zero")]
        seen = set()
        for dct in zeros:
            cs = builder.solve_mod2(r["G"], dct["b"])
            a, c = builder_index(dct["a"]), builder_index(cs[0])
            if a == c:
                continue
            p = tuple(sorted((a, c, a ^ c)))
            if p in seen:
                continue
            seen.add(p)
            nonsingular = [k for k in p if q_value(r["G"], k) == 1]
            if not nonsingular:
                # totally singular plane (plus type only): N_i would be even, outside the odd-pair construction
                rows.append({"det": r["determinant"], "type": builder.arf_type(r["G"])["type"], "plane": list(p),
                             "verdict": "totally singular plane, not certified"})
                continue
            a = nonsingular[0]
            c = next(k for k in p if k != a)
            cert = certify_plane(r["G"], list(BITS[a]), list(BITS[c]), r["automorphisms"])
            require(cert["theta"]["equal"] and cert["isometric"], "stored symmetry zeros give isometric pairs")
            rows.append({"det": r["determinant"], "type": cert["arf"]["type"], "plane": list(p),
                         "levels": cert["theta"]["levels"], "sturm": cert["theta"]["sturm_bound"],
                         "isometric": cert["isometric"], "verdict": cert["verdict"]})
    return rows


REGRESSION_SOURCES = (["results/rank4-parity/task-b-det201.json",
                       "results/rank4-minus-type/minus-det202-501.json.gz"]
                      + sorted(str(p).replace("\\", "/") for p in Path("results/rank4-gap-census/even-odd-202-279").glob("det*.json.gz"))
                      + sorted(str(p).replace("\\", "/") for p in Path("results/rank4-minus-type-ext/even-minus-509-997").glob("det*.json.gz")))


def _regress_file(path):
    start = time.perf_counter()
    out = regress_source(path)
    out["stored_zero_planes_certified"] = certify_stored_zeros(path)
    out["file"], out["seconds"] = path, time.perf_counter() - start
    return out


def regress_all(workers=16):
    """Every stored even odd-det census class: plane rejection vs stored verdicts, all six orderings;
    every stored zero plane through the full certification chain."""
    total, zeros, files = Counter(), [], 0
    with get_context("spawn").Pool(workers) as pool:
        for res in pool.imap_unordered(_regress_file, REGRESSION_SOURCES):
            files += 1
            zeros += res.pop("stored_zero_planes_certified")
            res.pop("file"); res.pop("seconds")
            total.update(res)
    return {"files": files, **dict(total), "stored_zero_planes_certified": len(zeros),
            "zero_planes_by_outcome": dict(Counter(f"{z['type']}: " + (f"certified, isometric={z['isometric']}"
                                                    if "isometric" in z else z["verdict"]) for z in zeros)),
            "zero_plane_examples": zeros[:5]}


def benchmark(d, workdir):
    """Steps 1-5 on an already-censused minus determinant, compared with the stored census chunk."""
    require(d % 8 == 5 and 509 <= d <= CENSUSED_MAX, "benchmark determinants come from the 509..997 census")
    stored = load_doc(f"results/rank4-minus-type-ext/even-minus-509-997/det{d}.json.gz")
    t0 = time.perf_counter()
    res = run_determinant(d, out_dir=workdir)
    t1 = time.perf_counter()
    s = stored["summary"]
    require(res["classes"] == s["classes"] and res["indecomposable"] == s["indecomposable_classes"],
            "class counts differ from the stored census")
    require(res["hits"] == 0 and res["survivors"] == 0, "the census found no zero here")
    return res | {"stored_classes": s["classes"], "stored_indecomposable": s["indecomposable_classes"],
                  "wall_seconds": t1 - t0}


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------

def _outcome(row):
    """no_hit / verified_hit / discrepancy / unresolved / incomplete / failed, from checked evidence only."""
    if row["status"] in ("failed", "incomplete"):
        return row["status"]
    rep = row["replay"]
    if rep["verdict"] != "verified":
        return "unresolved"
    if rep["discrepancies"]:
        return "discrepancy"
    return "verified_hit" if rep["hits"] else "no_hit"


def _run_task(task):
    """Build (or reuse) one chunk and replay it with the same limits; failures are returned, not raised."""
    d, out_dir, limits, ladder, resume, domain = task
    path = Path(out_dir) / f"det{d:05d}.json.gz"
    row = {"determinant": d, "file": str(path).replace("\\", "/")}
    start = time.perf_counter()
    try:
        if path.exists():
            require(resume, "chunk exists and --resume was not given")
            row["status"] = "reused"
        else:
            row |= run_determinant(d, out_dir=out_dir, ladder=tuple(ladder), limits=limits)
            row["status"] = "built"
        rep = replay_and_record(path, limits, "run", out_dir, domain)
        row["replay"] = {k: rep[k] for k in ("verdict", "classes", "indecomposable", "plane_orbits",
                                             "orbit_decisions_replayed", "survivors", "zeros", "hits",
                                             "discrepancies", "unresolved", "seconds")}
        row["sha256"] = sha256_file(path)
        row["locally_verified"] = rep["verdict"] == "verified"
    except DeterminantIncomplete as error:
        row |= {"status": "incomplete", "error": str(error),
                "checkpoint": str(checkpoint_path(out_dir, d)).replace("\\", "/"), "locally_verified": False}
    except Exception as error:                    # recorded in the manifest; the chunk is not counted
        row |= {"status": "failed", "error": f"{type(error).__name__}: {error}", "locally_verified": False}
    row["events"] = events_of(out_dir, d)
    row["outcome"] = _outcome(row)
    row["seconds_total"] = time.perf_counter() - start
    return row


def validate_dets(dets, allow_censused=False):
    require(bool(dets), "no determinants requested")
    require(len(set(dets)) == len(dets), "duplicate determinants requested")
    require(all(isinstance(d, int) and d % 8 == 5 for d in dets), "determinants must be 5 mod 8")
    require(allow_censused or all(d > CENSUSED_MAX for d in dets),
            f"determinants must be above {CENSUSED_MAX} (already censused)")


def run_search(dets, workers, out_dir=OUT, resume=False, limits=None, ladder=LADDER, allow_censused=False,
               _inline=False):
    """The search over dets; returns the manifest.  totals.complete is true only if every chunk verified.

    The manifest is written before any worker starts and replaced atomically after every determinant;
    a progress log is appended by this process only.  allow_censused and _inline (run the tasks in
    this process, so that test patches apply) are for tests only.
    """
    limits = limits_of(limits)
    validate_dets(dets, allow_censused)
    require(isinstance(workers, int) and workers >= 1, "an explicit worker count >= 1 is required")
    out_dir = Path(out_dir)
    existing = [d for d in dets if (out_dir / f"det{d:05d}.json.gz").exists()]
    require(resume or not existing, f"chunks already exist for {existing}; use --resume to validate and reuse them")
    stamp = time.strftime("%Y%m%d-%H%M%S")
    manifest_path = out_dir / f"manifest-{stamp}.json"
    log_path = out_dir / f"run-{stamp}.log"
    manifest = {"schema": SCHEMA + "-run-manifest", "started": stamp, "status": "running",
                "requested": sorted(dets), "workers": workers, "ladder": list(ladder), "limits": limits,
                "note": "the same limits are used to build and to replay each chunk; mass verification is "
                        "separate and asynchronous (task audit)",
                "resume": resume, "provenance": provenance(), "progress_log": str(log_path).replace("\\", "/"),
                "inputs": {"existing_chunks_sha256": {d: sha256_file(out_dir / f"det{d:05d}.json.gz") for d in existing},
                           "checkpoints_sha256": {d: sha256_file(checkpoint_path(out_dir, d)) for d in dets
                                                  if checkpoint_path(out_dir, d).exists()}},
                "chunks": []}
    start = time.perf_counter()

    def log(line):
        out_dir.mkdir(parents=True, exist_ok=True)
        with open(log_path, "a", encoding="utf-8") as f:
            f.write(time.strftime("%Y-%m-%d %H:%M:%S") + " " + line + "\n")

    def save(final=False):
        rows = sorted(manifest["chunks"], key=lambda r: r["determinant"])
        ok = [r for r in rows if r.get("locally_verified")]
        manifest["chunks"] = rows
        by = lambda o: [r["determinant"] for r in rows if r["outcome"] == o]      # noqa: E731
        manifest["totals"] = {
            "requested": len(dets), "done": len(rows),
            "built": sum(r["status"] == "built" for r in rows), "reused": sum(r["status"] == "reused" for r in rows),
            "locally_verified": len(ok), "failed": by("failed"), "incomplete": by("incomplete"),
            "unresolved": by("unresolved"),
            **{k: sum(r["replay"][k] for r in ok) for k in ("classes", "indecomposable", "plane_orbits",
                                                             "orbit_decisions_replayed", "survivors", "zeros", "hits",
                                                             "discrepancies")},
            "no_hit_determinants": by("no_hit"),
            "verified_hit_determinants": by("verified_hit"),
            "discrepancy_determinants": by("discrepancy"),
            "candidate_determinants": [r["determinant"] for r in rows if r["events"].get("candidate")],
            "disagreement_determinants": [r["determinant"] for r in rows if r["events"].get("disagreement")],
            "near_miss_determinants": [r["determinant"] for r in rows if r["events"].get("near_miss")],
            "complete": len(ok) == len(dets)}
        manifest["wall_seconds"] = time.perf_counter() - start
        if final:
            manifest["status"] = "finished"
        atomic_write_json(manifest_path, manifest, indent=1)

    save()
    log(f"start: {len(dets)} determinants {min(dets)}..{max(dets)}, {workers} workers, resume={resume}")
    domain = shared_domain(dets) if len(dets) > 1 else None   # one independent candidate scan for the batch
    if domain:
        manifest["shared_candidate_domain"] = domain.info
        save()
        log(f"shared candidate domain: {domain.info['forms_in_buckets']} forms for {len(dets)} determinants, "
            f"{domain.info['seconds']}s")
    tasks = [(d, str(out_dir), limits, list(ladder), resume, domain.restrict(d) if domain else None)
             for d in sorted(dets)]
    with (get_context("spawn").Pool(workers, maxtasksperchild=1) if not _inline else _InlinePool()) as pool:
        for row in pool.imap_unordered(_run_task, tasks):
            manifest["chunks"].append(row)
            save()
            rep = row.get("replay", {})
            line = (f"det {row['determinant']} {row['status']} outcome={row['outcome']} classes={rep.get('classes')} "
                    f"events={ {k: len(v) for k, v in row['events'].items()} } {row['seconds_total']:.0f}s"
                    + (f" error={row['error']}" if row.get("error") else ""))
            log(line)
            print(line, flush=True)
    save(final=True)
    log("finished: " + json.dumps(manifest["totals"]))
    manifest["file"] = str(manifest_path).replace("\\", "/")
    return manifest


class _InlinePool:
    """Test helper: the Pool interface, in this process."""

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def imap_unordered(self, f, tasks):
        return map(f, tasks)


def report_totals(manifest):
    """Print the totals.  Exit status 0 only for a complete local run; a verified hit is a successful
    search result and is reported in its own fields, so wrappers must read the outcome fields."""
    t = manifest["totals"]
    print(json.dumps(t))
    for key, label in (("verified_hit_determinants", "VERIFIED HITS"), ("discrepancy_determinants", "DISCREPANCIES"),
                       ("disagreement_determinants", "PROOF DISAGREEMENTS (investigate)"),
                       ("candidate_determinants", "CANDIDATE EVENTS"), ("near_miss_determinants", "near misses")):
        if t[key]:
            print(f"{label} at determinants {t[key]}; see results events/", flush=True)
    if not t["complete"]:
        print(f"INCOMPLETE: failed {t['failed']}, incomplete {t['incomplete']}, unresolved {t['unresolved']}; "
              f"manifest {manifest.get('file')}", file=sys.stderr, flush=True)
    return 0 if t["complete"] else 1


# ----------------------------------------------------------------------------
# verification receipts and the audit (RANK4-PAIRING-SEARCH-PREFLIGHT-VERIFICATION.md items 1-3)
# ----------------------------------------------------------------------------

def receipt_path(out_dir, d, sha):
    return Path(out_dir) / "verifications" / f"det{d:05d}-{sha[:20]}.json"


def record_verification(out_dir, d, path, rep, limits, by, domain=None):
    """A receipt of one independent replay, keyed by determinant and the chunk's SHA-256; written by
    run, replay, retry and audit --replay, read by audit."""
    sha = sha256_file(path)
    rec = {"schema": SCHEMA + "-verification", "determinant": d, "chunk": str(path).replace("\\", "/"),
           "chunk_sha256": sha, "verdict": rep["verdict"], "by": by, "limits": limits,
           "candidate_domain": ("per-determinant independent scan" if domain is None else
                                {"shared_independent_scan": domain.info}),
           "time": time.strftime("%Y-%m-%d %H:%M:%S"), "provenance": provenance(),
           "summary": {k: rep.get(k) for k in ("classes", "indecomposable", "plane_orbits", "orbit_decisions_replayed",
                                               "survivors", "zeros", "hits", "discrepancies", "unresolved")}}
    atomic_write_json(receipt_path(out_dir, d, sha), rec, indent=1)
    return rec


def replay_and_record(path, limits, by, out_dir=None, domain=None):
    path = Path(path)
    rep = replay_chunk(path, limits["max_nodes"], limits["enum_nodes"], domain)
    record_verification(out_dir or path.parent, rep["determinant"], path, rep, limits, by, domain)
    return rep


def shared_domain(dets):
    """One independent candidate scan for a batch of minus-type determinants (verify_rank4_parity.SharedDomain)."""
    import verify_rank4_parity as verifier
    return verifier.SharedDomain.build(sorted(dets), True, "odd")


def replay_batch(dets, out_dir, limits, by):
    """Replay several chunks with one shared candidate scan; yields (d, report)."""
    dets = sorted(dets)
    domain = shared_domain(dets) if len(dets) > 1 else None
    for d in dets:
        yield d, replay_and_record(Path(out_dir) / f"det{d:05d}.json.gz", limits, by, out_dir,
                                   domain.restrict(d) if domain else None)


def baseline_evidence():
    """The completed minus-type census through 997, from its own stored replays and the Magma mass
    cache; every hash is checked against the current files.  Returns through = 997 or None."""
    import magma_calc as C
    problems, evidence = [], []
    tb = Path("results/rank4-parity/replay-task-b-det201.json")
    m2 = Path("results/rank4-minus-type/replay-minus-det202-501.json")
    m3 = Path("results/rank4-minus-type-ext/even-minus-509-997-replay.json")
    try:
        t = json.loads(tb.read_text(encoding="utf-8"))
        ok = t["verdict"] == "verified" and sha256_file(t["input"].replace("\\", "/")) == t["input_sha256"]
        (evidence if ok else problems).append(f"{tb}: det <= 201 replay verdict {t['verdict']}, input hash "
                                              f"{'matches' if ok else 'does not match'}")
        r = json.loads(m2.read_text(encoding="utf-8"))
        ok = r["verdict"] == "verified" and sha256_file(r["source"]) == r["source_sha256"]
        (evidence if ok else problems).append(f"{m2}: 202..501 replay verdict {r['verdict']}, source hash "
                                              f"{'matches' if ok else 'does not match'}")
        x = json.loads(m3.read_text(encoding="utf-8"))
        bad = [c["determinant"] for c in x["chunks"] if c["verdict"] != "verified" or sha256_file(c["file"]) != c["sha256"]]
        ok = x["all_verified"] and x["chunks_done"] == x["chunks_expected"] == 62 and not bad
        (evidence if ok else problems).append(f"{m3}: 509..997, {x['chunks_done']} chunks, all verified and hashes "
                                              f"{'match' if not bad else 'differ at ' + str(bad)}")
    except (OSError, KeyError, ValueError) as error:
        problems.append(f"baseline replay evidence unreadable: {type(error).__name__}: {error}")
    need = [d for d in range(5, 998, 8)]
    have = set()
    if C.CACHE.exists():
        for f in C.CACHE.glob("mass-*.json"):
            rec = json.loads(f.read_text(encoding="utf-8"))
            if rec.get("ok"):
                have.add(rec["meta"]["det"])
    missing = [d for d in need if d not in have]
    (evidence if not missing else problems).append(
        f"Magma mass cache: {len(need) - len(missing)} of {len(need)} determinants 5..997 accepted"
        + (f"; missing {missing}" if missing else ""))
    return {"through": 997 if not problems else None, "evidence": evidence, "problems": problems}


def validate_audit_request(dets):
    require(bool(dets), "audit needs a nonempty request")
    require(all(isinstance(d, int) and d > 0 for d in dets), "determinants must be positive integers")
    require(len(set(dets)) == len(dets), "duplicate determinants requested")
    require(all(d % 8 == 5 for d in dets), "determinants must be 5 mod 8")


def audit(dets, out_dir=OUT, replay=False, limits=None):
    """Coverage audit of the admissible interval spanned by the request.

    Scope: every determinant = 5 mod 8 from min(dets) to max(dets), so a determinant missing from a
    sparse request is still examined and stops the frontier.  Per determinant: locally verified (a
    receipt or a manifest replay for this exact chunk hash; with replay=True a fresh replay, which
    writes a receipt) and mass verified (an accepted cached mass job whose class mass equals this
    chunk's).  contiguous_audited_through: the largest d such that every admissible determinant from
    the interval's start to d has both.  global_complete_through: that value joined to the census
    baseline through 997 (checked by baseline_evidence), only if the interval starts at 1005; otherwise
    null, and the result is coverage within the interval only.  complete: the whole interval is
    audited.  No network access.
    """
    import magma_calc as C
    from fractions import Fraction as F
    validate_audit_request(dets)
    lim = limits_of(limits)
    out_dir = Path(out_dir)
    lo, hi = min(dets), max(dets)
    interval = [d for d in range(lo, hi + 1) if d % 8 == 5]
    verified = set()
    for f in sorted((out_dir / "verifications").glob("det*.json")):
        try:
            r = json.loads(f.read_text(encoding="utf-8"))
        except ValueError:
            continue
        if r.get("verdict") == "verified":
            verified.add((r["determinant"], r["chunk_sha256"]))
    for mf in sorted(out_dir.glob("manifest-*.json")):
        try:
            m = json.loads(mf.read_text(encoding="utf-8"))
        except ValueError:
            continue
        for r in m.get("chunks", []):
            if r.get("sha256") and (r.get("locally_verified") or r.get("replay", {}).get("verdict") == "verified"):
                verified.add((r["determinant"], r["sha256"]))
    mass_cache = {}
    if C.CACHE.exists():
        for f in C.CACHE.glob("mass-*.json"):
            rec = json.loads(f.read_text(encoding="utf-8"))
            if rec.get("ok"):
                mass_cache.setdefault(rec["meta"]["det"], []).append(rec)
    rows = []
    domain = None
    if replay:
        todo = [d for d in interval if (out_dir / f"det{d:05d}.json.gz").exists()
                and (d, sha256_file(out_dir / f"det{d:05d}.json.gz")) not in verified]
        domain = shared_domain(todo) if len(todo) > 1 else None
    for d in interval:
        path = out_dir / f"det{d:05d}.json.gz"
        rec = {"determinant": d, "requested": d in dets, "chunk": str(path).replace("\\", "/")}
        if not path.exists():
            rows.append(rec | {"locally_verified": False, "mass_verified": False, "reason": "no chunk"})
            continue
        sha = sha256_file(path)
        local = (d, sha) in verified
        if not local and replay:
            local = replay_and_record(path, lim, "audit --replay", out_dir,
                                      domain.restrict(d) if domain else None)["verdict"] == "verified"
        doc = load_doc(path)
        ours = sum((F(1, r["aut_order"]) for r in doc["lattices"]), F(0))
        jobs = [r for r in mass_cache.get(d, [])
                if F(r["meta"]["ours"]) == ours and r["meta"]["classes"] == len(doc["lattices"])]
        rec |= {"chunk_sha256": sha, "locally_verified": local, "class_mass": str(ours),
                "mass_verified": bool(jobs), "mass_job": jobs[0]["id"] if jobs else None,
                "mass_magma": jobs[0].get("magma") if jobs else None}
        if not jobs:
            rec["mass_command"] = f"python -B rank4_pairing_magma.py submit mass --det {d}"
        if not local:
            rec["replay_command"] = f"python -B rank4_pairing_search.py replay --det {d}"
        rows.append(rec)
    through = None
    for rec in rows:
        if rec["locally_verified"] and rec["mass_verified"]:
            through = rec["determinant"]
        else:
            break
    complete = all(r["locally_verified"] and r["mass_verified"] for r in rows)
    base = baseline_evidence()
    joined = through if (base["through"] == 997 and lo == 1005 and through is not None) else None
    doc = {"schema": SCHEMA + "-audit", "time": time.strftime("%Y-%m-%d %H:%M:%S"),
           "requested": sorted(dets), "interval": [lo, hi], "admissible_in_interval": len(interval),
           "not_requested_but_in_interval": [d for d in interval if d not in dets],
           "complete": complete, "contiguous_audited_through": through,
           "scope_of_contiguous_value": f"within the interval {lo}..{hi} only",
           "global_complete_through": joined,
           "global_note": ("joined to the census baseline through 997" if joined is not None else
                           "not a global claim: needs the interval to start at 1005 and the baseline evidence"),
           "baseline": base, "rows": rows, "provenance": provenance()}
    name = f"audit-{time.strftime('%Y%m%d-%H%M%S')}-{sha256_json(doc)[:10]}.json"
    atomic_write_json(out_dir / name, doc, indent=1)
    doc["file"] = str(out_dir / name).replace("\\", "/")
    return doc


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("task", choices=("controls", "regress", "benchmark", "run", "replay", "retry", "audit"))
    parser.add_argument("--det", type=int, nargs="*", default=[])
    parser.add_argument("--dmin", type=int)
    parser.add_argument("--dmax", type=int)
    parser.add_argument("--workers", type=int, default=None, help="required for run (explicit worker count)")
    parser.add_argument("--go", action="store_true", help="required for the search on new determinants")
    parser.add_argument("--resume", action="store_true", help="validate (replay) and reuse existing chunks")
    parser.add_argument("--replay", action="store_true", help="audit: replay chunks with no recorded verification")
    parser.add_argument("--out", default=str(OUT), help="chunk directory")
    for name, value in DEFAULT_LIMITS.items():
        parser.add_argument("--" + name.replace("_", "-"), type=int, default=value)
    parser.add_argument("--workdir", default=None)
    args = parser.parse_args()
    limits = limits_of({k: getattr(args, k) for k in DEFAULT_LIMITS})

    def requested():
        if args.det:
            require(args.dmin is None and args.dmax is None, "give either --det or --dmin/--dmax, not both")
            return list(args.det)
        require(args.dmin is not None and args.dmax is not None, "give --det or both --dmin and --dmax")
        require(args.dmin <= args.dmax, "--dmin must not exceed --dmax")
        return [d for d in range(args.dmin, args.dmax + 1) if d % 8 == 5]
    status = 0
    if args.task == "controls":
        for f in (control_det21, control_shiota, control_det1729):
            print(json.dumps(f(), default=str))
    elif args.task == "regress":
        print(json.dumps(regress_all(args.workers or 16), indent=1))
    elif args.task == "benchmark":
        wd = Path(args.workdir) if args.workdir else None
        require(wd is not None, "--workdir (scratch) required")
        for d in args.det:
            print(json.dumps(benchmark(d, wd)))
    elif args.task == "run":
        require(args.go, "the search was not started: pass --go to run it")
        dets = requested()
        validate_dets(dets)
        require(args.workers is not None and args.workers >= 1, "--workers must be given explicitly (>= 1)")
        status = report_totals(run_search(dets, args.workers, Path(args.out), args.resume, limits))
    elif args.task == "audit":
        doc = audit(requested(), Path(args.out), args.replay, limits)
        print(json.dumps({k: doc[k] for k in ("interval", "complete", "contiguous_audited_through",
                                              "global_complete_through", "not_requested_but_in_interval", "file")}
                         | {"rows": [{k: r.get(k) for k in ("determinant", "locally_verified", "mass_verified")}
                                     for r in doc["rows"]]}, indent=1))
        status = 0 if doc["complete"] else 1        # a coverage gate: success only if the interval is audited
    elif args.task == "retry":
        for d in args.det:
            print(json.dumps(retry_unresolved(d, Path(args.out), limits)))
            rep = replay_and_record(Path(args.out) / f"det{d:05d}.json.gz", limits, "retry")
            print(json.dumps(rep))
            status |= rep["verdict"] != "verified"
    else:
        for d, rep in replay_batch(requested(), Path(args.out), limits, "replay"):
            print(json.dumps(rep))
            status |= rep["verdict"] != "verified"
    sys.exit(status)


if __name__ == "__main__":
    main()
