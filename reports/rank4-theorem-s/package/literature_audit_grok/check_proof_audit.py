"""Exact checks for PROOF_AUDIT_Grok_20261023.

Read-only against the project. Writes only next to this file.
Does not call results/rank4-secondary-fan/regenerate_v1.py and does not
write results/.

Checks:
  * secondary cones of the three central forms recorded in
    RANK4-SIX-CELLS-CLASSIFICATION.md, against Voronoi's rays R1..R12
    as printed in arXiv:2407.03513, section 2.1;
  * the v1 fan walk (same algorithm as regenerate_v1.fan_v1);
  * GL4(Z) match of fan-exact.json, if that file is present;
  * positive definite faces of the three cones;
  * extreme rays of the closed principal domain
    {q_ij <= 0 for i != j, row sums >= 0};
  * every even-diagonal 4x4 Gram over Z/4Z of odd determinant has
    det == 1 (mod 4);
  * the plane-invariance matrix Gamma against both standard symplectic J;
  * the determinant-21 involution witness;
  * the count of 210 markings tau = I+N with N^2 = 0 and rank 2.
"""
from __future__ import annotations

import json
import sys
import time
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import rank4_secondary_fan as SF  # noqa: E402
from verify_rank4_feasibility import EnumerationLimit, box_isometries  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "check_proof_audit.json"

RECORDED = [
    [[4, -1, -1, -1], [-1, 4, -1, -1], [-1, -1, 4, -1], [-1, -1, -1, 4]],
    [[5, -2, -2, -2], [-2, 6, -1, -1], [-2, -1, 6, -1], [-2, -1, -1, 6]],
    [[6, -2, -3, -3], [-2, 6, -1, -1], [-3, -1, 6, 1], [-3, -1, 1, 6]],
]

# Vallentin, Weissbach, Zimmermann, arXiv:2407.03513, section 2.1.
# Symmetric-matrix entries, project order applied by SF.to_q.
PUBLISHED = {
    "R1": [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
    "R2": [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
    "R3": [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 0]],
    "R4": [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 1]],
    "R5": [[1, -1, 0, 0], [-1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
    "R6": [[1, 0, -1, 0], [0, 0, 0, 0], [-1, 0, 1, 0], [0, 0, 0, 0]],
    "R7": [[1, 0, 0, -1], [0, 0, 0, 0], [0, 0, 0, 0], [-1, 0, 0, 1]],
    "R8": [[0, 0, 0, 0], [0, 1, -1, 0], [0, -1, 1, 0], [0, 0, 0, 0]],
    "R9": [[0, 0, 0, 0], [0, 1, 0, -1], [0, 0, 0, 0], [0, -1, 0, 1]],
    "R10": [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 1, -1], [0, 0, -1, 1]],
    "R11": [[4, 2, -2, -2], [2, 4, -2, -2], [-2, -2, 4, 0], [-2, -2, 0, 4]],
    "R12": [[1, 1, -1, -1], [1, 1, -1, -1], [-1, -1, 1, 1], [-1, -1, 1, 1]],
}
CONES_PUB = {
    "Delta": [f"R{i}" for i in range(1, 11)],
    "Delta_prime": [f"R{i}" for i in range(1, 5)] + [f"R{i}" for i in range(6, 12)],
    "Delta_double_prime": [f"R{i}" for i in range(1, 5)] + [f"R{i}" for i in range(6, 10)] + ["R11", "R12"],
}


def log(msg):
    print(msg, flush=True)


def key(q):
    return tuple(SF.primitive(q))


def mat_det_int(M):
    return int(SF.det(M))


def mul(A, B):
    n = len(A)
    m = len(B[0])
    p = len(B)
    return [[sum(A[i][k] * B[k][j] for k in range(p)) for j in range(m)] for i in range(n)]


def transpose(A):
    return [list(row) for row in zip(*A)]


def mat_eq(A, B):
    return all(a == b for ra, rb in zip(A, B) for a, b in zip(ra, rb))


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def psd(Q):
    for k in range(1, 5):
        for idx in combinations(range(4), k):
            sub = [[Q[i][j] for j in idx] for i in idx]
            if SF.det(sub) < 0:
                return False
    return True


def ray_matrix(names):
    rows = [SF.to_q(PUBLISHED[name]) for name in names]
    return rows


def independent(rows):
    return SF.det(rows) != 0


def pd_face_count(ray_rows):
    """Nonempty subsets whose ray sum is positive definite."""
    n = len(ray_rows)
    total = 0
    by_dim = [0] * (n + 1)
    for k in range(1, n + 1):
        for subset in combinations(range(n), k):
            acc = [0] * 10
            for i in subset:
                for t, v in enumerate(ray_rows[i]):
                    acc[t] += v
            if SF.positive_definite(SF.to_Q(acc)):
                total += 1
                by_dim[k] += 1
    return total, by_dim


def match_ray_sets(our_rows, pub_names):
    """Return U with U^t (our form) U running over the published rays, or None.

    box_isometries(G, H) returns T with T^t H T = G.  We want
    U^t our_central U = published_central, so G is the published central
    and H is our central.
    """
    pub_rows = [key(SF.to_q(PUBLISHED[name])) for name in pub_names]
    pub_set = set(pub_rows)
    our_keys = [key(r) for r in our_rows]
    if set(our_keys) == pub_set:
        names = []
        for r in our_keys:
            names.append(pub_names[pub_rows.index(r)])
        return {"U": eye(4), "det": 1, "published_ray_of_each_our_ray": names, "direct": True}
    our_c = SF.to_Q(SF.primitive([sum(r[t] for r in our_rows) for t in range(10)]))
    pub_c = SF.to_Q(SF.primitive([sum(SF.to_q(PUBLISHED[n])[t] for n in pub_names) for t in range(10)]))
    try:
        maps = box_isometries(pub_c, our_c, 50_000_000)
    except EnumerationLimit as exc:
        return {"error": f"EnumerationLimit: {exc}", "our_central": our_c, "published_central": pub_c}
    pub_index = {row: name for row, name in zip(pub_rows, pub_names)}
    for U in maps:
        mapped = [key(SF.act(U, r)) for r in our_rows]
        if set(mapped) == pub_set and all(m in pub_index for m in mapped):
            return {
                "U": U,
                "det": mat_det_int(U),
                "published_ray_of_each_our_ray": [pub_index[m] for m in mapped],
                "direct": False,
                "isometries_tested_until_hit": maps.index(U) + 1,
                "isometries_available": len(maps),
            }
    return {
        "matched": False,
        "isometries_available": len(maps),
        "our_central": our_c,
        "published_central": pub_c,
    }


def cone_from_central(Q):
    t = time.perf_counter()
    cone = SF.secondary_cone(Q)
    cone["central"] = SF.central(cone)
    cone["seconds"] = round(time.perf_counter() - t, 3)
    return cone


def fan_v1():
    """Depth-first walk from A4*, as in regenerate_v1.fan_v1. No file write."""
    Q0 = RECORDED[0]
    reps, queue = [], [Q0]
    while queue:
        Q = queue.pop()
        cone = SF.secondary_cone(Q)
        c = SF.central(cone)
        if any(SF.isometric(c, r["central"]) for r in reps):
            continue
        cone["central"] = c
        cone["generic_point"] = Q
        reps.append(cone)
        log(f"fan_v1 new orbit {len(reps)} central {c}")
        cq = SF.to_q(c)
        for fct in cone["facets"]:
            on = [r for r in cone["rays"] if sum(a * b for a, b in zip(fct, r)) == 0]
            fpt = [sum(r[k] for r in on) for k in range(10)]
            if not SF.positive_definite(SF.to_Q(fpt)):
                continue
            placed = False
            for t in (F(1, 8), F(1, 32), F(1, 128)):
                Qn = SF.to_Q(SF.primitive([x + t * (x - y) for x, y in zip(fpt, cq)]))
                if SF.positive_definite(Qn) and SF._generic(Qn):
                    queue.append(Qn)
                    placed = True
                    break
            if not placed:
                raise RuntimeError("no generic point across a facet")
    return reps


def inequality_cone_rays():
    """Closure of {q_ij < 0 (i != j), row sums > 0}, extreme rays via cdd."""
    import cdd.gmp as cdd

    # coordinates (q11,q22,q33,q44,q12,q13,q14,q23,q24,q34)
    ineq = []
    # -q_ij >= 0 for the six off-diagonal slots, indices 4..9
    for i in range(4, 10):
        row = [0] * 11
        row[1 + i] = -1
        ineq.append(row)
    # row sums >= 0
    # q11+q12+q13+q14
    ineq.append([0, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0])
    # q22+q12+q23+q24
    ineq.append([0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0])
    # q33+q13+q23+q34
    ineq.append([0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1])
    # q44+q14+q24+q34
    ineq.append([0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1])
    m = cdd.matrix_from_array(ineq, rep_type=cdd.RepType.INEQUALITY)
    cdd.matrix_canonicalize(m)
    gens = cdd.copy_generators(cdd.polyhedron_from_matrix(m))
    rays = []
    for row in gens.array:
        if row[0] == 0 and any(row[1:]):
            rays.append(SF.primitive(row[1:]))
    return rays


def det_mod4():
    """Even diagonal, arbitrary off-diagonal, entries in Z/4Z."""
    diag_vals = (0, 2)
    off_vals = (0, 1, 2, 3)
    odd = eq1 = eq3 = 0
    counters = []
    for d in product(diag_vals, repeat=4):
        for off in product(off_vals, repeat=6):
            q = [d[0], d[1], d[2], d[3], off[0], off[1], off[2], off[3], off[4], off[5]]
            Q = SF.to_Q(q)
            det = mat_det_int(Q) % 4
            if det % 2 == 0:
                continue
            odd += 1
            if det == 1:
                eq1 += 1
            elif det == 3:
                eq3 += 1
                if len(counters) < 5:
                    counters.append({"Q": Q, "det_mod4": det})
            else:
                raise RuntimeError(f"odd residue {det}")
    return {"matrices": 2 ** 4 * 4 ** 6, "odd_det": odd, "det_cong_1": eq1, "det_cong_3": eq3,
            "counterexamples": counters}


def gamma_check(G):
    D = mat_det_int(G)
    Ginv = SF.rational_inverse(G)
    H = [[int(D * Ginv[i][j]) for j in range(4)] for i in range(4)]
    if any(H[i][j] != H[j][i] for i in range(4) for j in range(4)):
        raise RuntimeError("adjugate not symmetric")
    delta = (3 * D + 1) // 2
    # Gamma = [[2I, 3G], [H, delta I]]
    A = [[2 if i == j else 0 for j in range(4)] for i in range(4)]
    B = [[3 * G[i][j] for j in range(4)] for i in range(4)]
    C = H
    Dm = [[delta if i == j else 0 for j in range(4)] for i in range(4)]
    Gamma = [A[i] + B[i] for i in range(4)] + [C[i] + Dm[i] for i in range(4)]
    # The blocks were concatenated wrong if I add lists. Build properly.
    Gamma = []
    for i in range(4):
        Gamma.append(A[i] + B[i])
    for i in range(4):
        Gamma.append(C[i] + Dm[i])
    I4 = eye(4)
    Z4 = [[0] * 4 for _ in range(4)]
    def block(TL, TR, BL, BR):
        out = []
        for i in range(4):
            out.append(TL[i] + TR[i])
        for i in range(4):
            out.append(BL[i] + BR[i])
        return out
    J_plus = block(Z4, I4, [[-x for x in row] for row in I4], Z4)   # [[0, I], [-I, 0]]
    J_minus = block(Z4, [[-x for x in row] for row in I4], I4, Z4)  # [[0, -I], [I, 0]]
    Gt = transpose(Gamma)
    images = {
        "J_0I_mI0": mul(mul(Gamma, J_plus), Gt),
        "J_0mI_I0": mul(mul(Gamma, J_minus), Gt),
    }
    # Characteristic map on every (a, b) in F2^4 x F2^4.
    # a_new = delta a - H b + diag(delta H)   (integer characteristic)
    # b_new = -3 G a + 2 b + diag(6 G)
    diag_dH = [delta * H[i][i] for i in range(4)]
    diag_6G = [6 * G[i][i] for i in range(4)]
    fails = []
    swaps = 0
    for a in product((0, 1), repeat=4):
        for b in product((0, 1), repeat=4):
            Hb = [sum(H[i][j] * b[j] for j in range(4)) for i in range(4)]
            Ga = [sum(G[i][j] * a[j] for j in range(4)) for i in range(4)]
            a_new = tuple((delta * a[i] - Hb[i] + diag_dH[i]) % 2 for i in range(4))
            b_new = tuple((-3 * Ga[i] + 2 * b[i] + diag_6G[i]) % 2 for i in range(4))
            # c = G^{-1} b mod 2. Since D is odd, c ≡ D G^{-1} b = H b (mod 2).
            c = tuple(x % 2 for x in Hb)
            expect_a = c
            expect_b = tuple(x % 2 for x in Ga)
            if a_new != expect_a or b_new != expect_b:
                if len(fails) < 5:
                    fails.append({"a": a, "b": b, "a_new": a_new, "b_new": b_new,
                                  "expect_a": expect_a, "expect_b": expect_b})
            else:
                swaps += 1
    return {
        "det_G": D,
        "det_G_mod4": D % 4,
        "delta": delta,
        "delta_even": delta % 2 == 0,
        "H": H,
        "det_Gamma": mat_det_int(Gamma),
        "Gamma_J_equals_J": {name: mat_eq(image, J)
                             for (name, image), J in zip(images.items(), (J_plus, J_minus))},
        "block_ATD_minus_CTB_eq_I": mat_eq(
            [[sum(A[k][i] * Dm[k][j] for k in range(4)) - sum(C[k][i] * B[k][j] for k in range(4))
              for j in range(4)] for i in range(4)],
            I4),
        "diag_corrections_even": all(x % 2 == 0 for x in diag_dH + diag_6G),
        "characteristic_pairs_checked": 256,
        "characteristic_swaps_ok": swaps,
        "characteristic_failures": fails,
    }


def det21_witness():
    G = [[2, 1, 0, 0], [1, 2, 0, 0], [0, 0, 2, 1], [0, 0, 1, 4]]
    T = [[-1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
    a = (1, 0, 1, 0)
    c = (0, 1, 0, 1)
    Gc = [sum(G[i][j] * c[j] for j in range(4)) for i in range(4)]
    Ta = [sum(T[i][j] * a[j] for j in range(4)) for i in range(4)]
    Tc = [sum(T[i][j] * c[j] for j in range(4)) for i in range(4)]
    half = [(Ta[i] - a[i]) // 2 for i in range(4)]
    # epsilon = B(half, c) = half^t G c
    eps = sum(half[i] * Gc[i] for i in range(4))
    TtGT = mul(mul(transpose(T), G), T)
    T2 = mul(T, T)
    return {
        "G": G,
        "det": mat_det_int(G),
        "T": T,
        "T_preserves_G": mat_eq(TtGT, G),
        "T_squared": T2,
        "T_fixes_a_mod2": all((Ta[i] - a[i]) % 2 == 0 for i in range(4)),
        "T_fixes_c_mod2": all((Tc[i] - c[i]) % 2 == 0 for i in range(4)),
        "Gc_mod2": [x % 2 for x in Gc],
        "epsilon": eps,
        "epsilon_mod2": eps % 2,
        "B_ac_mod2": sum(a[i] * Gc[i] for i in range(4)) % 2,
    }


def count_markings():
    """Every 4x4 matrix over F2. Count N with N^2 = 0 and rank 2.

    For those, dim ker N = 2 and im N is contained in ker N, so im N = ker N.
    """
    def rank(rows):
        a = [list(r) for r in rows]
        rnk = 0
        used = set()
        for col in range(4):
            piv = next((i for i in range(4) if i not in used and a[i][col] == 1), None)
            if piv is None:
                continue
            used.add(piv)
            rnk += 1
            for i in range(4):
                if i != piv and a[i][col] == 1:
                    a[i] = [(x + y) % 2 for x, y in zip(a[i], a[piv])]
        return rnk

    good = 0
    for entries in product((0, 1), repeat=16):
        N = [list(entries[4 * i:4 * i + 4]) for i in range(4)]
        N2 = [[sum(N[i][k] * N[k][j] for k in range(4)) % 2 for j in range(4)] for i in range(4)]
        if any(N2[i][j] for i in range(4) for j in range(4)):
            continue
        if rank(N) == 2:
            good += 1
    return {"matrices": 2 ** 16, "markings_N2_0_rank2": good}


def summarise_cone(cone, label):
    rays = cone["rays"]
    return {
        "label": label,
        "n_rays": len(rays),
        "n_facets": len(cone["facets"]),
        "central": cone["central"],
        "rays_independent": independent(rays) if len(rays) == 10 else False,
        "ray_det": str(SF.det(rays)) if len(rays) == 10 else None,
        "rays_psd": [psd(SF.to_Q(r)) for r in rays],
        "seconds": cone.get("seconds"),
    }


def main():
    t0 = time.perf_counter()
    doc = {}

    log("published ray sanity")
    pub_q = {name: SF.to_q(Q) for name, Q in PUBLISHED.items()}
    doc["published_rays_psd"] = {name: psd(Q) for name, Q in PUBLISHED.items()}
    doc["published_cones_independent"] = {
        name: bool(independent([pub_q[n] for n in names])) for name, names in CONES_PUB.items()
    }
    doc["sum_R1_to_R10"] = SF.to_Q(SF.primitive(
        [sum(pub_q[f"R{i}"][t] for i in range(1, 11)) for t in range(10)]))
    doc["sum_R1_to_R10_equals_recorded_0"] = doc["sum_R1_to_R10"] == RECORDED[0]

    log("principal inequality cone")
    ineq_rays = inequality_cone_rays()
    target = {key(pub_q[f"R{i}"]) for i in range(1, 11)}
    doc["inequality_cone"] = {
        "n_rays": len(ineq_rays),
        "equals_R1_through_R10": set(map(tuple, ineq_rays)) == target,
        "rays": ineq_rays,
    }

    log("det mod 4")
    doc["det_mod4"] = det_mod4()
    log(doc["det_mod4"])

    log("Gamma on two even odd-det Grams")
    A4 = [[2, -1, 0, 0], [-1, 2, -1, 0], [0, -1, 2, -1], [0, 0, -1, 2]]
    det21 = [[2, 1, 0, 0], [1, 2, 0, 0], [0, 0, 2, 1], [0, 0, 1, 4]]
    doc["gamma"] = {"A4_cartan": gamma_check(A4), "det21": gamma_check(det21)}
    log({k: {kk: vv for kk, vv in g.items() if kk != "H"} for k, g in doc["gamma"].items()})

    log("det 21 witness")
    doc["det21_witness"] = det21_witness()

    log("210 markings")
    doc["markings"] = count_markings()
    log(doc["markings"])

    log("secondary cones of the three recorded centrals")
    recorded_cones = []
    for i, Q in enumerate(RECORDED):
        log(f"secondary_cone recorded {i}")
        cone = cone_from_central(Q)
        recorded_cones.append(cone)
        log(summarise_cone(cone, f"recorded_{i}"))
    doc["recorded"] = []
    for i, cone in enumerate(recorded_cones):
        item = summarise_cone(cone, f"recorded_{i}")
        item["central_equals_input"] = cone["central"] == RECORDED[i]
        item["rays"] = cone["rays"]
        matches = {}
        for name, names in CONES_PUB.items():
            log(f"match recorded {i} to {name}")
            matches[name] = match_ray_sets(cone["rays"], names)
        item["match"] = matches
        n_pd, by_dim = pd_face_count(cone["rays"])
        item["pd_faces"] = n_pd
        item["pd_faces_by_number_of_rays"] = by_dim
        item["nonempty_ray_subsets"] = 2 ** len(cone["rays"]) - 1
        doc["recorded"].append(item)
        log(f"recorded {i} pd_faces {n_pd} matches "
            f"{ {k: ('U' in v) for k, v in matches.items()} }")

    doc["recorded_pd_faces_with_repetition"] = sum(it["pd_faces"] for it in doc["recorded"])
    doc["recorded_ray_subsets_with_repetition"] = sum(it["nonempty_ray_subsets"] for it in doc["recorded"])

    log("fan_v1 walk")
    t = time.perf_counter()
    reps = fan_v1()
    doc["fan_v1"] = {
        "orbits": len(reps),
        "centrals": [r["central"] for r in reps],
        "centrals_equal_recorded": [r["central"] for r in reps] == RECORDED,
        "seconds": round(time.perf_counter() - t, 3),
        "ray_sets_equal_recorded_cones": [
            [key(r) for r in rep["rays"]] == [key(x) for x in recorded_cones[i]["rays"]]
            for i, rep in enumerate(reps)
        ] if len(reps) == 3 else None,
    }
    log(doc["fan_v1"])

    exact_path = ROOT / "results" / "rank4-secondary-fan" / "fan-exact.json"
    if exact_path.exists():
        log("fan-exact.json match")
        exact = json.loads(exact_path.read_text(encoding="utf-8"))
        doc["fan_exact"] = {"orbits": exact.get("orbits"), "matches": []}
        for i, cone in enumerate(exact["cones"]):
            rows = cone["rays"]
            hit = None
            for j, rec in enumerate(recorded_cones):
                m = match_ray_sets(rows, [f"placeholder"]) if False else None
                # Match against the recorded cone by isometry of centrals, then ray transport.
                our_c = SF.to_Q(SF.primitive([sum(r[t] for r in rows) for t in range(10)]))
                # Prefer the stored central if present.
                if "central" in cone:
                    our_c = cone["central"]
                try:
                    maps = box_isometries(RECORDED[j], our_c, 50_000_000)
                except EnumerationLimit as exc:
                    doc["fan_exact"]["matches"].append({"exact_index": i, "recorded_index": j,
                                                       "error": str(exc)})
                    continue
                rec_keys = [key(r) for r in rec["rays"]]
                rec_set = set(rec_keys)
                for U in maps:
                    mapped = [key(SF.act(U, r)) for r in rows]
                    if set(mapped) == rec_set:
                        hit = {"exact_index": i, "recorded_index": j, "U": U,
                               "det": mat_det_int(U), "central": our_c}
                        break
                if hit:
                    break
            doc["fan_exact"]["matches"].append(hit if hit else {"exact_index": i, "matched": False,
                                                                "central": cone.get("central")})
            log(doc["fan_exact"]["matches"][-1])
    else:
        doc["fan_exact"] = {"present": False}

    # Drop bulky rays from the JSON only after matches are stored.
    # Keep rays: the report quotes maps, and the file is the certificate.
    doc["seconds"] = round(time.perf_counter() - t0, 3)
    OUT.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    log(f"wrote {OUT} in {doc['seconds']}s")


if __name__ == "__main__":
    main()
