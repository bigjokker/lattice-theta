"""Independent checks of the six-cell classification (analysis only).

  recheck   every saved cell: the generic point really has six balanced pairs (exact minima and multiplicities
            recomputed from scratch with an independent box enumeration), T is integral, unimodular, T = tau mod 2,
            and T preserves every saved generator; statistics of T (order, T^2).
  fan2      recompute the secondary fan from a different generic start and compare the three central forms up to
            isometry with the saved ones.
  controls  the certification step is not vacuous: (a) minimum equality only (multiplicities ignored) and
            (b) transvection markings tau = I + N with rank N = 1 must produce uncertified cells.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import json
from multiprocessing import get_context
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import rank4_pairing_search as S            # noqa: E402
import rank4_secondary_fan as SF            # noqa: E402
import rank4_six_cells as C                 # noqa: E402
from verify_lattice_census import box_vectors   # noqa: E402

HERE = Path(__file__).resolve().parent


def box_minima(Q):
    """Independent: class minima and multiplicities from one complete box (verify_lattice_census.box_vectors)."""
    R = max(SF.qform(Q, [((k >> i) & 1) for i in range(4)]) for k in range(1, 16))
    best = {}
    for v in box_vectors(Q, R, 50_000_000):
        if not any(v):
            continue
        k = C.label(v)
        n = SF.qform(Q, v)
        if k not in best or n < best[k][0]:
            best[k] = [n, 1]
        elif n == best[k][0]:
            best[k][1] += 1
    return best


def recheck_one(cell):
    Q = cell["generic_point"]
    tau = cell["tau"]
    mins = box_minima(Q)
    balanced = all(mins[a] == mins[b] for a, b in C.pairs_of(tau))
    T = cell["T"]
    ok_T = (T is not None and abs(S.determinant(T)) == 1 and list(S.mod2_perm(T)) == tau and
            all([[sum(T[a][i] * G[a][b] * T[b][j] for a in range(4) for b in range(4)) for j in range(4)]
                 for i in range(4)] == G for G in cell["generators"] + [Q]))
    T2 = [[sum(T[i][k] * T[k][j] for k in range(4)) for j in range(4)] for i in range(4)] if T else None
    I = [[int(i == j) for j in range(4)] for i in range(4)]
    kind = "T^2=I" if T2 == I else "T^2=-I" if T2 == [[-x for x in r] for r in I] else "other"
    return {"balanced": balanced, "T_ok": ok_T, "kind": kind, "dim": cell["cell_dim"]}


def recheck():
    doc = json.loads((ROOT / C.OUT / "cells.json").read_text())
    with get_context("spawn").Pool(12) as pool:
        res = pool.map(recheck_one, doc["cells"], chunksize=16)
    out = {"cells": len(res), "generic_point_six_balanced": sum(r["balanced"] for r in res),
           "T_verified": sum(r["T_ok"] for r in res), "T_kinds": dict(Counter(r["kind"] for r in res)),
           "kinds_by_dim": {f"dim {d} {k}": n for (d, k), n in sorted(Counter((r["dim"], r["kind"]) for r in res).items())}}
    print(json.dumps(out))
    return out


def fan2():
    start = [[7, 2, -3, 1], [2, 9, 1, -4], [-3, 1, 8, 2], [1, -4, 2, 11]]
    reps = SF.fan(start=start)
    saved = json.loads((ROOT / SF.OUT / "fan.json").read_text())["cones"]
    match = [[SF.isometric(r["central"], s["central"]) for s in saved] for r in reps]
    out = {"orbits": len(reps), "each_matches_exactly_one_saved": all(sum(m) == 1 for m in match),
           "saved_all_found": all(any(match[i][j] for i in range(len(reps))) for j in range(len(saved)))}
    print(json.dumps(out))
    return out


def transvections():
    out = []
    for bits in range(1 << 16):
        N = [[(bits >> (4 * i + j)) & 1 for j in range(4)] for i in range(4)]
        N2 = [[sum(N[i][k] * N[k][j] for k in range(4)) % 2 for j in range(4)] for i in range(4)]
        if any(any(r) for r in N2):
            continue
        img = {C.label([sum(N[i][j] * S.BITS[k][j] for j in range(4)) for i in range(4)]) for k in range(16)}
        if len(img) != 2:
            continue
        out.append(tuple(k ^ C.label([sum(N[i][j] * S.BITS[k][j] for j in range(4)) for i in range(4)]) for k in range(16)))
    return sorted(set(out))


def control_face(task):
    mode, cone_id, rays, subset = task
    R = [rays[i] for i in subset]
    Q1 = C.form_of(R, [1] * len(R))
    if not SF.positive_definite(Q1):
        return []
    import random
    rng = random.Random(7)
    m1 = C.minimal_sets(Q1)
    taus = C.TAUS if mode == "min_only" else TV
    out = []
    for tau in taus:
        prs = C.pairs_of(tau)
        if mode == "transvection" and any(len(m1[a][1]) != len(m1[b][1]) for a, b in prs):
            continue
        eqs = [[SF.qform(SF.to_Q(r), m1[a][1][0]) - SF.qform(SF.to_Q(r), m1[b][1][0]) for r in R] for a, b in prs]
        cell, _ = C.cell_rays(R, eqs)
        if not cell or not all(any(c[j] > 0 for c in cell) for j in range(len(R))):
            continue
        if mode == "min_only" and all(len(m1[a][1]) == len(m1[b][1]) for a, b in prs):
            continue                                  # only the cells the real classifier excludes
        T, _, _ = C.certify_cell(R, cell, tau, rng, tries=2)
        out.append(T is not None)
    return out


TV = transvections()


def controls():
    fan = json.loads((ROOT / SF.OUT / "fan.json").read_text())
    res = {}
    for mode in ("min_only", "transvection"):
        tasks = [(mode, cid, cone["rays"], tuple(i for i in range(10) if mask >> i & 1))
                 for cid, cone in enumerate(fan["cones"]) for mask in range(1, 1024)]
        with get_context("spawn").Pool(12) as pool:
            flat = [x for r in pool.map(control_face, tasks, chunksize=8) for x in r]
        res[mode] = {"cells": len(flat), "certified": sum(flat), "uncertified": len(flat) - sum(flat)}
    res["transvection_markings"] = len(TV)
    print(json.dumps(res))
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("task", choices=("recheck", "fan2", "controls"))
    a = ap.parse_args()
    t = time.perf_counter()
    r = {"recheck": recheck, "fan2": fan2, "controls": controls}[a.task]()
    (HERE / f"validate-{a.task}.json").write_text(json.dumps(dict(r, seconds=round(time.perf_counter() - t, 1)), indent=1) + "\n")
