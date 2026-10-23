"""Finite classification of six-balanced-minima cells (RANK4-SIX-MINIMA-COMMON-LATTICE.md section 3).

Universe: ALL real positive definite quaternary forms H on Z^4, together with every linear involution tau of
F2^4 with tau = I + N, N^2 = 0, rank N = 2 (210 markings).  For (H, tau) the twelve classes of Z^4/2Z^4 moved
by tau form six pairs {k, tau k}; "six balanced minima" means: for every pair, the two classes have equal
minimal norm AND equal minimal multiplicity.  (For H = half the Gram of the common lattice K of (L, P), this
is exactly "all six sub-coset minimal shells of the plane are balanced".)

Method (exact).  The secondary fan (rank4_secondary_fan.py) has three GL4(Z)-orbits of maximal cones, all
simplicial with 10 rays; every form lies in the closure of a cone of one of the three representatives, and
every face is a ray subset.  The condition is GL4(Z)-equivariant, so it suffices to treat, for each
representative cone, every ray subset S whose relative interior is positive definite, and every tau.
  1. On relint(S) the Delaunay subdivision is constant, hence so is the set M_k of minimal vectors of every
     class (recomputed at two different relative-interior points as a check).
  2. Multiplicities: |M_k| = |M_tau k| for the six pairs, else no balanced point in relint(S).
  3. Minimum norms: H(v_k) = H(v_tau k) is linear in the ray weights mu; the cell is
     {mu > 0 on S, six equalities}; its closure's rays are computed exactly (cdd); relint(S) meets the cell
     iff every ray of S has positive weight on some cell ray.
  4. Symmetry: at a generic point of the cell, every automorphism T with T = tau mod 2 is tested; a cell is
     CERTIFIED if one such T preserves every generator of the cell (hence the whole cell).  A nonempty cell
     with no such T is reported as a geometric counterexample to the symmetry lemma.

Usage: python -B rank4_six_cells.py run [--workers 12]
"""
import argparse
from collections import Counter
import hashlib
from fractions import Fraction as F
from itertools import product
import json
from multiprocessing import get_context
from pathlib import Path
import random
import time

import cdd.gmp as cdd

import rank4_pairing_search as S
import rank4_secondary_fan as SF
import rank4_six_minima as M

OUT = Path("results/rank4-six-cells")
BITS = S.BITS


def label(v):
    return sum((int(x) % 2) << i for i, x in enumerate(v))


def markings():
    """All tau = I + N over F2^4 with N^2 = 0 and rank N = 2, as permutations of the 16 class indices."""
    out = []
    for bits in range(1 << 16):
        N = [[(bits >> (4 * i + j)) & 1 for j in range(4)] for i in range(4)]
        N2 = [[sum(N[i][k] * N[k][j] for k in range(4)) % 2 for j in range(4)] for i in range(4)]
        if any(any(r) for r in N2):
            continue
        image = {label([sum(N[i][j] * BITS[k][j] for j in range(4)) for i in range(4)]) for k in range(16)}
        if len(image) != 4:                       # rank 2
            continue
        perm = [k ^ label([sum(N[i][j] * BITS[k][j] for j in range(4)) for i in range(4)]) for k in range(16)]
        out.append(tuple(perm))
    out = sorted(set(out))
    assert len(out) == 210, len(out)
    return out


def pairs_of(tau):
    return sorted({tuple(sorted((k, tau[k]))) for k in range(16) if tau[k] != k})


def minimal_sets(Q):
    return {k: M.class_minimum(Q, k) for k in range(1, 16)}


def form_of(rays, weights):
    q = [sum(w * r[i] for w, r in zip(weights, rays)) for i in range(10)]
    return SF.to_Q(q)


def cell_rays(rays, eqs):
    """Exact rays of {mu >= 0, eqs . mu = 0} in R^len(rays)."""
    n = len(rays)
    rows = [[0] + [int(i == j) for j in range(n)] for i in range(n)]
    lin = []
    for e in eqs:
        lin.append(len(rows))
        rows.append([0] + list(e))
    m = cdd.matrix_from_array(rows, rep_type=cdd.RepType.INEQUALITY, lin_set=lin)
    gens = cdd.copy_generators(cdd.polyhedron_from_matrix(m))
    out = []
    for row, is_lin in zip(gens.array, range(len(gens.array))):
        if row[0] == 0 and any(row[1:]):
            out.append([F(x) for x in row[1:]])
    lin_gens = getattr(gens, "lin_set", set())
    return out, lin_gens


def certify_cell(rays, cell, tau, rng, tries=4):
    """An automorphism T = tau mod 2 preserving every cell generator, or None."""
    from verify_rank4_feasibility import box_isometries
    gens = [SF.to_Q(SF.primitive([sum(c[j] * rays[j][i] for j in range(len(rays))) for i in range(10)])) for c in cell]
    for _ in range(tries):
        w = [rng.randint(1, 5) for _ in cell]
        mu = [sum(wi * c[j] for wi, c in zip(w, cell)) for j in range(len(rays))]
        Qs = SF.to_Q(SF.primitive([sum(mu[j] * rays[j][i] for j in range(len(rays))) for i in range(10)]))
        for T in box_isometries(Qs, Qs, 50_000_000):
            if tuple(S.mod2_perm(T)) == tau and all(
                    [[sum(T[a][i] * Gm[a][b] * T[b][j] for a in range(4) for b in range(4)) for j in range(4)]
                     for i in range(4)] == Gm for Gm in gens):
                return T, Qs, gens
    return None, Qs, gens


def dze_holds(R, eqs, T):
    """Astra's check: with D[:, j] = entries of T^t R_j T - R_j, D = Z E for a rational Z, i.e. every row of D lies
    in the row space of E.  Then T preserves H(mu) for EVERY mu with E mu = 0 (the whole equality subspace)."""
    D = []
    for r in R:
        Q = SF.to_Q(r)
        TQT = [[sum(T[a][i] * Q[a][b] * T[b][j] for a in range(4) for b in range(4)) for j in range(4)] for i in range(4)]
        D.append([x - y for x, y in zip(SF.to_q(TQT), r)])
    Drows = [[D[j][i] for j in range(len(R))] for i in range(10)]
    return SF_rank(eqs + Drows) == SF_rank(eqs) if eqs else all(not any(row) for row in Drows)


def process_face(task):
    """One face: an explicit disposition for each of the 210 markings, or a consistency failure for the face."""
    cone_id, rays, subset = task
    R = [rays[i] for i in subset]
    Q1 = form_of(R, [1] * len(R))
    if not SF.positive_definite(Q1):
        return {"cone": cone_id, "subset": subset, "pd": False}
    rng = random.Random(hash((cone_id, subset)) & 0xFFFF)
    Q2 = form_of(R, [rng.randint(2, 7) for _ in R])
    m1, m2 = minimal_sets(Q1), minimal_sets(Q2)
    same = all(set(m1[k][1]) == set(m2[k][1]) for k in m1)
    out = {"cone": cone_id, "subset": subset, "pd": True, "minimal_sets_constant": same,
           "multiplicities": {k: len(m1[k][1]) for k in m1}, "cells": [], "dispositions": Counter()}
    if not same:
        out["consistency_failure"] = True            # never skipped silently: the run cannot be proof-complete
        return out
    for tau in TAUS:
        prs = pairs_of(tau)
        if any(len(m1[a][1]) != len(m1[b][1]) for a, b in prs):
            out["dispositions"]["multiplicity_mismatch"] += 1
            continue
        eqs = []
        for a, b in prs:
            va, vb = m1[a][1][0], m1[b][1][0]
            eqs.append([SF.qform(SF.to_Q(r), va) - SF.qform(SF.to_Q(r), vb) for r in R])
        cell, _ = cell_rays(R, eqs)
        if not cell:
            out["dispositions"]["infeasible"] += 1
            continue
        if not all(any(c[j] > 0 for c in cell) for j in range(len(R))):
            out["dispositions"]["lower_face"] += 1     # the cell lies in a proper face: treated there
            continue
        T, Qs, gens = certify_cell(R, cell, tau, rng)
        dze = T is not None and dze_holds(R, eqs, T)
        ok = T is not None and dze
        out["dispositions"]["certified" if ok else "uncertified"] += 1
        dim = SF_rank([[sum(c[j] * R[j][i] for j in range(len(R))) for i in range(10)] for c in cell])
        out["cells"].append({"tau": list(tau), "cell_dim": dim, "cell_rays": len(cell), "generic_point": Qs,
                             "certified": ok, "T": T, "D_equals_ZE": dze, "equalities": eqs, "generators": gens})
    out["dispositions"] = dict(out["dispositions"])
    return out


def SF_rank(vectors):
    """Exact rank of a list of rational vectors (dimension of the cell in form space)."""
    rows = [[F(x) for x in v] for v in vectors]
    r = 0
    for col in range(len(rows[0]) if rows else 0):
        k = next((i for i in range(r, len(rows)) if rows[i][col] != 0), None)
        if k is None:
            continue
        rows[r], rows[k] = rows[k], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][col] != 0:
                f = rows[i][col] / rows[r][col]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[r])]
        r += 1
    return r


TAUS = markings()


def coverage_ok(fan):
    """Every facet of every representative has a coverage record (psd_boundary or certified neighbour)."""
    return all(len(c.get("coverage", [])) == len(c["facets"]) for c in fan["cones"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("task", choices=("run",))
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limit", type=int, default=None, help="test runs only: output marked incomplete")
    ap.add_argument("--fan", default=str(SF.OUT / "fan-exact.json"))
    ap.add_argument("--out", default="cells-v2.json", help="output file name in results/rank4-six-cells (never overwritten)")
    a = ap.parse_args()
    fan_path = Path(a.fan)
    fan = json.loads(fan_path.read_text())
    tasks = []
    for cid, cone in enumerate(fan["cones"]):
        for mask in range(1, 1 << len(cone["rays"])):
            tasks.append((cid, cone["rays"], tuple(i for i in range(len(cone["rays"])) if mask >> i & 1)))
    expected_faces = len(tasks)
    if a.limit:
        tasks = tasks[:a.limit]
    t = time.perf_counter()
    with get_context("spawn").Pool(a.workers) as pool:
        res = list(pool.imap_unordered(process_face, tasks, chunksize=8))
    pd = [r for r in res if r["pd"]]
    disp = Counter()
    for r in pd:
        disp.update(r.get("dispositions", {}))
    failures = [r for r in pd if r.get("consistency_failure")]
    cells = [dict(c, cone=r["cone"], subset=r["subset"]) for r in pd for c in r["cells"]]
    cases = sum(disp.values())
    summary = {"fan": str(fan_path).replace(chr(92), "/"), "fan_sha256": hashlib.sha256(fan_path.read_bytes()).hexdigest(),
               "source_sha256": {f: hashlib.sha256(Path(f).read_bytes()).hexdigest()
                                 for f in ("rank4_six_cells.py", "rank4_secondary_fan.py", "rank4_six_minima.py")},
               "faces_expected": expected_faces, "faces_processed": len(res), "pd_faces": len(pd),
               "markings": len(TAUS), "cases_expected": len(pd) * len(TAUS) if not failures else None,
               "cases_with_disposition": cases, "dispositions": dict(disp),
               "consistency_failures": len(failures), "fan_coverage_certificate": coverage_ok(fan),
               "cells": len(cells), "certified": sum(c["certified"] for c in cells),
               "uncertified": sum(not c["certified"] for c in cells),
               "cells_by_dim": dict(Counter(c["cell_dim"] for c in cells)), "seconds": round(time.perf_counter() - t, 1)}
    summary["proof_complete"] = (len(res) == expected_faces and not failures and summary["uncertified"] == 0
                                 and cases == len(pd) * len(TAUS) and summary["fan_coverage_certificate"]
                                 and all(c["D_equals_ZE"] for c in cells))
    OUT.mkdir(parents=True, exist_ok=True)
    name = a.out if not a.limit else f"cells-partial-limit{a.limit}.json"
    if not a.limit:
        S.require(not (OUT / name).exists(), f"{OUT / name} exists: refusing to overwrite saved classification results")
    summary["complete_run"] = not a.limit
    (OUT / name).write_text(json.dumps({"summary": summary, "cells": cells,
                                        "consistency_failures": [{k: r[k] for k in ("cone", "subset")} for r in failures]},
                                       indent=1, default=str) + "\n", encoding="utf-8")
    print(json.dumps(summary))
    for c in cells:
        if not c["certified"]:
            print("UNCERTIFIED", json.dumps({k: c[k] for k in ("cone", "subset", "tau", "cell_dim", "generic_point")}))
    raise SystemExit(0 if summary["proof_complete"] else 1)


if __name__ == "__main__":
    main()
