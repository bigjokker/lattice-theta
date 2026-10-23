"""Regenerate the first fan run (fan.json, overwritten by mistake on 2026-10-23) with the original algorithm.

The original `fan()` was deterministic: start at A4*, depth-first queue (pop from the end), a cone is new iff its
central form is not isometric to a known one, neighbours across every positive definite facet at the first of
t = 1/8, 1/32, 1/128 that gives a generic form.  It is re-run here with the current (exact) secondary_cone.
The regenerated file is then checked against cells.json: every saved generator and generic point of a cell must be
a nonnegative combination of the rays of its stated face (cone, subset), and the central forms must equal the three
recorded in RANK4-SIX-CELLS-CLASSIFICATION.md.
"""
from fractions import Fraction as F
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))
import rank4_secondary_fan as SF        # noqa: E402

RECORDED = [[[4, -1, -1, -1], [-1, 4, -1, -1], [-1, -1, 4, -1], [-1, -1, -1, 4]],
            [[5, -2, -2, -2], [-2, 6, -1, -1], [-2, -1, 6, -1], [-2, -1, -1, 6]],
            [[6, -2, -3, -3], [-2, 6, -1, -1], [-3, -1, 6, 1], [-3, -1, 1, 6]]]


def fan_v1():
    Q0 = [[4, -1, -1, -1], [-1, 4, -1, -1], [-1, -1, 4, -1], [-1, -1, -1, 4]]
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
        print(f"new cone orbit {len(reps)}: central {c}", flush=True)
        cq = SF.to_q(c)
        for fct in cone["facets"]:
            on = [r for r in cone["rays"] if sum(a * b for a, b in zip(fct, r)) == 0]
            fpt = [sum(r[k] for r in on) for k in range(10)]
            if not SF.positive_definite(SF.to_Q(fpt)):
                continue
            for t in (F(1, 8), F(1, 32), F(1, 128)):
                Qn = SF.to_Q(SF.primitive([x + t * (x - y) for x, y in zip(fpt, cq)]))
                if SF.positive_definite(Qn) and SF._generic(Qn):
                    queue.append(Qn)
                    break
            else:
                raise RuntimeError("no generic point across a facet")
    return reps


def in_cone(rays, Q):
    """Is the form Q a nonnegative combination of the given rays? (exact LP via cdd)."""
    import cdd.gmp as cdd
    q = SF.to_q(Q)
    n = len(rays)
    rows = [[0] + [int(i == j) for j in range(n)] for i in range(n)]       # mu >= 0
    lin = []
    for i in range(10):                                                      # sum mu_j r_j = q
        lin.append(len(rows))
        rows.append([-q[i]] + [r[i] for r in rays])
    m = cdd.matrix_from_array(rows, rep_type=cdd.RepType.INEQUALITY, lin_set=lin)
    gens = cdd.copy_generators(cdd.polyhedron_from_matrix(m))
    return any(row[0] == 1 for row in gens.array)                            # a point exists


def main():
    reps = fan_v1()
    centrals = [r["central"] for r in reps]
    doc = {"orbits": len(reps), "cones": reps, "note": __doc__.strip(), "centrals_match_recorded": centrals == RECORDED}
    out = HERE / "fan.json"
    assert not out.exists()
    out.write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
    print("centrals match the recorded ones:", centrals == RECORDED, flush=True)
    cells = json.loads((HERE.parent / "rank4-six-cells" / "cells.json").read_text())["cells"]
    bad = 0
    for c in cells:
        rays = [reps[c["cone"]]["rays"][i] for i in c["subset"]]
        for G in c["generators"] + [c["generic_point"]]:
            if not in_cone(rays, G):
                bad += 1
    print(f"cells.json: {len(cells)} cells, generators/generic points outside their stated face: {bad}", flush=True)


if __name__ == "__main__":
    main()
