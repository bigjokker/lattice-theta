"""Match every cone of the saved fan files to Voronoi's cones Delta, Delta', Delta'' (arXiv:2407.03513 section 2.1).

Published rays are taken from Grok's audit script bot_work/grok/proof_audit/match_primitive_rays.py (R11 replaced
by its primitive representative, as noted there).  For each saved cone, an integral U with det +-1 is searched
among the isometries between the central forms (sum of primitive rays); the match is accepted only if U^t R U runs
EXACTLY over the published cone's primitive rays as R runs over the saved cone's rays.
"""
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "bot_work" / "grok" / "proof_audit"))
import rank4_secondary_fan as SF                 # noqa: E402
import rank4_pairing_search as S                 # noqa: E402
from verify_rank4_feasibility import box_isometries   # noqa: E402
import match_primitive_rays as G                 # noqa: E402  (data only: PUB, CONES)


def published():
    out = {}
    for name, labels in G.CONES.items():
        rays = sorted(tuple(SF.primitive(SF.to_q(G.PUB[l]))) for l in labels)
        out[name] = {"rays": rays, "central": SF.to_Q(SF.primitive([sum(r[k] for r in rays) for k in range(10)]))}
    return out


def match(cone, pubs):
    rays = cone["rays"]
    c = SF.to_Q(SF.primitive([sum(r[k] for r in rays) for k in range(10)]))
    for name, p in pubs.items():
        for U in box_isometries(p["central"], c, 50_000_000):        # U^t c U = published central
            if sorted(tuple(SF.primitive(SF.act(U, r))) for r in rays) == p["rays"]:
                return name, U, S.determinant(U)
    return None


def main():
    pubs = published()
    report = {}
    for f in ("fan.json", "fan-exact.json"):
        fan = json.loads((HERE / f).read_text())
        rows = []
        for i, cone in enumerate(fan["cones"]):
            m = match(cone, pubs)
            rows.append({"cone": i, "matched": m[0] if m else None, "U": m[1] if m else None, "det_U": m[2] if m else None})
            print(f, i, "->", m[0] if m else "NO MATCH", "det U =", m[2] if m else None, flush=True)
        report[f] = rows
    report["all_matched_bijectively"] = all(
        sorted(r["matched"] for r in report[f]) == sorted(pubs) for f in ("fan.json", "fan-exact.json"))
    (HERE / "match-published.json").write_text(json.dumps(report, indent=1) + "\n", encoding="utf-8")
    print("each file's three cones match Delta, Delta', Delta'' bijectively:", report["all_matched_bijectively"])


if __name__ == "__main__":
    main()
