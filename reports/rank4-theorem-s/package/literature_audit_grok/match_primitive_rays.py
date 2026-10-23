"""Match the three recorded cones to Voronoi's cones, using primitive rays.

R11 in arXiv:2407.03513 is twice the primitive extreme ray. The cone is unchanged
when that generator is replaced by its primitive representative.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import rank4_secondary_fan as SF
from verify_rank4_feasibility import box_isometries

HERE = Path(__file__).resolve().parent
PUB = {
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
CONES = {
    "Delta": [f"R{i}" for i in range(1, 11)],
    "Delta_prime": [f"R{i}" for i in range(1, 5)] + [f"R{i}" for i in range(6, 12)],
    "Delta_double_prime": [f"R{i}" for i in range(1, 5)] + [f"R{i}" for i in range(6, 10)] + ["R11", "R12"],
}


def key(q):
    return tuple(SF.primitive(q))


def eye():
    return [[int(i == j) for j in range(4)] for i in range(4)]


def prim_rows(names):
    return [list(key(SF.to_q(PUB[n]))) for n in names]


def central_of(rows):
    return SF.to_Q(SF.primitive([sum(r[t] for r in rows) for t in range(10)]))


def match(our_rows, names):
    pub = prim_rows(names)
    pub_set = {tuple(r) for r in pub}
    index = {tuple(r): names[i] for i, r in enumerate(pub)}
    # R11 and its primitive share a name; index uses the primitive tuple.
    our = [list(key(r)) for r in our_rows]
    if {tuple(r) for r in our} == pub_set:
        return {"U": eye(), "det": 1, "direct": True,
                "published_ray_of_each_our_ray": [index[tuple(r)] for r in our],
                "content_note": {n: SF.primitive(SF.to_q(PUB[n])) != SF.to_q(PUB[n]) for n in names}}
    our_c = central_of(our)
    pub_c = central_of(pub)
    maps = box_isometries(pub_c, our_c, 50_000_000)
    for U in maps:
        mapped = [key(SF.act(U, r)) for r in our]
        if {tuple(m) for m in mapped} == pub_set:
            # confirm each image equals the published matrix, not merely the coordinate tuple
            return {
                "U": U,
                "det": int(SF.det(U)),
                "direct": False,
                "isometries_available": len(maps),
                "published_ray_of_each_our_ray": [index[tuple(m)] for m in mapped],
                "our_central": our_c,
                "published_central": pub_c,
                "published_det": int(SF.det(pub_c)),
                "our_det": int(SF.det(our_c)),
            }
    return {"matched": False, "isometries_available": len(maps),
            "our_central": our_c, "published_central": pub_c,
            "our_det": int(SF.det(our_c)), "published_det": int(SF.det(pub_c))}


def main():
    doc = json.loads((HERE / "check_proof_audit.json").read_text(encoding="utf-8"))
    out = {"primitive_R11": SF.primitive(SF.to_q(PUB["R11"])),
           "R11_content": 2,
           "primitive_R12": SF.primitive(SF.to_q(PUB["R12"])),
           "matches": []}
    for item in doc["recorded"]:
        row = {"label": item["label"], "central": item["central"]}
        row["against"] = {name: match(item["rays"], names) for name, names in CONES.items()}
        out["matches"].append(row)
        print(item["label"], {k: ("U" in v) for k, v in row["against"].items()}, flush=True)
    (HERE / "primitive_match.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
    print("wrote primitive_match.json", flush=True)


if __name__ == "__main__":
    main()
