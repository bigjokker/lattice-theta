"""Plus-type test of the combined conjecture (docs/rank-four-stage-two-reductions.md section 1).

Combined conjecture: every zero of an even, orthogonally indecomposable,
positive definite rank-four lattice of odd determinant has a symmetry
witness. By Theorem 5 (docs/rank-four-minus-type.md) the zero set is then
exactly the union, over all h in Aut(L) with h^2 = -I, of the six ordered
pairs (a, c) of distinct nonzero classes in the plane (h+I)L/2L.

Tasks
  run     complete class lists, every indecomposable class decided, for even
          lattices with det = 1 mod 8, 281 <= det <= 497 (one chunk per det),
          using rank4_gap_census.build_chunk unchanged.
  replay  independent replay of every preserved chunk (rank4_gap_census.replay_chunk).
  planes  compare, class by class, the set of zero characteristics with the
          union of the h-planes, on every even odd-det class list: Task B
          (det <= 201), the gap-census region even-odd-202-279, the minus-type
          extension (det = 5 mod 8, 202..501) and the plus-type region.
Below 281 the even odd-det universe is already complete (Task B and the gap
census), so it is not regenerated.

Release layout (run from the repository root): the preserved chunks, summary
and replay are under reports/rank4-plus-type/ (--chunks); fresh outputs go to
--out (default scratch/rank4-plus-type/). This is the research-workspace driver
of 2026-10-21 adapted to the release interface of rank4_gap_census.py
(explicit input and output directories); the preserved outputs were produced
by the workspace copy.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import rank4_gap_census as gc
import rank4_parity_search as builder
import verify_rank4_parity as verifier

REGION = "even-plus-281-497"
CHUNKS = Path("reports/rank4-plus-type")
gc.FILTERS["1_mod_8"] = lambda d: d % 8 == 1
builder.DET_FILTERS["1_mod_8"] = lambda d: d % 8 == 1
verifier.DET_FILTERS["1_mod_8"] = lambda d: d % 8 == 1
gc.REGIONS[REGION] = dict(even_only=True, odd_lattices_only=False, det_filter="1_mod_8", dmin=281, dmax=497)

BITS = [tuple((n >> k) & 1 for k in range(4)) for n in range(16)]


def sources(chunks):
    return [
        ("Task B: even, odd det <= 201", ["reports/rank4-parity/task-b-det201.json"]),
        ("gap census: even, odd det 202..279",
         sorted(p.as_posix() for p in Path("reports/rank4-gap-census/even-odd-202-279").glob("det*.json.gz"))),
        ("minus-type extension: even, det = 5 mod 8, 202..501", ["reports/rank4-minus-type/minus-det202-501.json.gz"]),
        ("plus-type region: even, det = 1 mod 8, 281..497",
         sorted(p.as_posix() for p in (Path(chunks) / REGION).glob("det*.json.gz"))),
    ]


def raw_bytes(path):
    """File bytes; a .json present only as .json.gz is read decompressed (same hash)."""
    p = Path(path)
    if p.exists():
        return p.read_bytes()
    return gzip.decompress(Path(str(p) + ".gz").read_bytes())


def load(path):
    raw = raw_bytes(path)
    return json.loads(gzip.decompress(raw) if str(path).endswith(".gz") else raw)


def mod2(v):
    return tuple(x % 2 for x in v)


def matvec(A, v):
    return [sum(A[i][j] * v[j] for j in range(4)) for i in range(4)]


def solve_mod2(G, b):
    """The unique c in F2^4 with G c = b mod 2 (det G odd)."""
    sols = [c for c in BITS if mod2(matvec(G, list(c))) == mod2(b)]
    assert len(sols) == 1
    return sols[0]


def h_plane_pairs(G, auts):
    I = [[int(i == j) for j in range(4)] for i in range(4)]
    neg = [[-x for x in r] for r in I]
    pairs, count = set(), 0
    for h in auts:
        if [[sum(h[i][k] * h[k][j] for k in range(4)) for j in range(4)] for i in range(4)] != neg:
            continue
        count += 1
        hI = [[h[i][j] + I[i][j] for j in range(4)] for i in range(4)]
        plane = {mod2(matvec(hI, list(x))) for x in BITS} - {(0, 0, 0, 0)}
        assert len(plane) == 3, plane
        pairs |= {(u, v) for u in plane for v in plane if u != v}
    return pairs, count


def compare(doc):
    out = dict(classes=0, indecomposable=0, with_sqrt_minus_I=0, zeros=0, predicted=0,
               unexplained_zeros=0, predicted_but_nonzero=0, unresolved=0, cutoff_zeros=0, examples=[])
    for r in doc["lattices"]:
        out["classes"] += 1
        if r["decomposable"]:
            continue
        out["indecomposable"] += 1
        G = r["G"]
        decisions = r["parity"]["decisions"]
        assert len(decisions) == 105
        zeros = {(mod2(d["a"]), solve_mod2(G, d["b"])) for d in decisions if d["verdict"].startswith("zero")}
        out["unresolved"] += sum(d["verdict"] == "unresolved" for d in decisions)
        out["cutoff_zeros"] += sum(d["verdict"] == "zero_by_cutoff" for d in decisions)
        predicted, n = h_plane_pairs(G, r["automorphisms"])
        assert (n > 0) == bool(r["has_square_root_of_minus_identity"])
        out["with_sqrt_minus_I"] += n > 0
        out["zeros"] += len(zeros)
        out["predicted"] += len(predicted)
        bad, missing = zeros - predicted, predicted - zeros
        out["unexplained_zeros"] += len(bad)
        out["predicted_but_nonzero"] += len(missing)
        if bad or missing or zeros:
            out["examples"].append(dict(G=G, det=r["determinant"], type=r["parity"]["decisions"][0]["f2"]["type"],
                                        zeros=len(zeros), predicted=len(predicted), unexplained=len(bad),
                                        predicted_but_nonzero=len(missing), sqrt_minus_I_elements=n))
    return out


def planes(chunks, output):
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    report = {"schema": "work7-rank4-plus-type-planes-v1",
              "test": "per indecomposable class: zero set == union over h in Aut(L), h^2 = -I, of ordered pairs of "
                      "distinct nonzero classes of (h+I)L/2L; c = G^{-1} b mod 2",
              "sources": []}
    for label, paths in sources(chunks):
        if not paths:
            raise FileNotFoundError(f"no chunks for {label}")
        total = None
        for p in paths:
            res = compare(load(p))
            total = res if total is None else {k: total[k] + v for k, v in res.items()}
        report["sources"].append({"label": label, "files": len(paths),
                                  "sha256": {p: hashlib.sha256(raw_bytes(p)).hexdigest() for p in paths}
                                  if len(paths) <= 2 else "per chunk in the region summary", **total})
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes((json.dumps(report, indent=2) + "\n").encode("utf-8"))
    for s in report["sources"]:
        print({k: v for k, v in s.items() if k not in ("sha256", "examples")})
        for e in s["examples"]:
            print("   ", e)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("task", choices=("run", "replay", "planes"))
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--chunks", type=Path, default=CHUNKS, help="preserved chunks, summary and replay")
    parser.add_argument("--out", type=Path, default=Path("scratch/rank4-plus-type"), help="fresh outputs")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.task == "run":
        gc.run([REGION], args.workers, args.out)
    elif args.task == "replay":
        gc.replay([REGION], args.workers, args.chunks, args.out)
    else:
        planes(args.chunks, args.out / "planes.json")


if __name__ == "__main__":
    main()
