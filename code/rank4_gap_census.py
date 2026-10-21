"""Closing the determinant gap 25 <= det <= 279 for the rank-four symmetry converse.

Universe: every isometry class of positive definite integral rank-four
lattices with 25 <= det <= 279 that is not already in a complete decided
list of release 1.1.0 (Task E: odd lattices, 4 !| det, 25..168; Task B: even
lattices, odd det <= 201; Task C: even lattices, 4 | det <= 64). The four
remaining regions are listed in REGIONS; two further regions census
determinant 280 itself. Decomposable classes are not decided: they satisfy
the converse by the uniform orthogonal-sum corollary
(docs/rank-four-decomposition-and-fibres.md). Even lattices with det = 2 or 3
mod 4 do not exist; `summary` confirms the generator produces none.

Each region is split into one chunk per determinant value. A chunk is a
complete class list (rank4_parity_search.class_list with D0 = dmin = det)
with every indecomposable class decided (decide_class_list), in the standard
schema, so verify_rank4_parity.replay_decided replays it unchanged. A zero
with no symmetry witness is proved by the modular cutoff and given its level
(2 to 4); any such zero is a counterexample to the symmetry converse.

Release layout (run from the repository root): the preserved chunks, region
summaries, replays and closure summary are under reports/rank4-gap-census/,
read by `replay` and `summary`; fresh outputs go to --out (default
scratch/rank4-gap-census/). This is the workspace script of 2026-10-21 with
the input and output directories made explicit and the prior-list paths
relocated; the preserved outputs were produced by the workspace copy.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from multiprocessing import get_context
from pathlib import Path
import time

CHUNKS = Path("reports/rank4-gap-census")
REGIONS = {
    "odd-notdiv4-169-279": dict(even_only=False, odd_lattices_only=True, det_filter="not_divisible_by_4", dmin=169, dmax=279),
    "odd-div4-25-279": dict(even_only=False, odd_lattices_only=True, det_filter="divisible_by_4", dmin=25, dmax=279),
    "even-odd-202-279": dict(even_only=True, odd_lattices_only=False, det_filter="odd", dmin=202, dmax=279),
    "even-div4-65-279": dict(even_only=True, odd_lattices_only=False, det_filter="divisible_by_4", dmin=65, dmax=279),
    # determinant 280 itself (not part of the closure): which classes are counterexamples
    "odd-div4-280": dict(even_only=False, odd_lattices_only=True, det_filter="divisible_by_4", dmin=280, dmax=280),
    "even-div4-280": dict(even_only=True, odd_lattices_only=False, det_filter="divisible_by_4", dmin=280, dmax=280),
}
CLOSURE_REGIONS = ("odd-notdiv4-169-279", "odd-div4-25-279", "even-odd-202-279", "even-div4-65-279")
# The complete lists the closure relies on besides the new regions (release paths;
# a .json.gz beside the path is accepted).  The hashes are of the uncompressed bytes.
PRIOR = {
    "census det <= 24 (all lattices)": "data/census-rank4/det24-bound32.json",
    "Task E (odd lattices, 4 !| det, 25..168)": "reports/rank4-parity/task-e-det168.json",
    "Task B (even lattices, odd det <= 201)": "reports/rank4-parity/task-b-det201.json",
    "Task C (even lattices, 4 | det <= 64)": "reports/rank4-parity/task-c-det64.json",
}
FILTERS = {"odd": lambda d: d % 2 == 1, "divisible_by_4": lambda d: d % 4 == 0, "not_divisible_by_4": lambda d: d % 4 != 0}


def admitted(region):
    spec = REGIONS[region]
    return [d for d in range(spec["dmin"], spec["dmax"] + 1) if FILTERS[spec["det_filter"]](d)]


def chunk_path(base, region, det):
    return Path(base) / region / f"det{det:03d}.json.gz"


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        return json.load(f)


def raw_bytes(path):
    """Bytes of a preserved file, expanded if only the .json.gz is present."""
    path = Path(path)
    if path.exists():
        return path.read_bytes()
    return gzip.decompress(Path(str(path) + ".gz").read_bytes())


def chunk_row(region, det, doc, path, seconds):
    s = doc["summary"]
    indec = [r for r in doc["lattices"] if not r["decomposable"]]
    return {"region": region, "determinant": det, "file": str(path).replace("\\", "/"), "sha256": sha256(path),
            "candidates": doc["classes_doc"]["candidates"], "classes": s["classes"],
            "decomposable_classes": s["decomposable_classes"], "indecomposable_classes": s["indecomposable_classes"],
            "P_holds": s["P_holds"], "P_fails": s["P_fails"], "P_unresolved": s["P_unresolved"],
            # keys are ints in memory and strings after a JSON round trip; normalize before summing
            "zero_levels": {str(k): v for k, v in s["zero_levels"].items()},
            "even_zeros": sum(r["parity"]["even_zeros"] for r in indec),
            "zeros_without_witness": sum(1 for r in indec for d in r["parity"]["decisions"] if d["verdict"] == "zero_by_cutoff"),
            "unresolved_characteristics": sum(1 for r in indec for d in r["parity"]["decisions"] if d["verdict"] == "unresolved"),
            "largest_first_nonzero_norm4": s["largest_first_nonzero_norm4"],
            "largest_aut_order": max((r["aut_order"] for r in doc["lattices"]), default=0),
            "seconds": seconds}


def build_chunk(task):
    region, det, out = task
    path = chunk_path(out, region, det)
    if path.exists():
        doc = load(path)
        return chunk_row(region, det, doc, path, doc.get("chunk_seconds", 0.0)) | {"cached": True}
    import rank4_parity_search as builder
    spec = REGIONS[region]
    start = time.perf_counter()
    classes = builder.class_list(det, spec["even_only"], spec["det_filter"], dmin=det,
                                 odd_lattices_only=spec["odd_lattices_only"])
    doc = builder.decide_class_list(classes)
    doc["region"] = region
    doc["chunk_seconds"] = time.perf_counter() - start
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    with gzip.open(tmp, "wt", encoding="utf-8") as f:
        json.dump(doc, f)
        f.write("\n")
    tmp.replace(path)
    return chunk_row(region, det, doc, path, doc["chunk_seconds"])


def replay_chunk(task):
    region, det, chunks = task
    import verify_rank4_parity as verifier
    verifier.FAST_BOX[0] = True
    path = chunk_path(chunks, region, det)
    doc = load(path)
    start = time.perf_counter()
    result = verifier.replay_decided(doc, 5_000_000, 50_000_000)
    return {"region": region, "determinant": det, "file": str(path).replace("\\", "/"), "sha256": sha256(path),
            "seconds": time.perf_counter() - start, **result}


def aggregate(rows):
    keys = ("candidates", "classes", "decomposable_classes", "indecomposable_classes", "P_holds", "P_fails",
            "P_unresolved", "even_zeros", "zeros_without_witness", "unresolved_characteristics")
    out = {k: sum(r[k] for r in rows) for k in keys}
    out["zero_levels"] = dict(sum((Counter(r["zero_levels"]) for r in rows), Counter()))
    out["failing_determinants"] = sorted(r["determinant"] for r in rows if r["P_fails"])
    out["largest_first_nonzero_norm4"] = max((r["largest_first_nonzero_norm4"] or 0 for r in rows), default=None)
    out["largest_aut_order"] = max((r["largest_aut_order"] for r in rows), default=0)
    out["cpu_seconds"] = sum(r["seconds"] for r in rows)
    return out


def run_pool(func, tasks, workers, label, sink):
    rows, start = [], time.perf_counter()
    with get_context("spawn").Pool(workers, maxtasksperchild=1) as pool:
        for row in pool.imap_unordered(func, tasks):
            rows.append(row)
            print(f"  {label} {row['region']} det {row['determinant']} {len(rows)}/{len(tasks)} "
                  f"{row.get('verdict', '')} cl {row.get('classes')} indec {row.get('indecomposable_classes')} "
                  f"fails {row.get('P_fails')} unres {row.get('P_unresolved')} "
                  f"{row['seconds']:.0f}s (wall {time.perf_counter() - start:.0f}s)", flush=True)
            sink(rows)
    return rows


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2) + "\n", encoding="utf-8")


def run(regions, workers, out):
    for region in regions:
        spec = REGIONS[region]
        tasks = [(region, d, out) for d in admitted(region)]
        print(f"region {region}: {len(tasks)} determinant chunks, {workers} workers", flush=True)
        start = time.perf_counter()

        def sink(rows, region=region, spec=spec, tasks=tasks, start=start):
            rows = sorted(rows, key=lambda r: r["determinant"])
            write_json(out / f"{region}-summary.json", {
                "schema": "work7-rank4-gap-census-v1-region", "region": region, **spec,
                "universe": f"all isometry classes of positive definite integral rank-four lattices with "
                            f"{spec['dmin']} <= det <= {spec['dmax']}, det {spec['det_filter']}, "
                            + ("even lattices only" if spec["even_only"] else "odd lattices only"),
                "dedup": "exact isometry within each determinant chunk (rank4_parity_search.class_list); "
                         "distinct determinants are distinct classes",
                "decided": "every indecomposable class, all 105 even characteristics (decide_class_list); "
                           "decomposable classes by the uniform orthogonal-sum corollary",
                "chunks_expected": len(tasks), "chunks_done": len(rows), "complete": len(rows) == len(tasks),
                "totals": aggregate(rows), "wall_seconds": time.perf_counter() - start, "chunks": rows})
        run_pool(build_chunk, tasks, workers, "built", sink)


def replay(regions, workers, chunks, out):
    for region in regions:
        tasks = [(region, d, chunks) for d in admitted(region)]
        missing = [d for _, d, _ in tasks if not chunk_path(chunks, region, d).exists()]
        if missing:
            raise FileNotFoundError(f"{region}: chunks missing for determinants {missing}")
        print(f"replay {region}: {len(tasks)} chunks, {workers} workers", flush=True)
        start = time.perf_counter()

        def sink(rows, region=region, tasks=tasks, start=start):
            rows = sorted(rows, key=lambda r: r["determinant"])
            write_json(out / f"{region}-replay.json", {
                "schema": "work7-rank4-gap-census-v1-replay", "region": region,
                "replay": "verify_rank4_parity.replay_decided per chunk (independent loose-box candidate domain, "
                          "box isometries, box shell enumeration, cutoff replay), FAST_BOX",
                "chunks_expected": len(tasks), "chunks_done": len(rows),
                "all_verified": len(rows) == len(tasks) and all(r["verdict"] == "verified" for r in rows),
                "totals": {k: sum(r[k] for r in rows) for k in ("candidates", "classes", "pairwise_searches",
                                                                  "indecomposable_classes", "P_holds", "P_fails", "P_unresolved")},
                "wall_seconds": time.perf_counter() - start, "chunks": rows})
        run_pool(replay_chunk, tasks, workers, "replayed", sink)


def summary(chunks, out):
    """Consistency of the preserved region summaries and replays, the prior-list hashes, and the closure."""
    import rank4_parity_search as builder
    even_all = builder.candidate_grams(279, even_only=True, det_filter="all", dmin=25)
    mod4 = Counter(builder.determinant(G) % 4 for G in even_all)
    assert mod4[2] == 0 and mod4[3] == 0
    preserved = json.loads((Path(chunks) / "closure-25-279.json").read_text(encoding="utf-8"))
    regions = {}
    for region in CLOSURE_REGIONS:
        s = json.loads((Path(chunks) / f"{region}-summary.json").read_text(encoding="utf-8"))
        r = json.loads((Path(chunks) / f"{region}-replay.json").read_text(encoding="utf-8"))
        assert s["complete"] and r["all_verified"], region
        assert [c["sha256"] for c in s["chunks"]] == [c["sha256"] for c in r["chunks"]], region
        assert [c["sha256"] for c in s["chunks"]] == [sha256(chunk_path(chunks, region, c["determinant"])) for c in s["chunks"]], region
        assert s["totals"]["classes"] == r["totals"]["classes"] and s["totals"]["candidates"] == r["totals"]["candidates"]
        regions[region] = {"universe": s["universe"], "totals": s["totals"], "replay_totals": r["totals"],
                           "replay_all_verified": r["all_verified"], "build_wall_seconds": s["wall_seconds"],
                           "replay_wall_seconds": r["wall_seconds"]}
    prior = {}
    for name, path in PRIOR.items():
        digest = hashlib.sha256(raw_bytes(path)).hexdigest()
        recorded = preserved["prior_complete_lists"][name]
        assert recorded["sha256"] == digest, f"prior list hash differs: {name}"
        prior[name] = {**recorded, "file": path, "hash_verified_here": True}
    totals = {k: sum(r["totals"][k] for r in regions.values())
              for k in ("candidates", "classes", "decomposable_classes", "indecomposable_classes", "P_holds", "P_fails",
                        "P_unresolved", "even_zeros", "zeros_without_witness", "unresolved_characteristics")}
    totals["zero_levels"] = dict(sum((Counter(r["totals"]["zero_levels"]) for r in regions.values()), Counter()))
    assert totals == preserved["new_totals"], "new totals differ from the preserved closure"
    closed = (totals["zeros_without_witness"] == 0 and totals["unresolved_characteristics"] == 0
              and all(set(p["zero_levels"] or {}) <= {"1"} and not p["P_unresolved"] for p in prior.values() if p["zero_levels"] is not None))
    doc = {"schema": "work7-rank4-gap-census-v1-closure", "claim_universe":
           "all positive definite integral rank-four lattices with 25 <= det <= 279",
           "even_lattice_det_mod_4_counts_25_279": dict(sorted(mod4.items())),
           "decomposable_classes": "not decided; uniform orthogonal-sum corollary",
           "prior_complete_lists": prior, "new_regions": regions, "new_totals": totals,
           "symmetry_converse_holds_for_25_to_279": closed,
           "consequence": ("the minimum determinant of a rank-four counterexample to the symmetry converse is 280"
                           if closed else "gap not closed; see zeros_without_witness / unresolved"),
           "preserved_closure_agrees": True}
    write_json(Path(out) / "closure-25-279.json", doc)
    print(json.dumps({k: v for k, v in doc.items() if k not in ("prior_complete_lists", "new_regions")}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("task", choices=("run", "replay", "summary"))
    parser.add_argument("--regions", nargs="*", default=list(REGIONS))
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--chunks", type=Path, default=CHUNKS, help="preserved chunks, summaries and replays")
    parser.add_argument("--out", type=Path, default=Path("scratch/rank4-gap-census"), help="fresh outputs")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.task == "run":
        run(args.regions, args.workers, args.out)
    elif args.task == "replay":
        replay(args.regions, args.workers, args.chunks, args.out)
    else:
        summary(args.chunks, args.out)


if __name__ == "__main__":
    main()
