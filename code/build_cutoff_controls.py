"""Write separate cutoff control evidence and metadata audit; no census rewrite."""

import hashlib
import json
from pathlib import Path

from exact_theta import identity, require
from rank3_cutoff import build, independent_modular_check, metadata, replay


def save(path, document):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


def main():
    G20 = [[2, -1, -1], [-1, 4, -1], [-1, -1, 4]]
    controls = [("orthogonal modular zero", identity(3), [1,1,0], [1,1,0], {}),
                ("inclusive endpoint missing", identity(3), [1,1,0], [1,1,0], {"norm4_bound": 5}),
                ("enumeration interrupted", identity(3), [1,1,0], [1,1,0], {"max_nodes": 1}),
                ("I3 nonzero", identity(3), [0,1,0], [1,0,0], {}),
                ("A3 nonzero", [[2,-1,0],[-1,2,-1],[0,-1,2]], [0,1,0], [1,0,0], {}),
                ("2I3 nonzero", [[2,0,0],[0,2,0],[0,0,2]], [0,1,0], [1,0,0], {}),
                ("det20 cancelled first shell", G20, [0,1,1], [1,0,0], {"norm4_bound": 6}),
                ("det20 complete cutoff", G20, [0,1,1], [1,0,0], {}),
                ("constant term nonzero", identity(3), [0,0,0], [1,0,0], {"norm4_bound": 0})]
    pack = {"schema": "work7-rank3-cutoff-controls-v1", "certificates": [
        {"name": name, **build(G, a, b, **kwargs)} for name, G, a, b, kwargs in controls]}
    pack_path = Path("data/census-rank3/cutoff-controls.json")
    save(pack_path, pack)
    digest = hashlib.sha256(pack_path.read_bytes()).hexdigest()
    for backend in ("box", "ldl"):
        results = [{"name": cert["name"], **replay(cert, backend=backend)} for cert in pack["certificates"]]
        report = {"schema": "work7-rank3-cutoff-controls-result-v1", "backend": backend,
                  "pack_sha256": digest, "all_claims_replayed": all(r["certificate_replay_complete"] for r in results),
                  "reports": results}
        save(Path(f"reports/census-rank3/cutoff-controls-{backend}.json"), report)
    census_path = Path("data/census-rank3/det24-bound32.json")
    raw = census_path.read_bytes()
    census = json.loads(raw)
    rows = []
    for i, row in enumerate(census["lattices"]):
        for stored in row["modular_inputs"]:
            item = metadata(row["G"], stored["b"])
            require(all(stored[k] == item[k] for k in stored if k != "cutoff"), "historical modular metadata mismatch")
            B = independent_modular_check(row["G"], stored["b"], item)
            rows.append({"lattice": i, "b": stored["b"], "N0": item["N0"],
                         "unsigned_level": item["unsigned_level"], "index": item["index"], "cutoff": B})
    audit = {"schema": "work7-rank3-cutoff-index-audit-v1", "census_pack_sha256": hashlib.sha256(raw).hexdigest(),
             "records_checked": len(rows), "minimum_cutoff": min(r["cutoff"] for r in rows),
             "maximum_cutoff": max(r["cutoff"] for r in rows),
             "coefficient_enumeration_performed": False, "rows": rows}
    save(Path("reports/census-rank3/cutoff-index-audit.json"), audit)
    print(json.dumps({"controls": len(controls), "metadata_records_checked": len(rows),
                      "minimum_cutoff": audit["minimum_cutoff"], "maximum_cutoff": audit["maximum_cutoff"],
                      "pack_sha256": digest}))


if __name__ == "__main__":
    main()
