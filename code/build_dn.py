"""Build and verify compact complete D2..D10 certificate streams."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from time import perf_counter

from dn_exact import (coefficient_tables, expand_record, header, records,
                      verify_records)
from exact_theta import verify

ROOT = Path(__file__).resolve().parent.parent


def build(rank):
    start = perf_counter()
    document = header(rank)
    coefficients, norms, histograms, stats = coefficient_tables(document)
    path = ROOT / "data" / "dn" / f"d{rank}.jsonl.gz"
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    samples = {}
    with temporary.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as stream:
        stream.write((json.dumps(document, separators=(",", ":")) + "\n").encode("ascii"))

        def saved_rows():
            for row in records(document, coefficients, norms):
                stream.write((json.dumps(row, separators=(",", ":")) + "\n").encode("ascii"))
                kind = row[2]
                if kind not in samples or (kind == 0 and row[3] > samples[kind][3]):
                    samples[kind] = row
                yield row

        result = verify_records(document, saved_rows(), coefficients, norms, histograms, stats, "ldl")
    temporary.replace(path)
    result["certificate_stream_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    result["certificate_stream"] = f"data/dn/d{rank}.jsonl.gz"
    result["seconds"] = round(perf_counter() - start, 4)
    result["compressed_bytes"] = path.stat().st_size
    result_path = ROOT / "reports" / f"dn-d{rank}-ldl.json"
    result_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if rank in (2, 3, 4, 10):
        for kind, row in samples.items():
            doc = expand_record(document, row)
            verify(doc)
            label = {0: "nonzero", 1: "product-zero", 2: "balanced-zero"}[kind]
            (path.parent / f"d{rank}-{label}.json").write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"rank": rank, "even_classes": result["even_characteristics"],
                      "zeros": result["even_vanishing"], "nonzeros": result["even_nonvanishing"],
                      "seconds": result["seconds"], "bytes": result["compressed_bytes"]}), flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ranks", type=int, nargs="+", choices=range(2, 11))
    args = parser.parse_args()
    results = [build(rank) for rank in args.ranks]
    (ROOT / "reports" / "dn-audit-summary.json").write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
