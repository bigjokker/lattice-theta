"""Stream-verifier for the complete Dn finite audit; choose LDL or ambient shells."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from time import perf_counter

from dn_exact import coefficient_tables, validate_header, verify_records
from exact_theta import EnumerationLimit, InvalidCertificate


def verify_stream(path, backend="ambient", max_nodes=1_000_000):
    start = perf_counter()
    with gzip.open(path, "rt", encoding="ascii") as stream:
        document = json.loads(stream.readline())
        validate_header(document)
        coefficients, norms, histograms, stats = coefficient_tables(document, backend, max_nodes)
        result = verify_records(document, (json.loads(line) for line in stream), coefficients, norms,
                                histograms, stats, backend)
    result["certificate_stream_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    result["seconds"] = round(perf_counter() - start, 4)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stream", type=Path)
    parser.add_argument("--backend", choices=("ldl", "ambient"), default="ambient")
    parser.add_argument("--max-nodes", type=int, default=1_000_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        result = verify_stream(args.stream, args.backend, args.max_nodes)
    except EnumerationLimit as error:
        result = {"verdict": "unresolved", "reason": str(error)}
    except (InvalidCertificate, OSError, ValueError, EOFError) as error:
        result = {"verdict": "invalid_certificate", "reason": str(error)}
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result))
    return {"classification_verified": 0, "invalid_certificate": 1, "unresolved": 2}[result["verdict"]]


if __name__ == "__main__":
    raise SystemExit(main())
