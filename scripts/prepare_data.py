"""Expand lossless compressed inputs and reports, or remove verified generated copies."""

import argparse
import gzip
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MANIFESTS = {
    "data": "data/compressed-files.json",
    "reports": "reports/compressed-reports.json",
}


def prepare(clean=False, kinds=("data",)):
    changed = 0
    for kind in kinds:
        base = ROOT / kind
        manifest = json.loads((ROOT / MANIFESTS[kind]).read_text(encoding="utf-8"))
        for entry in manifest["files"]:
            target = (ROOT / entry["path"]).resolve()
            compressed = (ROOT / entry["gzip_path"]).resolve()
            if not target.is_relative_to(base) or not compressed.is_relative_to(base):
                raise ValueError(f"Dataset path escapes {kind} directory")
            if target.exists():
                if hashlib.sha256(target.read_bytes()).hexdigest() != entry["sha256"]:
                    raise ValueError(f"Existing dataset differs from certificate: {entry['path']}")
                if clean:
                    target.unlink()
                    changed += 1
                continue
            if clean:
                continue
            raw = gzip.decompress(compressed.read_bytes())
            if len(raw) != entry["size"] or hashlib.sha256(raw).hexdigest() != entry["sha256"]:
                raise ValueError(f"Compressed certificate hash mismatch: {entry['gzip_path']}")
            target.write_bytes(raw)
            changed += 1
    return changed


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clean", action="store_true", help="Remove only unchanged expanded copies")
    parser.add_argument("--reports", action="store_true",
                        help="Also expand the large preserved reports listed in reports/compressed-reports.json")
    args = parser.parse_args()
    kinds = ("data", "reports") if args.reports else ("data",)
    count = prepare(args.clean, kinds)
    print(f"{'Removed' if args.clean else 'Prepared'} {count} expanded file(s).")
