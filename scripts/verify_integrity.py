"""Verify all distributed files against checksums.sha256."""

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def main():
    count = 0
    for line in (ROOT / "checksums.sha256").read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        target = (ROOT / name).resolve()
        if not target.is_relative_to(ROOT):
            raise ValueError("Checksum path escapes repository")
        if hashlib.sha256(target.read_bytes()).hexdigest() != digest:
            raise ValueError(f"File integrity check failed: {name}")
        count += 1
    print(f"Verified {count} distributed files.")


if __name__ == "__main__":
    main()
