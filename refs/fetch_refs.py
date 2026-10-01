#!/usr/bin/env python3
"""Download the reference data listed in refs.tsv and verify md5s (SPEC §6).

Usage: fetch_refs.py OUTDIR [NAME ...]

Exits non-zero on any md5 mismatch; a mismatching download is removed.
Already present files with the correct md5 are not downloaded again.
Standard library only: runs in the image build before any env exists.
"""

import csv
import hashlib
import sys
import tarfile
import urllib.request
from pathlib import Path

TABLE = Path(__file__).with_name("refs.tsv")


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def fetch(row: dict, outdir: Path) -> None:
    name, url, expected = row["name"], row["url"], row["md5"]
    if row["unpack"] == "tar.gz":
        target = outdir / "downloads" / f"{name}.tar.gz"
    else:
        target = outdir / row["dest"]
    target.parent.mkdir(parents=True, exist_ok=True)
    if not (target.exists() and md5(target) == expected):
        print(f"fetching {name}", file=sys.stderr)
        tmp = target.with_suffix(target.suffix + ".part")
        with urllib.request.urlopen(url, timeout=300) as resp, tmp.open("wb") as out:
            while block := resp.read(1 << 20):
                out.write(block)
        tmp.rename(target)
    got = md5(target)
    if got != expected:
        target.unlink()
        sys.exit(f"md5 mismatch for {name}: expected {expected}, got {got}")
    if row["unpack"] == "tar.gz":
        dest = outdir / row["dest"]
        dest.mkdir(parents=True, exist_ok=True)
        with tarfile.open(target) as tar:
            tar.extractall(dest, filter="data")
    print(f"{name}\t{got}\tOK", file=sys.stderr)


def main(argv: list[str]) -> None:
    if not argv:
        sys.exit(__doc__)
    outdir, names = Path(argv[0]), set(argv[1:])
    with TABLE.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    unknown = names - {r["name"] for r in rows}
    if unknown:
        sys.exit(f"unknown reference(s): {', '.join(sorted(unknown))}")
    for row in rows:
        if not names or row["name"] in names:
            fetch(row, outdir)


if __name__ == "__main__":
    main(sys.argv[1:])
