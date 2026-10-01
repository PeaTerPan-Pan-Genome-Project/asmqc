"""Shared `prep` stage (SPEC §5.2): working copy of the assembly and AGP.

Decompresses the assembly once, applies the --chromosomes renaming to the
headers (sequence bytes unchanged) and records the md5 of the decompressed
original bytes (`assembly_md5`). The AGP object names are renamed the same way.

Usage: python -m asmqc.prep --assembly A --map MAP.json --outdir DIR [--agp AGP]
"""

import argparse
import hashlib
import json
from pathlib import Path

from asmqc.seqio import open_bytes


def rename_header(line: bytes, mapping: dict[bytes, bytes]) -> bytes:
    name, sep, rest = line[1:].partition(b" ")
    name = name.rstrip(b"\r\n")
    new = mapping.get(name)
    if new is None:
        return line
    eol = b"\r\n" if line.endswith(b"\r\n") else b"\n" if line.endswith(b"\n") else b""
    rest = rest.rstrip(b"\r\n")
    return b">" + new + (b" " + rest if sep else b"") + eol


def prep_fasta(src: Path, dest: Path, mapping: dict[str, str]) -> str:
    """Write the renamed copy; return the md5 of the decompressed source."""
    bmap = {k.encode(): v.encode() for k, v in mapping.items()}
    h = hashlib.md5()
    with open_bytes(src) as fin, dest.open("wb") as fout:
        for line in fin:
            h.update(line)
            fout.write(rename_header(line, bmap) if line.startswith(b">") and bmap else line)
    return h.hexdigest()


def prep_agp(src: Path, dest: Path, mapping: dict[str, str]) -> None:
    with open_bytes(src) as fin, dest.open("wb") as fout:
        for line in fin:
            if line.startswith(b"#") or not line.strip():
                fout.write(line)
                continue
            obj, sep, rest = line.partition(b"\t")
            fout.write(mapping.get(obj.decode(), obj.decode()).encode() + sep + rest)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--assembly", type=Path, required=True)
    ap.add_argument("--agp", type=Path)
    ap.add_argument("--map", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    args = ap.parse_args(argv)

    mapping = json.loads(args.map.read_text())
    args.outdir.mkdir(parents=True, exist_ok=True)
    md5 = prep_fasta(args.assembly, args.outdir / "asm.fa", mapping)
    if args.agp:
        prep_agp(args.agp, args.outdir / "asm.agp", mapping)
    (args.outdir / "asm.md5").write_text(md5 + "\n")


if __name__ == "__main__":
    main()
