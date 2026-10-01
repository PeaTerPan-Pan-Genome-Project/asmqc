"""Shared `scan` stage: one pass over the prepared FASTA (SPEC §8.1, §8.2).

Writes to OUTDIR:
  sequences.tsv  per-sequence composition, N ends, line widths
  gaps.tsv       N-runs >= 10 bp (seq_id, start, end, length; 1-based, inclusive)
  file.json      file-level facts: CRLF, line widths, bad names, duplicates

M1 interprets these; M2 uses them for the N-split contigs. Holds one sequence
in memory at a time.

Usage: python -m asmqc.scan FASTA OUTDIR
"""

import csv
import json
import re
import sys
from collections import Counter
from dataclasses import astuple, dataclass, fields
from pathlib import Path

from asmqc.params import PARAMS

IUPAC = b"RYSWKMBDHVryswkmbdhv"
KNOWN = b"ACGTNacgtn" + IUPAC
LOWER = bytes(range(ord("a"), ord("z") + 1))
GAP_RE = re.compile(rb"[Nn]{%d,}" % PARAMS["m01"]["gap_min_bp"])


@dataclass
class SeqStats:
    seq_id: str
    length: int
    a: int
    c: int
    g: int
    t: int
    n: int
    iupac: int
    other: int
    lower: int
    leading_n: int
    trailing_n: int
    line_width: int  # 0 for a single-line sequence
    width_consistent: bool


SEQ_COLUMNS = [f.name for f in fields(SeqStats)]


def _count(seq: bytes, base: bytes) -> int:
    return seq.count(base) + seq.count(base.lower())


def seq_stats(name: str, lines: list[bytes]) -> SeqStats:
    seq = b"".join(lines)
    other = len(seq.translate(None, KNOWN))
    iupac = sum(seq.count(bytes([b])) for b in IUPAC)
    width = len(lines[0]) if len(lines) > 1 else 0
    consistent = all(len(x) == width for x in lines[:-1]) and (
        len(lines) < 2 or 0 < len(lines[-1]) <= width)
    return SeqStats(
        seq_id=name, length=len(seq),
        a=_count(seq, b"A"), c=_count(seq, b"C"), g=_count(seq, b"G"), t=_count(seq, b"T"),
        n=_count(seq, b"N"), iupac=iupac, other=other,
        lower=len(seq) - len(seq.translate(None, LOWER)),
        leading_n=len(seq) - len(seq.lstrip(b"Nn")),
        trailing_n=len(seq) - len(seq.rstrip(b"Nn")) if seq.strip(b"Nn") else 0,
        line_width=width, width_consistent=consistent)


def bad_name(name: str) -> bool:
    return not name or any(not 33 <= ord(ch) <= 126 for ch in name)


def scan(fasta: Path, outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    crlf = False
    names: list[str] = []
    with (fasta.open("rb") as fh, (outdir / "sequences.tsv").open("w", newline="") as fs,
          (outdir / "gaps.tsv").open("w", newline="") as fg):
        ws = csv.writer(fs, delimiter="\t", lineterminator="\n")
        wg = csv.writer(fg, delimiter="\t", lineterminator="\n")
        ws.writerow(SEQ_COLUMNS)
        wg.writerow(["seq_id", "start", "end", "length"])

        def flush(name, lines):
            ws.writerow(astuple(seq_stats(name, lines)))
            for m in GAP_RE.finditer(b"".join(lines)):
                wg.writerow([name, m.start() + 1, m.end(), m.end() - m.start()])

        name, lines = None, []
        for raw in fh:
            if raw.endswith(b"\r\n"):
                crlf = True
            line = raw.rstrip(b"\r\n")
            if line.startswith(b">"):
                if name is not None:
                    flush(name, lines)
                fields_ = line[1:].split(None, 1)
                name = fields_[0].decode(errors="replace") if fields_ else ""
                names.append(name)
                lines = []
            elif name is not None:
                lines.append(line)
        if name is not None:
            flush(name, lines)

    dup = {n: c for n, c in Counter(names).items() if c > 1}
    facts = {
        "crlf": crlf,
        "bad_names": sorted({n for n in names if bad_name(n)}),
        "duplicate_names": dup,
    }
    (outdir / "file.json").write_text(json.dumps(facts, indent=1) + "\n")


def read_sequences(scan_dir: Path) -> list[SeqStats]:
    out = []
    with (scan_dir / "sequences.tsv").open(newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            vals = {k: (r[k] == "True" if k == "width_consistent" else
                        r[k] if k == "seq_id" else int(r[k])) for k in SEQ_COLUMNS}
            out.append(SeqStats(**vals))
    return out


def read_gaps(scan_dir: Path) -> list[tuple[str, int, int]]:
    with (scan_dir / "gaps.tsv").open(newline="") as fh:
        return [(r["seq_id"], int(r["start"]), int(r["end"]))
                for r in csv.DictReader(fh, delimiter="\t")]


if __name__ == "__main__":
    scan(Path(sys.argv[1]), Path(sys.argv[2]))
