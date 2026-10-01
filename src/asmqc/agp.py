"""AGP 2.0 parsing, consistency check and cut coordinates (SPEC §8.1, §8.2)."""

import re
from dataclasses import dataclass
from pathlib import Path

GAP_TYPES = ("N", "U")
SUBSEQ_RE = re.compile(r"_subseq_(\d+):(\d+)$")


@dataclass(frozen=True)
class Row:
    obj: str
    obj_beg: int
    obj_end: int
    part: int
    type: str
    comp: str  # component id, or gap length for N/U lines
    comp_beg: int  # 0 for gaps
    comp_end: int  # 0 for gaps

    @property
    def is_gap(self) -> bool:
        return self.type in GAP_TYPES

    @property
    def span(self) -> int:
        return self.obj_end - self.obj_beg + 1


def parse(path: Path) -> tuple[list[Row], list[str]]:
    """Rows in file order, and problems with malformed lines."""
    rows, problems = [], []
    for i, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\r").split("\t")
        try:
            if f[4] in GAP_TYPES:
                rows.append(Row(f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5], 0, 0))
            else:
                rows.append(Row(f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5],
                                int(f[6]), int(f[7])))
        except (IndexError, ValueError):
            problems.append(f"line {i}: malformed AGP line")
    return rows, problems


def check(rows: list[Row], lengths: dict[str, int]) -> list[tuple[str, str]]:
    """(object, problem) for every inconsistency between the AGP and the FASTA."""
    problems: list[tuple[str, str]] = []
    by_obj: dict[str, list[Row]] = {}
    for r in rows:
        by_obj.setdefault(r.obj, []).append(r)
    for obj in by_obj:
        if obj not in lengths:
            problems.append((obj, "AGP object not in the FASTA"))
    for name in lengths:
        if name not in by_obj:
            problems.append((name, "FASTA sequence not in the AGP"))
    for obj, rs in by_obj.items():
        if obj not in lengths:
            continue
        pos = 0
        for r in rs:
            if r.obj_beg != pos + 1:
                problems.append((obj, f"part {r.part} starts at {r.obj_beg}, expected {pos + 1}"))
                break
            if r.obj_end < r.obj_beg:
                problems.append((obj, f"part {r.part} ends before it starts"))
                break
            if r.is_gap:
                if not r.comp.isdigit() or int(r.comp) != r.span:
                    problems.append((obj, f"gap part {r.part} length {r.comp} != span {r.span}"))
            elif r.comp_beg < 1 or r.comp_end - r.comp_beg + 1 != r.span:
                n = r.comp_end - r.comp_beg + 1
                problems.append((obj, f"component part {r.part} length {n} != span {r.span}"))
            pos = r.obj_end
        else:
            if pos != lengths[obj]:
                problems.append((obj, f"AGP length {pos} != FASTA length {lengths[obj]}"))
    return problems


def cut_coordinates(rows: list[Row]) -> list[tuple[Row, int]]:
    """Cut coordinates for the round-number diagnostic (SPEC §8.1)."""
    cuts = []
    for r in rows:
        if r.is_gap:
            continue
        if r.comp_beg > 1:
            cuts.append((r, r.comp_beg - 1))
        m = SUBSEQ_RE.search(r.comp)
        if m:
            a, b = int(m[1]), int(m[2])
            if a > 1:
                cuts.append((r, a - 1))
            cuts.append((r, b))
        else:
            cuts.append((r, r.comp_end))
    return cuts
