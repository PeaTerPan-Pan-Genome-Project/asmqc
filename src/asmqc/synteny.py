"""BUSCO-anchored synteny against the reference (Caméor v2), report only.

The anchors (refs/cameor_v2_busco_anchors.tsv) are the reference's Complete
single-copy BUSCOs on chr1-chr7. Joining them with the tested assembly's
Complete BUSCOs gives one point per shared gene: enough to see swapped,
reversed or translocated chromosomes without aligning the genomes. Nothing
here enters qc_summary.tsv (SPEC §8.2, §8.6).
"""

from collections import Counter
from pathlib import Path

from asmqc.validate import CHROMOSOMES

POINT_COLUMNS = ["busco_id", "ref_chromosome", "ref_pos", "seq_id", "pos", "same_strand"]
TABLE_COLUMNS = ["chromosome", "length", "n_busco", "best_ref_chromosome", "frac_on_best",
                 "orientation", "rho"]


def anchors_path() -> Path:
    from asmqc.runner import asmqc_home

    return asmqc_home() / "refs" / "cameor_v2_busco_anchors.tsv"


def read_anchors(path: Path) -> tuple[dict[str, int], dict[str, tuple[str, float, str]]]:
    """(reference chromosome lengths, busco_id -> (chromosome, midpoint, strand))."""
    lengths, anchors = {}, {}
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if line.startswith("# length"):
            lengths[f[1]] = int(f[2])
        elif not line.startswith("#") and f[0] != "busco_id":
            anchors[f[0]] = (f[1], (int(f[2]) + int(f[3])) / 2, f[4])
    return lengths, anchors


def points(rows: list[dict], anchors: dict) -> list[tuple]:
    """Shared Complete BUSCOs: (busco_id, ref_chr, ref_pos, seq, pos, same_strand)."""
    out = []
    for r in rows:
        if r["status"] != "Complete" or r["busco_id"] not in anchors:
            continue
        s, e = sorted((int(r["Gene Start"]), int(r["Gene End"])))
        ref_chr, ref_pos, ref_strand = anchors[r["busco_id"]]
        out.append((r["busco_id"], ref_chr, round(ref_pos), r["sequence"], (s + e) // 2,
                    "yes" if r["Strand"] == ref_strand else "no"))
    return out


def _ranks(xs: list[float]) -> list[float]:
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    for rank, i in enumerate(order):
        ranks[i] = rank
    return ranks


def spearman(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 3:
        return None
    rx, ry = _ranks(xs), _ranks(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    cov = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    vx = sum((a - mx) ** 2 for a in rx)
    vy = sum((b - my) ** 2 for b in ry)
    return cov / (vx * vy) ** 0.5 if vx and vy else None


def table(pts: list[tuple], lengths: dict[str, int]) -> list[tuple]:
    """Per chromosome: best-matching reference chromosome and orientation.

    Orientation is the sign of the rank correlation of positions on the best
    reference chromosome: forward (rho >= 0.5), reverse (<= -0.5), else mixed.
    Unplaced sequences are pooled in one row.
    """
    out = []
    for c in [*CHROMOSOMES, "unplaced"]:
        mine = [p for p in pts if (p[3] == c if c != "unplaced" else p[3] not in CHROMOSOMES)]
        if not mine:
            out.append((c, lengths.get(c, "NA"), 0, "NA", "NA", "NA", "NA"))
            continue
        best, n_best = Counter(p[1] for p in mine).most_common(1)[0]
        on_best = [p for p in mine if p[1] == best]
        rho = spearman([p[2] for p in on_best], [p[4] for p in on_best]) if c != "unplaced" \
            else None
        orient = ("NA" if rho is None else "forward" if rho >= 0.5 else
                  "reverse" if rho <= -0.5 else "mixed")
        length = lengths.get(c, sum(lengths.get(s, 0) for s in {p[3] for p in mine}))
        out.append((c, length, len(mine), best, f"{n_best / len(mine):.3f}", orient,
                    "NA" if rho is None else f"{rho:.3f}"))
    return out
