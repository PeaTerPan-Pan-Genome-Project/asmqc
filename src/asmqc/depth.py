"""Read depth helpers for M9 and M11.

median_depth: median of `samtools depth -a` sampled every 1 kb over chr1-chr7
(SPEC §8.9, §8.11).
callable_bed: positions with depth in [low, high] as merged BED intervals,
one samtools depth | gawk pipeline per sequence, run in parallel.
"""

import statistics
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from asmqc.validate import CHROMOSOMES

# Merge consecutive positions with depth in [lo, hi] into BED intervals.
GAWK = r"""
$3 >= lo && $3 <= hi {
    if ($1 == c && $2 == e + 1) { e = $2; next }
    if (c != "") print c "\t" s - 1 "\t" e
    c = $1; s = $2; e = $2; next
}
END { if (c != "") print c "\t" s - 1 "\t" e }
"""


def median_depth(bam: Path, lengths: dict[str, int], work: Path, step: int = 1000) -> float:
    work.mkdir(parents=True, exist_ok=True)
    sample = work / "depth_sample.bed"
    with sample.open("w") as fh:
        for c in CHROMOSOMES:
            for pos in range(1, lengths[c] + 1, step):
                fh.write(f"{c}\t{pos - 1}\t{pos}\n")
    out = subprocess.run(["samtools", "depth", "-a", "-b", str(sample), str(bam)],
                         capture_output=True, text=True, check=True).stdout
    depths = [int(line.split("\t")[2]) for line in out.splitlines()]
    return float(statistics.median(depths)) if depths else 0.0


def callable_bed(bam: Path, lengths: dict[str, int], low: float, high: float, out: Path,
                 threads: int) -> int:
    """Write the callable BED; return callable bp."""
    parts = out.parent / "callable_parts"
    parts.mkdir(parents=True, exist_ok=True)

    def one(i_seq: tuple[int, str]) -> Path:
        i, seq = i_seq
        part = parts / f"{i}.bed"
        cmd = (f"samtools depth -a -r {seq!r}:1-{lengths[seq]} {bam}"
               f" | gawk -v lo={low} -v hi={high} '{GAWK}' > {part}")
        subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)
        return part

    with ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        files = list(pool.map(one, enumerate(lengths)))
    total = 0
    with out.open("w") as fh:
        for f in files:
            for line in f.read_text().splitlines():
                _, s, e = line.split("\t")
                total += int(e) - int(s)
                fh.write(line + "\n")
            f.unlink()
    parts.rmdir()
    return total
