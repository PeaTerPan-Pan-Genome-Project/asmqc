"""Module 4: telomeres (SPEC §8.4).

tidk convention (checked on the synthetic data, §13.3): forward_repeat_number
counts the search string TTTAGGG, reverse_repeat_number its reverse
complement CCCTAAA; `window` is the window end, capped at the sequence length.
A capped start arm is CCCTAAA-dominated, a capped end arm TTTAGGG-dominated.

Usage: python -m asmqc.m04_telomeres --windows TSV --fai FAI --outdir DIR --work DIR
"""

import argparse
import csv
import shutil
from dataclasses import dataclass
from pathlib import Path

from asmqc import flags as fl
from asmqc import module
from asmqc.params import PARAMS
from asmqc.validate import CHROMOSOMES

P = PARAMS["m04"]
M = "m04"
MOTIF_LEN = len(P["motif"])


@dataclass
class Band:
    seq: str
    start: int  # 1-based, first base of the first window
    end: int  # last base of the last window
    fwd: int
    rev: int

    @property
    def approx_bp(self) -> int:
        return (self.fwd + self.rev) * MOTIF_LEN


def read_windows(path: Path) -> list[tuple[str, int, int, int]]:
    """(seq, window_end, forward, reverse), sorted by sequence order then position."""
    with path.open(newline="") as fh:
        rows = [(r["id"], int(r["window"]), int(r["forward_repeat_number"]),
                 int(r["reverse_repeat_number"])) for r in csv.DictReader(fh, delimiter="\t")]
    order = {s: i for i, s in enumerate(dict.fromkeys(r[0] for r in rows))}
    return sorted(set(rows), key=lambda r: (order[r[0]], r[1]))


def bands(windows: list[tuple[str, int, int, int]]) -> dict[str, list[Band]]:
    """Merge adjacent telomeric windows (forward + reverse >= min_repeats)."""
    w = P["window"]
    out: dict[str, list[Band]] = {}
    for seq, end, f, r in windows:
        if f + r < P["min_repeats"]:
            continue
        start = (end - 1) // w * w + 1
        bs = out.setdefault(seq, [])
        if bs and bs[-1].end == start - 1:
            b = bs[-1]
            b.end, b.fwd, b.rev = end, b.fwd + f, b.rev + r
        else:
            bs.append(Band(seq, start, end, f, r))
    return out


def classify_arm(band: Band | None, distance: int | None, capped_motif: str) -> str:
    if band is None or distance > P["terminal_bp"]:
        return "absent"
    rev_dominant = band.rev > band.fwd
    capped = rev_dominant if capped_motif == "rev" else band.fwd > band.rev
    return "capped" if capped else "wrong_orientation"


def evaluate(windows, lengths: dict[str, int]):
    by_seq = bands(windows)
    term = P["terminal_bp"]
    arms, flags, interstitial, unplaced = [], [], [], []
    for c in CHROMOSOMES:
        n, bs = lengths[c], by_seq.get(c, [])
        first, last = (bs[0], bs[-1]) if bs else (None, None)
        for arm, band, dist, motif in (
                ("start", first, first.start - 1 if first else None, "rev"),
                ("end", last, n - last.end if last else None, "fwd")):
            status = classify_arm(band, dist, motif)
            if status == "absent":
                arms.append((c, arm, status, "NA", "NA", "NA", "NA"))
                continue
            arms.append((c, arm, status, dist, band.approx_bp, band.fwd, band.rev))
            if status == "wrong_orientation":
                flags.append(fl.make(M, "wrong_orientation_telomere", c, band.start, band.end,
                                     arm, f"{arm} arm band is dominated by the "
                                          f"{'TTTAGGG' if motif == 'rev' else 'CCCTAAA'} strand"))
        for b in bs:
            if b.start - 1 > term and n - b.end > term:
                interstitial.append((c, b.start, b.end, b.fwd, b.rev, b.approx_bp))
                flags.append(fl.make(M, "interstitial_telomere", c, b.start, b.end,
                                     b.approx_bp, "telomeric band > 50 kb from both ends"))
    for seq, bs in by_seq.items():
        if seq in CHROMOSOMES:
            continue
        n = lengths[seq]
        for end_name, b, dist in (("start", bs[0], bs[0].start - 1),
                                  ("end", bs[-1], n - bs[-1].end)):
            if dist <= term:
                unplaced.append((seq, n, end_name, dist, b.fwd, b.rev, b.approx_bp))
    status_of = {(a[0], a[1]): a[2] for a in arms}
    values = {
        "m04_capped_arms": sum(a[2] == "capped" for a in arms),
        "m04_t2t_chromosomes": sum(status_of[(c, "start")] == status_of[(c, "end")] == "capped"
                                   for c in CHROMOSOMES),
        "m04_wrong_orientation_arms": sum(a[2] == "wrong_orientation" for a in arms),
        "m04_interstitial_arrays": len(interstitial),
        "m04_unplaced_with_telomere": len({u[0] for u in unplaced}),
    }
    return values, flags, arms, interstitial, unplaced, by_seq


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--windows", type=Path, required=True)
    ap.add_argument("--fai", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    args = ap.parse_args(argv)

    lengths = module.read_fai(args.fai)
    values, flags, arms, inter, unpl, _ = evaluate(read_windows(args.windows), lengths)
    out = args.outdir
    module.write_tsv(out / "telomeres.tsv",
                     ["chromosome", "arm", "status", "distance_from_end_bp", "approx_array_bp",
                      "fwd_repeats", "rev_repeats"], arms)
    module.write_tsv(out / "interstitial.tsv",
                     ["chromosome", "start", "end", "fwd_repeats", "rev_repeats",
                      "approx_array_bp"], inter)
    module.write_tsv(out / "telomeres_unplaced.tsv",
                     ["seq_id", "length", "end", "distance_from_end_bp", "fwd_repeats",
                      "rev_repeats", "approx_array_bp"], unpl)
    shutil.copyfile(args.windows, out / "tidk_windows.tsv")
    module.finish(args.work, M, values, flags)


if __name__ == "__main__":
    main()
