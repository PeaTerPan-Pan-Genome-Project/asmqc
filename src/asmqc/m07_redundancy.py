"""Module 7: redundancy of unplaced scaffolds (SPEC §8.7). A measurement, never a purge.

`duplicate` and `partial_overlap` need one primary PAF record (one collinear
chain to one locus): summed coverage would call purely repetitive scaffolds
duplicates in an ~85 % repetitive genome.

Usage:
  python -m asmqc.m07_redundancy lists --fai FAI --work DIR
      writes chroms.txt, unplaced_ge1kb.txt (region lists for samtools faidx -r)
  python -m asmqc.m07_redundancy classify --fai FAI --paf PAF --total-bp N
      --outdir DIR --work DIR
"""

import argparse
from pathlib import Path

from asmqc import module
from asmqc.m05_organelle_rdna import covered
from asmqc.params import PARAMS
from asmqc.validate import CHROMOSOMES

P = PARAMS["m07"]
M = "m07"
CLASSES = ("duplicate", "partial_overlap", "repeat_like", "unique", "short")


def classify(paf_lines: list[str], lengths: dict[str, int]) -> dict[str, tuple]:
    """unplaced seq -> (class, chromosome, tstart, tend, identity, mapq, coverage)."""
    best: dict[str, tuple] = {}
    hq: dict[str, list] = {}
    for line in paf_lines:
        f = line.split("\t")
        q, qlen, qs, qe = f[0], int(f[1]), int(f[2]), int(f[3])
        ident = int(f[9]) / int(f[10]) if int(f[10]) else 0.0
        mapq = int(f[11])
        if ident < P["min_identity"]:
            continue
        hq.setdefault(q, []).append((qs, qe))
        primary = "tp:A:P" in f[12:]
        cov = (qe - qs) / qlen
        if primary and mapq >= P["min_mapq"]:
            rec = (cov, f[5], int(f[7]) + 1, int(f[8]), ident, mapq)
            if q not in best or rec[0] > best[q][0]:
                best[q] = rec

    out = {}
    for seq, n in lengths.items():
        if seq in CHROMOSOMES:
            continue
        if n < P["min_len_bp"]:
            out[seq] = ("short", "", "", "", "", "", "")
            continue
        b = best.get(seq)
        if b and b[0] >= P["duplicate_min_cov"]:
            cls = "duplicate"
        elif b and b[0] >= P["partial_min_cov"]:
            cls = "partial_overlap"
        elif covered(hq.get(seq, [])) >= P["repeat_min_cov"] * n:
            cls = "repeat_like"
        else:
            cls = "unique"
        if b and cls in ("duplicate", "partial_overlap"):
            cov, chrom, ts, te, ident, mapq = b
            out[seq] = (cls, chrom, ts, te, f"{ident:.4f}", mapq, f"{cov:.4f}")
        else:
            out[seq] = (cls, "", "", "", "", "", "")
    return out


def values_from(classes: dict[str, tuple], lengths: dict[str, int], total_bp: int) -> dict:
    v = {}
    for cls in CLASSES:
        seqs = [s for s, c in classes.items() if c[0] == cls]
        v[f"m07_n_{cls}"] = len(seqs)
        v[f"m07_{cls}_bp"] = sum(lengths[s] for s in seqs)
    v["m07_total_minus_duplicate_bp"] = total_bp - v["m07_duplicate_bp"]
    return v


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=("lists", "classify"))
    ap.add_argument("--fai", type=Path, required=True)
    ap.add_argument("--paf", type=Path)
    ap.add_argument("--outdir", type=Path)
    ap.add_argument("--work", type=Path, required=True)
    args = ap.parse_args(argv)
    lengths = module.read_fai(args.fai)
    args.work.mkdir(parents=True, exist_ok=True)

    if args.step == "lists":
        unplaced = [s for s, n in lengths.items()
                    if s not in CHROMOSOMES and n >= P["min_len_bp"]]
        (args.work / "chroms.txt").write_text("".join(f"{c}\n" for c in CHROMOSOMES))
        (args.work / "unplaced_ge1kb.txt").write_text("".join(f"{s}\n" for s in unplaced))
        return

    total = sum(lengths.values())
    if not (args.work / "unplaced_ge1kb.txt").read_text().strip():
        classes = classify([], lengths)  # only short ones, if any
        module.finish(args.work, M, {}, status="skipped_no_input",
                      reason="no unplaced sequence >= 1 kb")
        _write_table(args.outdir, classes, lengths)
        return
    classes = classify(args.paf.read_text().splitlines(), lengths)
    _write_table(args.outdir, classes, lengths)
    module.finish(args.work, M, values_from(classes, lengths, total))


def _write_table(outdir: Path, classes: dict[str, tuple], lengths: dict[str, int]) -> None:
    module.write_tsv(outdir / "redundancy.tsv",
                     ["scaffold", "length", "class", "chromosome", "target_start",
                      "target_end", "identity", "mapq", "coverage_frac"],
                     ((s, lengths[s], *c) for s, c in classes.items()))


if __name__ == "__main__":
    main()
