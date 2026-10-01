"""Module 2: contiguity and anchoring (SPEC §8.2).

Contigs are the AGP components when the AGP is consistent with the FASTA,
otherwise the pieces between N-runs >= 10 bp. The N-split numbers must agree
with QUAST (--split-scaffolds); a disagreement is a bug and fails the module.

Usage: python -m asmqc.m02_contiguity --scan DIR --quast-report TSV [--agp AGP]
           --outdir RESULT_DIR --work WORK_DIR
"""

import argparse
import shutil
from pathlib import Path

from asmqc import agp, module
from asmqc.scan import read_gaps, read_sequences
from asmqc.stats import nx, pieces_between
from asmqc.validate import CHROMOSOMES

M = "m02"


def parse_quast(path: Path) -> dict[str, tuple[str, str]]:
    """report.tsv rows -> (assembly value, broken-scaffolds value)."""
    out = {}
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 3:
            out[f[0]] = (f[1], f[2])
    return out


def quast_disagreements(quast: dict, scaffolds: list[int], nsplit: list[int]) -> list[str]:
    problems = []
    for col, lengths, label in ((0, scaffolds, "scaffold"), (1, nsplit, "nsplit10 contig")):
        n50, l50 = nx(lengths, 0.5)
        n90, _ = nx(lengths, 0.9)
        for key, ours in (("N50", n50), ("L50", l50), ("N90", n90)):
            theirs = quast.get(key, ("NA", "NA"))[col]
            if theirs != str(ours):
                problems.append(f"{label} {key}: asmqc {ours}, QUAST {theirs}")
    return problems


def evaluate(scan_dir: Path, agp_path: Path | None) -> tuple[dict, list, list[int], list[int]]:
    seqs = read_sequences(scan_dir)
    lengths = {s.seq_id: s.length for s in seqs}
    gaps_by_seq: dict[str, list[tuple[int, int]]] = {}
    for name, start, end in read_gaps(scan_dir):
        gaps_by_seq.setdefault(name, []).append((start, end))

    scaffolds = [s.length for s in seqs]
    nsplit_by_seq = {s.seq_id: pieces_between(s.length, gaps_by_seq.get(s.seq_id, []))
                     for s in seqs}
    nsplit = [n for p in nsplit_by_seq.values() for n in p]

    rows: list[agp.Row] = []
    method = "nsplit10"
    if agp_path is not None:
        rows, malformed = agp.parse(agp_path)
        if not malformed and not agp.check(rows, lengths):
            method = "agp"

    contigs_by_seq: dict[str, list[int]] = {}
    gaps_per_seq: dict[str, list[int]] = {}
    if method == "agp":
        for r in rows:
            target = gaps_per_seq if r.is_gap else contigs_by_seq
            target.setdefault(r.obj, []).append(r.span)
    else:
        contigs_by_seq = nsplit_by_seq
        gaps_per_seq = {n: [e - s + 1 for s, e in g] for n, g in gaps_by_seq.items()}
    contigs = [n for p in contigs_by_seq.values() for n in p]
    gap_lengths = [n for g in gaps_per_seq.values() for n in g]

    total = sum(scaffolds)
    chrom_bp = sum(lengths.get(c, 0) for c in CHROMOSOMES)
    unplaced = [s.length for s in seqs if s.seq_id not in CHROMOSOMES]
    sn50, sl50 = nx(scaffolds, 0.5)
    cn50, cl50 = nx(contigs, 0.5)
    values = {
        "m02_scaffold_n50": sn50,
        "m02_scaffold_l50": sl50,
        "m02_scaffold_n90": nx(scaffolds, 0.9)[0],
        "m02_longest_bp": max(scaffolds, default=0),
        "m02_contig_method": method,
        "m02_contig_n50": cn50,
        "m02_contig_l50": cl50,
        "m02_contig_n90": nx(contigs, 0.9)[0],
        "m02_n_contigs": len(contigs),
        "m02_contig_n50_nsplit10": nx(nsplit, 0.5)[0],
        "m02_chrom_bp": chrom_bp,
        "m02_anchored_pct": 100 * chrom_bp / total if total else 0.0,
        "m02_n_unplaced": len(unplaced),
        "m02_unplaced_bp": sum(unplaced),
        "m02_n_gaps": len(gap_lengths),
        "m02_gap_bp": sum(gap_lengths),
    }
    per_chrom = []
    for c in CHROMOSOMES:
        cg = gaps_per_seq.get(c, [])
        per_chrom.append((c, lengths.get(c, 0), len(contigs_by_seq.get(c, [])), len(cg),
                          sum(cg), nx(contigs_by_seq.get(c, []), 0.5)[0]))
    return values, per_chrom, scaffolds, nsplit


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", type=Path, required=True)
    ap.add_argument("--quast-report", type=Path, required=True)
    ap.add_argument("--agp", type=Path)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    args = ap.parse_args(argv)

    values, per_chrom, scaffolds, nsplit = evaluate(args.scan, args.agp)
    problems = quast_disagreements(parse_quast(args.quast_report), scaffolds, nsplit)
    if problems:
        raise SystemExit("asmqc and QUAST disagree (a bug):\n" + "\n".join(problems))

    out = args.outdir
    module.write_tsv(out / "contiguity.tsv", ["metric", "value"],
                     [(k.removeprefix("m02_"), v) for k, v in values.items()])
    module.write_tsv(out / "per_chromosome.tsv",
                     ["chromosome", "length", "contigs", "gaps", "gap_bp", "contig_n50"],
                     per_chrom)
    shutil.copyfile(args.quast_report, out / "quast_report.tsv")
    module.finish(args.work, M, values)


if __name__ == "__main__":
    main()
