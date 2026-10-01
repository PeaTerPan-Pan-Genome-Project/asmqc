#!/usr/bin/env python3
"""Build a ~50 Mb real-data fixture from a full pea assembly and its asmqc result.

The fixture keeps, per chromosome, both ends (telomeres) and blocks around
the rDNA arrays and interstitial telomeric arrays found by asmqc on the full
assembly, joined by 100 N gaps; plus a selection of unplaced scaffolds
(organelle-like, rDNA-only, telomere-carrying, duplicates whose source lies in
the fixture, then a deterministic sample). Sequence names stay as in the
source FASTA, so runs need the same --chromosomes mapping.

Written for Cameor v2 (GCA_977071245.1); nothing here is specific to it.
Only published data go in; the outputs are never committed.

Usage:
  make_cameor_fixture.py --fasta FULL.fa --result RESULT_DIR --chromosomes MAP --out DIR
      [--end-bp 3000137] [--flank-bp 100000] [--max-block-bp 4000000]
      [--min-array-copies 10] [--unplaced-bp 5000000]

  FULL.fa       uncompressed source FASTA with a samtools .fai
  RESULT_DIR    asmqc result of FULL.fa (<outdir>/<label>/), modules 4, 5 and 7
  MAP           the --chromosomes string of that run (old=chr1,...)

Writes to DIR:
  fixture.fa         fixture assembly
  fixture.agp        AGP 2.0; components named <source>_subseq_<start>:<end>
  regions.bed        source intervals kept (0-based), for read extraction
  chromosomes.txt    the --chromosomes value for asmqc runs on the fixture
  fixture.json       parameters, blocks and totals
"""

import argparse
import csv
import json
import subprocess
import zlib
from pathlib import Path

GAP = 100
WIDTH = 60


def read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def merge(blocks: list[tuple[int, int, str]]) -> list[tuple[int, int, str]]:
    """Merge overlapping 1-based inclusive blocks; reasons are joined."""
    out: list[list] = []
    for s, e, why in sorted(blocks):
        if out and s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
            out[-1][2] += f"+{why}"
        else:
            out.append([s, e, why])
    return [tuple(b) for b in out]


def off_grid(s: int, e: int, length: int) -> tuple[int, int]:
    """Move block ends off multiples of 1,000: tool windows (tidk, 10 kb) would
    otherwise make round AGP cuts that M1 reports as curation artefacts."""
    if s > 1 and (s - 1) % 1000 == 0:
        s += 137
    if e < length and e % 1000 == 0:
        e -= 137
    return s, e


def chromosome_blocks(chrom: str, length: int, res: Path, args) -> list[tuple[int, int, str]]:
    end = min(args.end_bp, length // 2)
    blocks = [(1, end, "start"), (length - end + 1, length, "end")]
    # Arrays below min_array_copies are mostly dispersed subunit fragments
    # (BLAST hits of any length count as copies in M5).
    for a in read_tsv(res / "m05_organelle_rdna" / "rdna_arrays.tsv"):
        if a["seq_id"] == chrom and int(a["copies"]) >= args.min_array_copies:
            s = max(1, int(a["start"]) - args.flank_bp)
            e = min(length, int(a["end"]) + args.flank_bp, s + args.max_block_bp - 1)
            blocks.append((*off_grid(s, e, length), f"rdna{a['family']}"))
    for t in read_tsv(res / "m04_telomeres" / "interstitial.tsv"):
        if t["chromosome"] == chrom:
            s = max(1, int(t["start"]) - 50_000)
            e = min(length, int(t["end"]) + 50_000)
            blocks.append((*off_grid(s, e, length), "interstitial_telomere"))
    return merge(blocks)


def unplaced_selection(unplaced: dict[str, int], res: Path, kept: dict[str, list], args):
    """Must-include scaffolds first (by reason), then a deterministic sample."""
    def inside(chrom, s, e):
        return any(a <= s and e <= b for a, b, _ in kept.get(chrom, []))

    reasons: dict[str, str] = {}
    for r in read_tsv(res / "m05_organelle_rdna" / "organelle_scaffolds.tsv"):
        if r["organelle"] != "none" and r["seq_id"] in unplaced:
            reasons.setdefault(r["seq_id"], f"organelle_{r['organelle']}")
    for a in read_tsv(res / "m05_organelle_rdna" / "rdna_arrays.tsv"):
        if a["context"] == "rdna_only_scaffold" and a["seq_id"] in unplaced:
            reasons.setdefault(a["seq_id"], f"rdna_only_{a['family']}")
    for t in read_tsv(res / "m04_telomeres" / "telomeres_unplaced.tsv"):
        if t["seq_id"] in unplaced:
            reasons.setdefault(t["seq_id"], "telomere")
    for r in read_tsv(res / "m07_redundancy" / "redundancy.tsv"):
        if (r["class"] in ("duplicate", "partial_overlap") and r["scaffold"] in unplaced
                and inside(r["chromosome"], int(r["target_start"]), int(r["target_end"]))):
            reasons.setdefault(r["scaffold"], r["class"])

    chosen, total = {}, 0
    budget = args.unplaced_bp
    for name in sorted(reasons, key=lambda n: (reasons[n], unplaced[n], n)):
        if total + unplaced[name] <= budget or not chosen:
            chosen[name] = reasons[name]
            total += unplaced[name]
    for name in sorted(unplaced, key=lambda n: zlib.crc32(n.encode())):
        if name in chosen or total + unplaced[name] > budget:
            continue
        chosen[name] = "sample"
        total += unplaced[name]
    return chosen


def fetch(fasta: Path, seq: str, start: int, end: int) -> str:
    out = subprocess.run(["samtools", "faidx", str(fasta), f"{seq}:{start}-{end}"],
                         capture_output=True, text=True, check=True).stdout
    return "".join(out.splitlines()[1:])


def write_record(fh, name: str, seq: str) -> None:
    fh.write(f">{name}\n")
    for i in range(0, len(seq), WIDTH):
        fh.write(seq[i:i + WIDTH] + "\n")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--fasta", type=Path, required=True)
    ap.add_argument("--result", type=Path, required=True)
    ap.add_argument("--chromosomes", required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--end-bp", type=int, default=3_000_137)
    ap.add_argument("--flank-bp", type=int, default=100_000)
    ap.add_argument("--min-array-copies", type=int, default=10)
    ap.add_argument("--max-block-bp", type=int, default=4_000_000)
    ap.add_argument("--unplaced-bp", type=int, default=5_000_000)
    args = ap.parse_args(argv)

    fai = Path(f"{args.fasta}.fai")
    if not fai.exists():
        subprocess.run(["samtools", "faidx", str(args.fasta)], check=True)
    lengths = {f[0]: int(f[1]) for f in (line.split("\t") for line in fai.read_text().splitlines())}
    to_chr = dict(item.split("=", 1) for item in args.chromosomes.split(","))
    from_chr = {v: k for k, v in to_chr.items()}
    unplaced = {n: L for n, L in lengths.items() if n not in to_chr}

    kept = {f"chr{i}": chromosome_blocks(f"chr{i}", lengths[from_chr[f"chr{i}"]], args.result,
                                         args)
            for i in range(1, 8)}
    chosen = unplaced_selection(unplaced, args.result, kept, args)

    args.out.mkdir(parents=True, exist_ok=True)
    agp, bed, blocks_out = [], [], []
    with (args.out / "fixture.fa").open("w") as fa:
        for i in range(1, 8):
            chrom = f"chr{i}"
            src = from_chr[chrom]
            parts, pos, part_no = [], 0, 0
            for j, (s, e, why) in enumerate(kept[chrom]):
                if j:
                    parts.append("N" * GAP)
                    part_no += 1
                    agp.append([src, pos + 1, pos + GAP, part_no, "N", GAP, "contig", "no", "na"])
                    pos += GAP
                seq = fetch(args.fasta, src, s, e)
                parts.append(seq)
                part_no += 1
                agp.append([src, pos + 1, pos + len(seq), part_no, "W",
                            f"{src}_subseq_{s}:{e}", 1, len(seq), "+"])
                bed.append([src, s - 1, e, f"{chrom}:{why}"])
                blocks_out.append({"source": src, "chromosome": chrom, "start": s, "end": e,
                                   "reason": why, "fixture_start": pos + 1})
                pos += len(seq)
            write_record(fa, src, "".join(parts))
        for name, why in chosen.items():
            seq = fetch(args.fasta, name, 1, unplaced[name])
            write_record(fa, name, seq)
            agp.append([name, 1, len(seq), 1, "W", name, 1, len(seq), "+"])
            bed.append([name, 0, len(seq), f"unplaced:{why}"])
            blocks_out.append({"source": name, "chromosome": None, "start": 1,
                               "end": len(seq), "reason": why, "fixture_start": 1})

    with (args.out / "fixture.agp").open("w") as fh:
        fh.write("##agp-version\t2.0\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in agp)
    with (args.out / "regions.bed").open("w") as fh:
        fh.writelines("\t".join(map(str, r)) + "\n" for r in bed)
    (args.out / "chromosomes.txt").write_text(args.chromosomes + "\n")
    chrom_bp = sum(e - s + 1 for b in kept.values() for s, e, _ in b)
    unpl_bp = sum(unplaced[n] for n in chosen)
    (args.out / "fixture.json").write_text(json.dumps({
        "source_fasta": args.fasta.name,
        "parameters": {k: v for k, v in vars(args).items()
                       if k.endswith(("_bp", "_copies"))},
        "chromosome_bp": chrom_bp, "unplaced_bp": unpl_bp, "n_unplaced": len(chosen),
        "unplaced_reasons": {r: sum(1 for v in chosen.values() if v == r)
                             for r in sorted(set(chosen.values()))},
        "blocks": blocks_out,
    }, indent=1) + "\n")
    print(f"fixture: {chrom_bp:,} bp in chr1-chr7 blocks, {unpl_bp:,} bp in "
          f"{len(chosen)} unplaced scaffolds -> {args.out}")


if __name__ == "__main__":
    main()
