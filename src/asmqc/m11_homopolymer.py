"""Module 11: homopolymer and short-STR errors (SPEC §8.11).

Steps (each a subcommand, run by workflow/rules/m11.smk):
  callable  median depth over chromosomes (1 kb sampling) -> callable.bed
  call      bcftools mpileup | call per region in parallel, concat, norm
  classify  hom-alt calls -> HP / STR2 / other indel / SNV, tables, summary

Heterozygous calls in an inbred line mostly reflect mismapping in repeats;
they are counted (m11_het_calls) but never as errors.
"""

import argparse
import bisect
import gzip
import re
import subprocess
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from asmqc import depth, module
from asmqc import flags as fl
from asmqc.params import PARAMS
from asmqc.validate import CHROMOSOMES

P = PARAMS["m11"]
M = "m11"
RUN_BINS = [(4, 8), (9, 10), (11, 12), (13, 15), (16, 20), (21, None)]
CHANGE_BINS = ["1", "2", "3", "4", ">=5"]
CHUNK_BP = 50_000_000  # unplaced sequences are called in chunks of about this size


# --- callable ------------------------------------------------------------------
def cmd_callable(args) -> None:
    lengths = module.read_fai(args.fai)
    med = depth.median_depth(args.bam, lengths, args.work)
    lo, hi = P["callable_depth_low"] * med, P["callable_depth_high"] * med
    n = depth.callable_bed(args.bam, lengths, lo, hi, args.work / "callable.bed", args.threads)
    (args.work / "callable.txt").write_text(f"median_depth\t{med}\ncallable_bp\t{n}\n")


# --- call ------------------------------------------------------------------------
def regions(lengths: dict[str, int]) -> list[list[str]]:
    """Chromosomes one by one, unplaced sequences grouped into chunks."""
    out = [[c] for c in CHROMOSOMES]
    chunk, size = [], 0
    for s, n in lengths.items():
        if s in CHROMOSOMES:
            continue
        chunk.append(s)
        size += n
        if size >= CHUNK_BP:
            out.append(chunk)
            chunk, size = [], 0
    if chunk:
        out.append(chunk)
    return out


def cmd_call(args) -> None:
    lengths = module.read_fai(args.fai)
    parts_dir = args.work / "call_parts"
    parts_dir.mkdir(parents=True, exist_ok=True)
    pacbio = "-X pacbio-ccs" if args.read_type == "hifi" else ""

    def one(i_seqs):
        i, seqs = i_seqs
        bed = parts_dir / f"{i}.regions.bed"
        bed.write_text("".join(f"{s}\t0\t{lengths[s]}\n" for s in seqs))
        vcf = parts_dir / f"{i}.vcf.gz"
        cmd = (f"bcftools mpileup -f {args.fa} -a AD,DP -q {P['min_mapq']} -Q {P['min_baseq']}"
               f" -d {P['max_depth']} {pacbio} -T {args.work}/callable.bed -R {bed} {args.bam}"
               f" 2>> {parts_dir}/{i}.log"
               f" | bcftools call -m -v --ploidy 2 -Oz -o {vcf} 2>> {parts_dir}/{i}.log")
        subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)
        return vcf

    with ThreadPoolExecutor(max_workers=max(1, args.threads)) as pool:
        vcfs = list(pool.map(one, enumerate(regions(lengths))))
    out = args.work / "calls.norm.vcf.gz"
    cmd = (f"bcftools concat -Ou {' '.join(map(str, vcfs))}"
           f" | bcftools norm -f {args.fa} -m -both -Oz -o {out} && bcftools index -c {out}")  # CSI: TBI has the 2^29 limit
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)


# --- classify ----------------------------------------------------------------------
@dataclass
class Call:
    seq: str
    pos: int  # 1-based VCF POS (anchor base for indels)
    ref: str
    alt: str
    qual: float
    gt: str
    dp: int


def read_calls(vcf: Path) -> tuple[list[Call], list[str]]:
    calls, header = [], []
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line)
                continue
            f = line.rstrip("\n").split("\t")
            fmt, sample = f[8].split(":"), f[9].split(":")
            gt = sample[fmt.index("GT")].replace("|", "/")
            dp = re.search(r"(?:^|;)DP=(\d+)", f[7])
            qual = float(f[5]) if f[5] != "." else 0.0
            calls.append(Call(f[0], int(f[1]), f[3].upper(), f[4].upper(), qual, gt,
                              int(dp[1]) if dp else 0))
    return calls, header


def run_length(seq: str, start: int, base: str) -> int:
    """Length of the run of `base` starting at 0-based `start`."""
    n = 0
    while start + n < len(seq) and seq[start + n] == base:
        n += 1
    return n


def tandem_copies(seq: str, start: int, unit: str) -> int:
    n = 0
    while seq[start + 2 * n:start + 2 * n + 2] == unit:
        n += 1
    return n


def classify_call(c: Call, seq: str) -> dict:
    """Class of one hom-alt call; `seq` is the uppercase assembly sequence."""
    if len(c.ref) == len(c.alt):
        return {"class": "snv"}
    ins = len(c.alt) > len(c.ref)
    indel = c.alt[1:] if ins else c.ref[1:]
    k = len(indel)
    after = c.pos  # 0-based index of the base right after the anchor
    if len(set(indel)) == 1:
        b = indel[0]
        run = run_length(seq, after, b)
        change = k if ins else -k
        if max(run, run + change) >= P["hp_min_len"]:
            return {"class": "hp", "base": b, "run": run, "size": k,
                    "dir": "ins" if ins else "del"}
    if k % 2 == 0 and indel[0] != indel[1] and indel == indel[:2] * (k // 2):
        unit = indel[:2]
        copies = tandem_copies(seq, after, unit)
        if copies >= P["str2_min_copies"]:
            return {"class": "str2", "unit": unit, "copies": copies, "size": k // 2,
                    "dir": "ins" if ins else "del"}
    return {"class": "other_indel"}


def run_bin(n: int) -> str:
    for lo, hi in RUN_BINS:
        if hi is None:
            return f">{lo - 1}"
        if lo <= n <= hi:
            return f"{lo}-{hi}"
    return f"<{RUN_BINS[0][0]}"


def change_bin(k: int) -> str:
    return str(k) if k < 5 else ">=5"


def read_fasta(path: Path):
    name, chunks = None, []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks).upper()


def read_bed(path: Path) -> dict[str, tuple[list[int], list[int]]]:
    """seq -> (starts, ends) of sorted, non-overlapping intervals."""
    out: dict[str, tuple[list[int], list[int]]] = {}
    for line in path.read_text().splitlines():
        s, a, b = line.split("\t")[:3]
        starts, ends = out.setdefault(s, ([], []))
        starts.append(int(a))
        ends.append(int(b))
    return out


def inside(bed, seq: str, start: int, end: int) -> bool:
    """Is the half-open [start, end) entirely within one callable interval?"""
    if seq not in bed:
        return False
    starts, ends = bed[seq]
    i = bisect.bisect_right(starts, start) - 1
    return i >= 0 and end <= ends[i]


HP_RE = re.compile("|".join(f"{b}{{{P['hp_min_len']},}}" for b in "ACGT"))


def cmd_classify(args) -> None:
    calls, header = read_calls(args.work / "calls.norm.vcf.gz")
    bed = read_bed(args.work / "callable.bed")
    callable_bp = sum(sum(e) - sum(s) for s, e in bed.values())
    by_seq: dict[str, list[Call]] = {}
    for c in calls:
        if c.qual >= P["min_qual"]:
            by_seq.setdefault(c.seq, []).append(c)

    errors, het, hp_runs = [], 0, 0
    for name, seq in read_fasta(args.fa):
        hp_runs += sum(inside(bed, name, m.start(), m.end()) for m in HP_RE.finditer(seq))
        for c in by_seq.get(name, []):
            if c.gt in ("0/1", "1/0"):
                het += 1
            elif c.gt == "1/1" and inside(bed, c.seq, c.pos - 1, c.pos):
                errors.append((c, classify_call(c, seq)))

    counts = Counter(e["class"] for _, e in errors)
    hp = [e for _, e in errors if e["class"] == "hp"]
    str2 = [e for _, e in errors if e["class"] == "str2"]
    n_ins = sum(e["dir"] == "ins" for e in hp)
    n_del = len(hp) - n_ins
    mb = callable_bp / 1e6
    values = {
        "m11_read_type": args.read_type,
        "m11_reads_independent": "no" if args.reads_used_in_assembly == "yes" else "yes",
        "m11_callable_bp": callable_bp,
        "m11_hp_errors": len(hp),
        "m11_hp_errors_per_mb": len(hp) / mb if mb else None,
        "m11_hp_errors_per_10k_runs": 1e4 * len(hp) / hp_runs if hp_runs else None,
        "m11_dinuc_errors": len(str2),
        "m11_dinuc_errors_per_mb": len(str2) / mb if mb else None,
        "m11_other_indels_per_mb": counts["other_indel"] / mb if mb else None,
        "m11_snv_per_mb": counts["snv"] / mb if mb else None,
        "m11_hp_pct_of_errors": 100 * len(hp) / len(errors) if errors else None,
        "m11_hp_ins_del_ratio": n_ins / n_del if n_del else None,
        "m11_hp_at_pct": 100 * sum(e["base"] in "AT" for e in hp) / len(hp) if hp else None,
        "m11_het_calls": het,
    }
    flags = []
    if values["m11_reads_independent"] == "no":
        flags.append(fl.make(M, "reads_not_independent", value=args.read_type,
                             message="reads were used to build the assembly"))
    write_outputs(args, errors, header, hp_runs, callable_bp)
    module.finish(args.work, M, values, flags)


def write_outputs(args, errors, header, hp_runs, callable_bp) -> None:
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    table = Counter()
    for _, e in errors:
        if e["class"] == "hp":
            table[("hp", run_bin(e["run"]), change_bin(e["size"]),
                   "AT" if e["base"] in "AT" else "GC", e["dir"])] += 1
        elif e["class"] == "str2":
            table[("str2", e["unit"], change_bin(e["size"]), "", e["dir"])] += 1
    rows = [(*k, n) for k, n in sorted(table.items())]
    rows += [("denominator", "callable_bp", "", "", "", callable_bp),
             ("denominator", "hp_runs_ge4_in_callable", "", "", "", hp_runs)]
    module.write_tsv(out / "errors.tsv",
                     ["class", "run_length_or_unit", "change_size", "base", "direction",
                      "count"], rows)
    with gzip.GzipFile(out / "errors.bed.gz", "wb", mtime=0) as gz:
        for c, e in errors:
            if e["class"] not in ("hp", "str2"):
                continue
            k = len(c.alt) - len(c.ref)
            motif = e.get("base") or e.get("unit")
            change = f"{'+' if k > 0 else '-'}{abs(k)}{motif if e['class'] == 'hp' else ''}"
            run = e.get("run", e.get("copies"))
            end = c.pos + (abs(k) if k < 0 else 0)
            gz.write(f"{c.seq}\t{c.pos}\t{end}\t{e['class']}\t{motif}\t{run}\t{change}"
                     f"\t{c.qual:g}\t{c.dp}\n".encode())
    with gzip.GzipFile(out / "hom_calls.vcf.gz", "wb", mtime=0) as gz:
        # header lines with local paths (reference, command lines) are dropped:
        # this file is sent to the consortium (SPEC §7.3, §10)
        gz.write("".join(h for h in header
                         if not h.startswith(("##reference=", "##bcftools"))).encode())
        for c, _ in errors:
            gz.write(f"{c.seq}\t{c.pos}\t.\t{c.ref}\t{c.alt}\t{c.qual:g}\t.\tDP={c.dp}"
                     f"\tGT\t1/1\n".encode())


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=("callable", "call", "classify"))
    ap.add_argument("--bam", type=Path)
    ap.add_argument("--fa", type=Path)
    ap.add_argument("--fai", type=Path)
    ap.add_argument("--read-type", choices=("illumina", "hifi"))
    ap.add_argument("--reads-used-in-assembly", choices=("yes", "no"))
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--outdir", type=Path)
    ap.add_argument("--work", type=Path, required=True)
    args = ap.parse_args(argv)
    args.work.mkdir(parents=True, exist_ok=True)
    {"callable": cmd_callable, "call": cmd_call, "classify": cmd_classify}[args.step](args)


if __name__ == "__main__":
    main()
