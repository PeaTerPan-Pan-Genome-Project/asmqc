"""Module 8: QV, k-mer completeness, spectra-cn (SPEC §8.8).

Reads: Illumina, else HiFi; never ONT. If the reads were used to build the
assembly, the QV is a self-consistency QV (INFO reads_not_independent).

Usage: python -m asmqc.m08_merqury --merqury-dir DIR --prefix P --read-type T
           --reads-used-in-assembly yes|no --outdir DIR --work DIR
"""

import argparse
import gzip
import shutil
from itertools import pairwise
from pathlib import Path

from asmqc import flags as fl
from asmqc import module
from asmqc.params import PARAMS

P = PARAMS["m08"]
M = "m08"


def coverage_peak(hist: list[tuple[int, int]]) -> int:
    """Multiplicity of the first major peak after the error trough.

    The trough is the first multiplicity whose count rises again; the peak is
    the highest count beyond it. 0 when there is no trough (no usable peak).
    """
    h = dict(hist)
    ms = sorted(h)
    trough = next((m for m, nxt in pairwise(ms) if h[nxt] > h[m]), None)
    if trough is None:
        return 0
    return max((m for m in ms if m > trough), key=lambda m: (h[m], -m))


def read_hist(path: Path) -> list[tuple[int, int]]:
    rows = []
    for line in path.read_text().splitlines():
        f = line.split()
        if len(f) >= 2 and f[0].isdigit():
            rows.append((int(f[0]), int(f[1])))
    return rows


def parse(mq: Path, prefix: str) -> dict:
    asm, _only, _total, qv, err = (mq / f"{prefix}.qv").read_text().split()[:5]
    comp = (mq / f"{prefix}.completeness.stats").read_text().split()
    return {"qv": float(qv), "error_rate": float(err), "completeness_pct": float(comp[4]),
            "asm": asm}


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--merqury-dir", type=Path, required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--read-type", choices=("illumina", "hifi"), required=True)
    ap.add_argument("--reads-used-in-assembly", choices=("yes", "no"), required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    args = ap.parse_args(argv)

    mq, pre, out = args.merqury_dir, args.prefix, args.outdir
    res = parse(mq, pre)
    cov = coverage_peak(read_hist(mq / "reads.hist"))
    independent = "no" if args.reads_used_in_assembly == "yes" else "yes"
    flags = []
    if cov < P["low_coverage"]:
        flags.append(fl.make(M, "low_coverage", value=cov,
                             message=f"k-mer coverage {cov}x < {P['low_coverage']}x"))
    if independent == "no":
        flags.append(fl.make(M, "reads_not_independent", value=args.read_type,
                             message="reads were used to build the assembly: "
                                     "self-consistency QV"))

    out.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(mq / f"{pre}.qv", out / "merqury.qv")
    shutil.copyfile(mq / f"{pre}.completeness.stats", out / "completeness.stats")
    module.write_tsv(out / "per_chromosome_qv.tsv",
                     ["seq_id", "asm_only_kmers", "total_kmers", "qv", "error_rate"],
                     (line.split("\t") for line in
                      (mq / f"{pre}.{res['asm']}.qv").read_text().splitlines()))
    with (mq / "asm_only.bed").open("rb") as src, \
            gzip.GzipFile(out / "asm_only_kmers.bed.gz", "wb", mtime=0) as dst:
        shutil.copyfileobj(src, dst)
    shutil.copyfile(mq / f"{pre}.{res['asm']}.spectra-cn.fl.png", out / "spectra-cn.png")
    shutil.copyfile(mq / f"{pre}.spectra-asm.fl.png", out / "spectra-asm.png")

    values = {
        "m08_read_type": args.read_type,
        "m08_reads_independent": independent,
        "m08_k": P["k"],
        "m08_kmer_coverage": cov,
        "m08_low_coverage": "yes" if cov < P["low_coverage"] else "no",
        "m08_qv": res["qv"],
        "m08_error_rate": res["error_rate"],
        "m08_completeness_pct": res["completeness_pct"],
    }
    module.finish(args.work, M, values, flags)


if __name__ == "__main__":
    main()
