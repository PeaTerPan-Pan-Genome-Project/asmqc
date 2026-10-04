"""Module 9: read-back structural validation with CRAQ (SPEC §8.9).

CRAQ 1.10 run on the shared BAMs writes (relative to its run directory):
  output/runAQI_out/out_final.Report          per-sequence and Genome rows:
      #Chr Covered.Rate Low-confident.Rate Avg.CRH Avg.CSH Avg.CRE(R-AQI) Avg.CSE(S-AQI)
  output/runAQI_out/locER_out/out_final.CRE.bed
  output/runAQI_out/strER_out/out_final.CSE.bed

Usage: python -m asmqc.m09_craq --craq-dir DIR --bam LONG_BAM --fai FAI [--agp AGP]
           --long-read-type T --short-reads-used yes|no --outdir DIR --work DIR
"""

import argparse
import re
import shutil
from pathlib import Path

from asmqc import agp, depth, module
from asmqc import flags as fl
from asmqc.params import PARAMS

P = PARAMS["m09"]
M = "m09"
AQI_RE = re.compile(r"^([\d.eE+-]+)\(([\d.eE+-]+)\)$")


def parse_report(path: Path) -> dict[str, float]:
    """R-AQI and S-AQI from the Genome row of out_final.Report."""
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if f[0] == "Genome":
            r = AQI_RE.match(f[5])
            s = AQI_RE.match(f[6])
            if not (r and s):
                raise ValueError(f"{path}: unexpected Genome row {line!r}")
            return {"r_aqi": float(r[2]), "s_aqi": float(s[2]),
                    "avg_cre_per_mb": float(r[1]), "avg_cse_per_mb": float(s[1])}
    raise ValueError(f"{path}: no Genome row")


def read_bed(path: Path) -> list[tuple[str, int, int]]:
    out = []
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 3:
            out.append((f[0], int(f[1]), int(f[2])))
    return out


def junctions(rows: list[agp.Row]) -> dict[str, list[tuple[int, int]]]:
    """AGP gaps (start, end) and W-W junctions (pos, pos) per object, 1-based."""
    out: dict[str, list[tuple[int, int]]] = {}
    prev = None
    for r in rows:
        js = out.setdefault(r.obj, [])
        if r.is_gap:
            js.append((r.obj_beg, r.obj_end))
        elif prev is not None and prev.obj == r.obj and not prev.is_gap:
            js.append((prev.obj_end, r.obj_beg))
        prev = r
    return out


def junction_crosscheck(cse: list[tuple[str, int, int]], js: dict, window: int):
    """Per CSE: distance to the nearest gap or W-W junction (None if none)."""
    out = []
    for seq, s, e in cse:
        d = min((max(0, a - e, s - b) for a, b in js.get(seq, [])), default=None)
        out.append((seq, s, e, "NA" if d is None else d,
                    "yes" if d is not None and d <= window else "no"))
    return out


def sorted_report(text: str, order: list[str]) -> str:
    """out_final.Report with the per-sequence rows in assembly order. CRAQ
    writes them in Perl hash order, which differs between runs; the header
    and Genome rows stay first."""
    lines = text.splitlines()
    rank = {s: i for i, s in enumerate(order)}
    head = [ln for ln in lines if ln.split("\t")[0] not in rank]
    rows = sorted((ln for ln in lines if ln.split("\t")[0] in rank),
                  key=lambda ln: rank[ln.split("\t")[0]])
    return "\n".join(head + rows) + "\n"


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--craq-dir", type=Path, required=True)
    ap.add_argument("--bam", type=Path, required=True)
    ap.add_argument("--fai", type=Path, required=True)
    ap.add_argument("--agp", type=Path)
    ap.add_argument("--long-read-type", choices=("hifi", "ont_r9", "ont_r10"), required=True)
    ap.add_argument("--short-reads-used", choices=("yes", "no"), required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    args = ap.parse_args(argv)

    run = args.craq_dir / "output" / "runAQI_out"
    report = parse_report(run / "out_final.Report")
    cre = read_bed(run / "locER_out" / "out_final.CRE.bed")
    cse = read_bed(run / "strER_out" / "out_final.CSE.bed")
    lengths = module.read_fai(args.fai)
    cov = depth.median_depth(args.bam, lengths, args.work)

    near = inside = None
    rows: list[tuple] = [(s, a, b, "NA", "NA") for s, a, b in cse]
    if args.agp:
        agp_rows, malformed = agp.parse(args.agp)
        if not malformed and not agp.check(agp_rows, lengths):
            rows = junction_crosscheck(cse, junctions(agp_rows), P["junction_window_bp"])
            near = sum(r[4] == "yes" for r in rows)
            inside = len(rows) - near

    flags = []
    if cov < P["low_coverage"]:
        flags.append(fl.make(M, "low_coverage", value=f"{cov:g}",
                             message=f"median long-read depth {cov:g}x < {P['low_coverage']}x"))
    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)
    (out / "out_final.Report").write_text(
        sorted_report((run / "out_final.Report").read_text(), list(module.read_fai(args.fai))))
    shutil.copyfile(run / "locER_out" / "out_final.CRE.bed", out / "CRE.bed")
    shutil.copyfile(run / "strER_out" / "out_final.CSE.bed", out / "CSE.bed")
    module.write_tsv(out / "craq_derived.tsv",
                     ["seq_id", "start", "end", "distance_to_agp_junction_bp",
                      "near_agp_junction"], rows)
    values = {
        "m09_long_read_type": args.long_read_type,
        "m09_short_reads_used": args.short_reads_used,
        "m09_coverage": cov,
        "m09_low_coverage": "yes" if cov < P["low_coverage"] else "no",
        "m09_aqi": None,  # CRAQ 1.10 reports no single AQI; see IMPLEMENTATION_PLAN §6
        "m09_r_aqi": report["r_aqi"],
        "m09_s_aqi": report["s_aqi"],
        "m09_n_cre": len(cre),
        "m09_n_cse": len(cse),
        "m09_cse_near_agp_junction": near,
        "m09_cse_inside_contig": inside,
    }
    module.finish(args.work, M, values, flags)


if __name__ == "__main__":
    main()
