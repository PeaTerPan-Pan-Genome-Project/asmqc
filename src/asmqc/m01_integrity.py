"""Module 1: integrity, format and ENA rules (SPEC §8.1).

Interprets the shared scan (asmqc.scan) and checks the AGP.

Usage: python -m asmqc.m01_integrity --scan DIR --md5 FILE --assembly-name NAME
           [--agp AGP --agp-original AGP] --outdir RESULT_DIR --work WORK_DIR
"""

import argparse
import json
from dataclasses import astuple
from pathlib import Path

from asmqc import agp, module
from asmqc import flags as fl
from asmqc.params import PARAMS
from asmqc.scan import SEQ_COLUMNS, read_gaps, read_sequences
from asmqc.seqio import file_md5

P = PARAMS["m01"]
M = "m01"


def evaluate(scan_dir: Path, agp_path: Path | None) -> tuple[dict, list[fl.Flag], dict]:
    """Summary values, flags and extra tables (AGP problems) for module 1."""
    seqs = read_sequences(scan_dir)
    gaps = read_gaps(scan_dir)
    facts = json.loads((scan_dir / "file.json").read_text())
    flags: list[fl.Flag] = []

    def flag(code, seq_id="", start="", end="", value="", message=""):
        flags.append(fl.make(M, code, seq_id, start, end, value, message))

    total = sum(s.length for s in seqs)
    for s in seqs:
        if s.length == 0:
            flag("empty_sequence", s.seq_id, message="sequence has no bases")
        elif s.length < P["min_len_ena"]:
            flag("seq_lt_20bp", s.seq_id, value=s.length,
                 message=f"{s.length} bp; ENA requires >= {P['min_len_ena']} bp")
        elif s.length < P["min_len_warning"]:
            flag("seq_lt_200bp", s.seq_id, value=s.length, message=f"{s.length} bp")
        if s.length and (s.leading_n or s.trailing_n):
            if s.leading_n:
                flag("terminal_n", s.seq_id, 1, s.leading_n, s.leading_n, "leading N")
            if s.trailing_n:
                flag("terminal_n", s.seq_id, s.length - s.trailing_n + 1, s.length,
                     s.trailing_n, "trailing N")
        if s.other:
            flag("invalid_char", s.seq_id, value=s.other,
                 message=f"{s.other} non-IUPAC character(s)")
        if s.iupac:
            flag("iupac_present", s.seq_id, value=s.iupac,
                 message=f"{s.iupac} IUPAC ambiguity code(s)")
        if s.length and s.n > P["n_fraction_warning"] * s.length:
            flag("n_fraction_gt_50pct", s.seq_id, value=f"{s.n / s.length:.4f}",
                 message=f"{100 * s.n / s.length:.1f} % N")
    for name, count in facts["duplicate_names"].items():
        flag("duplicate_name", name, value=count, message=f"name occurs {count} times")
    for name in facts["bad_names"]:
        flag("invalid_char", name, message="name empty or with non-printable characters")
    if facts["crlf"]:
        flag("crlf_line_endings", message="CRLF line endings")
    widths = {s.line_width for s in seqs if s.line_width}
    if len(widths) > 1 or not all(s.width_consistent for s in seqs):
        flag("inconsistent_line_width", value=",".join(map(str, sorted(widths))),
             message="line width differs between or within sequences")

    unplaced = [s for s in seqs if not module.is_chromosome(s.seq_id)]
    unit = P["round_unit_bp"]
    for s in unplaced:
        if s.length and s.length % unit == 0:
            flag("round_length", s.seq_id, value=s.length,
                 message=f"length is a multiple of {unit:,}")

    agp_consistent, round_cuts, n_cuts, agp_problems = None, None, None, []
    if agp_path is not None:
        rows, malformed = agp.parse(agp_path)
        lengths = {s.seq_id: s.length for s in seqs}
        agp_problems = [("", p) for p in malformed] + agp.check(rows, lengths)
        agp_consistent = "no" if agp_problems else "yes"
        for obj, problem in agp_problems:
            flag("agp_mismatch", obj, message=problem)
        cuts = agp.cut_coordinates(rows)
        n_cuts = len(cuts)
        round_cuts = 0
        for row, coord in cuts:
            if coord % unit == 0:
                round_cuts += 1
                flag("round_agp_cut", row.obj, row.obj_beg, row.obj_end, coord,
                     f"component {row.comp}: cut at {coord:,}")

    gap_lengths = [e - s + 1 for _, s, e in gaps]
    values = {
        "m01_n_seq": len(seqs),
        "m01_total_bp": total,
        "m01_n_lt20bp": sum(s.length < P["min_len_ena"] for s in seqs),
        "m01_n_lt200bp": sum(s.length < P["min_len_warning"] for s in seqs),
        "m01_n_terminal_n": sum(bool(s.length and (s.leading_n or s.trailing_n)) for s in seqs),
        "m01_n_duplicate_names": len(facts["duplicate_names"]),
        "m01_n_invalid_chars": sum(s.other for s in seqs),
        "m01_n_iupac": sum(s.iupac for s in seqs),
        "m01_n_gaps_ge10": len(gap_lengths),
        "m01_gap_bp_ge10": sum(gap_lengths),
        "m01_n_gaps_ge100": sum(n >= P["gap_long_bp"] for n in gap_lengths),
        "m01_n_seq_gt50pct_n": sum(bool(s.length) and s.n > P["n_fraction_warning"] * s.length
                                   for s in seqs),
        "m01_softmask_pct": 100 * sum(s.lower for s in seqs) / total if total else 0.0,
        "m01_agp_consistent": agp_consistent,
        "m01_n_round_lengths": sum(bool(s.length) and s.length % unit == 0 for s in unplaced),
        "m01_round_lengths_expected": len(unplaced) / unit,
        "m01_n_round_agp_cuts": round_cuts,
        "m01_n_agp_cuts": n_cuts,
    }
    return values, flags, {"agp_problems": agp_problems, "seqs": seqs, "gaps": gaps}


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan", type=Path, required=True)
    ap.add_argument("--md5", type=Path, required=True)
    ap.add_argument("--assembly-name", required=True)
    ap.add_argument("--agp", type=Path)
    ap.add_argument("--agp-original", type=Path)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--work", type=Path, required=True)
    args = ap.parse_args(argv)

    values, flags, extra = evaluate(args.scan, args.agp)
    out = args.outdir
    module.write_tsv(out / "integrity.tsv", ["metric", "value"],
                     [(k.removeprefix("m01_"), "NA" if v is None else v)
                      for k, v in values.items()])
    module.write_tsv(out / "sequences.tsv", SEQ_COLUMNS, (astuple(s) for s in extra["seqs"]))
    module.write_tsv(out / "gaps.tsv", ["seq_id", "start", "end", "length"],
                     ((n, s, e, e - s + 1) for n, s, e in extra["gaps"]))
    lines = [f"{args.md5.read_text().strip()}  {args.assembly_name}"]
    if args.agp_original:
        lines.append(f"{file_md5(args.agp_original)}  {args.agp_original.name}")
    (out / "checksums.md5").write_text("\n".join(lines) + "\n")
    module.finish(args.work, M, values, flags)


if __name__ == "__main__":
    main()
