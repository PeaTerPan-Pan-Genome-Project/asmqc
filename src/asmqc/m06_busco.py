"""Module 6: gene-space completeness with BUSCO (SPEC §8.6).

Usage:
  python -m asmqc.m06_busco check-lineage --lineage-dir DIR
  python -m asmqc.m06_busco summarise --busco-dir DIR --fai FAI --outdir DIR --work DIR
"""

import argparse
import json
import shutil
from pathlib import Path

from asmqc import module, synteny
from asmqc.params import PARAMS
from asmqc.validate import CHROMOSOMES

P = PARAMS["m06"]
M = "m06"


def read_cfg(lineage_dir: Path) -> dict[str, str]:
    cfg = {}
    for line in (lineage_dir / "dataset.cfg").read_text().splitlines():
        k, sep, v = line.partition("=")
        if sep:
            cfg[k.strip()] = v.strip()
    return cfg


def check_lineage(lineage_dir: Path) -> list[str]:
    """Problems with the installed lineage; empty when it is the pinned one."""
    if not (lineage_dir / "dataset.cfg").exists():
        return [f"{lineage_dir}: no dataset.cfg"]
    cfg = read_cfg(lineage_dir)
    want = {"name": P["lineage"], "creation_date": P["lineage_date"],
            "number_of_BUSCOs": str(P["n_markers"])}
    return [f"lineage {k} is {cfg.get(k)!r}, expected {v!r}"
            for k, v in want.items() if cfg.get(k) != v]


def check_summary(s: dict) -> list[str]:
    """Problems with a short_summary JSON: wrong lineage, predictor or version."""
    lin, par = s["lineage_dataset"], s["parameters"]
    problems = []
    if lin["name"] != P["lineage"] or lin["creation_date"] != P["lineage_date"]:
        problems.append(f"lineage {lin['name']} {lin['creation_date']}")
    if int(s["results"]["n_markers"]) != P["n_markers"]:
        problems.append(f"{s['results']['n_markers']} markers, expected {P['n_markers']}")
    if par.get("gene_predictor") != "miniprot":
        problems.append(f"gene predictor {par.get('gene_predictor')}, expected miniprot")
    if s["versions"]["busco"] != "6.1.0":
        problems.append(f"BUSCO {s['versions']['busco']}, expected 6.1.0")
    return problems


def without_paths(s: dict) -> dict:
    """short_summary with path parameters reduced to their last component.

    The result directory is shared with the consortium; local paths stay out
    (SPEC §7.3).
    """
    params = {k: Path(v).name if isinstance(v, str) and v.startswith("/") else v
              for k, v in s["parameters"].items()}
    return s | {"parameters": params}


def read_full_table(path: Path) -> list[dict]:
    """Rows of full_table.tsv keyed busco_id, status, sequence (plus the rest).

    BUSCO 6.1.0 writes 10 columns for found BUSCOs (incl. OrthoDB url and
    description) and 2 for missing ones; column names come from the
    "# Busco id" header line.
    """
    keys = {"Busco id": "busco_id", "Status": "status", "Sequence": "sequence"}
    cols: list[str] = []
    rows = []
    with path.open() as fh:
        for line in fh:
            if line.startswith("# Busco id"):
                cols = [keys.get(c, c) for c in line[2:].rstrip("\n").split("\t")]
                continue
            if line.startswith("#") or not line.strip():
                continue
            if not cols:
                raise ValueError(f"{path}: no '# Busco id' header line")
            f = line.rstrip("\n").split("\t")
            rows.append(dict(zip(cols, f + [""] * (len(cols) - len(f)))))
    return rows


def derived(rows: list[dict]) -> tuple[dict, list[tuple]]:
    """Placement of Complete and Duplicated BUSCOs on chromosomes vs unplaced."""
    by_id: dict[str, dict] = {}
    for r in rows:
        if r["status"] not in ("Complete", "Duplicated"):
            continue
        d = by_id.setdefault(r["busco_id"], {"status": r["status"], "seqs": []})
        d["seqs"].append(r["sequence"])
    table, on_unplaced, dup_chrom, dup_unpl = [], 0, 0, 0
    for bid, d in sorted(by_id.items()):
        n_unpl = sum(s not in CHROMOSOMES for s in d["seqs"])
        on_unplaced += n_unpl > 0
        if d["status"] == "Duplicated":
            dup_unpl += n_unpl > 0
            dup_chrom += n_unpl == 0
        table.append((bid, d["status"], len(d["seqs"]), len(d["seqs"]) - n_unpl, n_unpl,
                      ";".join(d["seqs"])))
    counts = {"m06_complete_on_unplaced": on_unplaced, "m06_dup_both_on_chrom": dup_chrom,
              "m06_dup_any_on_unplaced": dup_unpl}
    return counts, table


def summarise(busco_dir: Path, outdir: Path, work: Path, fai: Path | None = None) -> None:
    lineage = P["lineage"]
    summary_json = busco_dir / f"short_summary.specific.{lineage}.busco.json"
    full_table = busco_dir / f"run_{lineage}" / "full_table.tsv"
    s = json.loads(summary_json.read_text())
    problems = check_summary(s)
    if problems:
        raise SystemExit("BUSCO run does not match the pinned setup: " + "; ".join(problems))
    res = s["results"]
    rows = read_full_table(full_table)
    counts, table = derived(rows)
    values = {
        "m06_complete_pct": res["Complete percentage"],
        "m06_single_pct": res["Single copy percentage"],
        "m06_duplicated_pct": res["Multi copy percentage"],
        "m06_fragmented_pct": res["Fragmented percentage"],
        "m06_missing_pct": res["Missing percentage"],
        "m06_n_markers": res["n_markers"],
        "m06_internal_stop_pct": res["internal_stop_codon_percent"],
        **counts,
        "m06_lineage": f"{lineage} {P['lineage_date']}",
    }
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "short_summary.json").write_text(json.dumps(without_paths(s), indent=4) + "\n")
    shutil.copyfile(full_table, outdir / "full_table.tsv")
    module.write_tsv(outdir / "busco_derived.tsv",
                     ["busco_id", "status", "n_copies", "n_on_chromosomes", "n_on_unplaced",
                      "sequences"], table)
    anchors = synteny.anchors_path()
    if fai is not None and anchors.exists():  # report-only synteny against Caméor v2
        _, ref = synteny.read_anchors(anchors)
        pts = synteny.points(rows, ref)
        module.write_tsv(outdir / "synteny_points.tsv", synteny.POINT_COLUMNS, pts)
        module.write_tsv(outdir / "synteny.tsv", synteny.TABLE_COLUMNS,
                         synteny.table(pts, module.read_fai(fai)))
    module.finish(work, M, values)


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=("check-lineage", "summarise"))
    ap.add_argument("--lineage-dir", type=Path)
    ap.add_argument("--busco-dir", type=Path)
    ap.add_argument("--outdir", type=Path)
    ap.add_argument("--work", type=Path)
    ap.add_argument("--fai", type=Path, help="tested assembly .fai, for the synteny table")
    args = ap.parse_args(argv)
    if args.step == "check-lineage":
        problems = check_lineage(args.lineage_dir)
        if problems:
            raise SystemExit("; ".join(problems))
        print(json.dumps(read_cfg(args.lineage_dir)))
    else:
        summarise(args.busco_dir, args.outdir, args.work, args.fai)


if __name__ == "__main__":
    main()
