"""Merge module results into qc_summary.tsv, flags.tsv and run_manifest.json.

Runs in the wrapper after Snakemake returns, not as a rule, so the outputs
exist even when a module failed (IMPLEMENTATION_PLAN §2).

Module contract: a module that finished writes
  <work>/<mNN>/summary.json  {"status": "ok" | "skipped_no_input",
                              "values": {<every module column>: value or null},
                              "reason": str | null}
  <work>/<mNN>/flags.tsv     (optional) flags in flags.tsv format
A planned module without summary.json has failed.
"""

import csv
import json
from dataclasses import dataclass
from pathlib import Path

from asmqc import flags as fl
from asmqc import schema


@dataclass
class ModuleResult:
    module: str
    status: str
    values: dict
    reason: str | None
    flags: list[fl.Flag]


def collect(work: Path, plan: dict[str, dict]) -> dict[str, ModuleResult]:
    results = {}
    for m in schema.MODULES:
        p = plan[m]
        if p["status"] != "run":
            results[m] = ModuleResult(m, p["status"], {}, p.get("reason"), [])
            continue
        summary = work / m / "summary.json"
        if not summary.exists():
            results[m] = ModuleResult(m, "failed", {}, "module did not finish; see logs/", [])
            continue
        data = json.loads(summary.read_text())
        status = data.get("status", "ok")
        if status not in ("ok", "skipped_no_input"):
            raise ValueError(f"{summary}: invalid status {status!r}")
        values = data.get("values", {})
        expected = set(schema.module_columns(m)) if status == "ok" else set()
        if set(values) != expected:
            raise ValueError(
                f"{summary}: columns differ from schema; "
                f"missing {sorted(expected - set(values))}, extra {sorted(set(values) - expected)}")
        flag_file = work / m / "flags.tsv"
        flags = fl.read(flag_file) if flag_file.exists() else []
        if any(f.module != m for f in flags):
            raise ValueError(f"{flag_file}: flags of another module")
        results[m] = ModuleResult(m, status, values, data.get("reason"), flags)
    return results


def all_flags(results: dict[str, ModuleResult]) -> list[fl.Flag]:
    return [f for m in schema.MODULES for f in results[m].flags]


def summary_row(identity: dict, results: dict[str, ModuleResult]) -> dict[str, str]:
    flags = all_flags(results)
    n = {s: sum(f.severity == s for f in flags) for s in fl.SEVERITIES}
    raw = dict(identity)
    m01_ran = results["m01"].status == "ok"
    raw["ena_rules"] = ("FAIL" if n[fl.ENA_BLOCKING] else "PASS") if m01_ran else None
    raw["n_flags_ena_blocking"] = n[fl.ENA_BLOCKING]
    raw["n_flags_warning"] = n[fl.WARNING]
    raw["n_flags_info"] = n[fl.INFO]
    for m, r in results.items():
        raw[f"{m}_status"] = r.status
        raw.update(r.values)
    return {c: schema.format_value(c, raw.get(c)) for c in schema.HEADER}


def write_summary(path: Path, row: dict[str, str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(schema.HEADER)
        w.writerow([row[c] for c in schema.HEADER])


def read_summary(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    header, body = rows[0], rows[1:]
    return header, [dict(zip(header, r, strict=True)) for r in body]


def wall_seconds(work: Path) -> dict[str, float]:
    """Sum Snakemake benchmark seconds per key (<key>.<rule>.tsv)."""
    out: dict[str, float] = {}
    for f in sorted((work / "benchmarks").glob("*.tsv")):
        key = f.name.split(".", 1)[0]
        with f.open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        out[key] = round(out.get(key, 0.0) + sum(float(r["s"]) for r in rows), 1)
    return out
