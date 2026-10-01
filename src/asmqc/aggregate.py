"""`asmqc aggregate`: merge result directories (SPEC §4.1, §7.2, §7.4).

Refuses inputs whose qc_summary.tsv headers differ, whose asmqc versions
differ in MAJOR.MINOR (results are comparable only within one minor version,
SPEC principle E), or whose labels repeat.
"""

import csv
import sys
from pathlib import Path

from asmqc import report, schema


class AggregateError(Exception):
    pass


def major_minor(version: str) -> str:
    return ".".join(version.split(".")[:2])


def collect(dirs: list[Path]) -> list[dict]:
    results, problems = [], []
    for d in dirs:
        if not (d / "qc_summary.tsv").exists():
            problems.append(f"{d}: no qc_summary.tsv")
            continue
        results.append(report.load(d))
    if problems:
        raise AggregateError("\n".join(problems))
    if not results:
        raise AggregateError("no result directories given")
    ref = results[0]
    for r in results[1:]:
        if r["header"] != ref["header"]:
            raise AggregateError(f"{r['dir']}: qc_summary.tsv header differs from {ref['dir']}")
    if ref["header"] != schema.HEADER:
        raise AggregateError(f"{ref['dir']}: header differs from this asmqc version's schema")
    versions = {major_minor(r["row"]["asmqc_version"]) for r in results}
    if len(versions) > 1:
        raise AggregateError("results from different asmqc minor versions are not comparable: "
                             + ", ".join(sorted(versions)))
    labels = [r["row"]["label"] for r in results]
    dup = sorted({x for x in labels if labels.count(x) > 1})
    if dup:
        raise AggregateError("duplicate labels: " + ", ".join(dup))
    return results


def run(dirs: list[Path], out: Path) -> int:
    try:
        results = collect(dirs)
    except AggregateError as e:
        print(f"asmqc aggregate: error: {e}", file=sys.stderr)
        return 1
    out.mkdir(parents=True, exist_ok=True)
    with (out / "qc_summary.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(schema.HEADER)
        for r in results:
            w.writerow([r["row"][c] for c in schema.HEADER])
    report.render_combined(results, out)
    print(f"asmqc aggregate: {len(results)} assemblies -> {out}", file=sys.stderr)
    return 0

