"""Helpers every module uses to report to the wrapper (see summary.py)."""

import csv
import json
from pathlib import Path

from asmqc import flags as fl
from asmqc import schema
from asmqc.validate import CHROMOSOMES


def finish(work: Path, module: str, values: dict, flags: list[fl.Flag] = (),
           status: str = "ok", reason: str | None = None) -> None:
    """Write <work>/summary.json and flags.tsv; checks the column set first."""
    expected = set(schema.module_columns(module)) if status == "ok" else set()
    if set(values) != expected:
        raise ValueError(f"{module}: missing {sorted(expected - set(values))}, "
                         f"extra {sorted(set(values) - expected)}")
    for k, v in values.items():
        schema.format_value(k, v)  # fail here, not in the merge
    work.mkdir(parents=True, exist_ok=True)
    fl.write(work / "flags.tsv", list(flags))
    (work / "summary.json").write_text(
        json.dumps({"status": status, "values": values, "reason": reason}, indent=1) + "\n")


def is_chromosome(name: str) -> bool:
    return name in CHROMOSOMES


def write_tsv(path: Path, header: list[str], rows) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def read_fai(path: Path) -> dict[str, int]:
    """Sequence lengths from a samtools .fai, in file order."""
    with path.open() as fh:
        return {f[0]: int(f[1]) for f in (line.split("\t") for line in fh)}
