"""Flag codes, severities and `flags.tsv` I/O (SPEC §7.5)."""

import csv
from dataclasses import astuple, dataclass, fields
from pathlib import Path

ENA_BLOCKING, WARNING, INFO = "ENA_BLOCKING", "WARNING", "INFO"
SEVERITIES = (ENA_BLOCKING, WARNING, INFO)

# code -> (severity, modules that may raise it)
CODES: dict[str, tuple[str, tuple[str, ...]]] = {
    "seq_lt_20bp": (ENA_BLOCKING, ("m01",)),
    "terminal_n": (ENA_BLOCKING, ("m01",)),
    "duplicate_name": (ENA_BLOCKING, ("m01",)),
    "invalid_char": (ENA_BLOCKING, ("m01",)),
    "empty_sequence": (ENA_BLOCKING, ("m01",)),
    "seq_lt_200bp": (WARNING, ("m01",)),
    "n_fraction_gt_50pct": (WARNING, ("m01",)),
    "agp_mismatch": (WARNING, ("m01",)),
    "crlf_line_endings": (WARNING, ("m01",)),
    "iupac_present": (INFO, ("m01",)),
    "inconsistent_line_width": (INFO, ("m01",)),
    "round_length": (INFO, ("m01",)),
    "round_agp_cut": (INFO, ("m01",)),
    "organelle_scaffold": (WARNING, ("m05",)),
    "wrong_orientation_telomere": (INFO, ("m04",)),
    "interstitial_telomere": (INFO, ("m04",)),
    "low_coverage": (WARNING, ("m08", "m09")),
    "reads_not_independent": (INFO, ("m08", "m11")),
}


@dataclass(frozen=True)
class Flag:
    module: str
    code: str
    severity: str
    seq_id: str = ""
    start: int | str = ""
    end: int | str = ""
    value: str = ""
    message: str = ""


COLUMNS = [f.name for f in fields(Flag)]


def make(module: str, code: str, seq_id: str = "", start="", end="", value="",
         message: str = "") -> Flag:
    """Create a flag; the severity comes from CODES, never from the caller."""
    severity, modules = CODES[code]
    if module not in modules:
        raise ValueError(f"flag {code} cannot be raised by {module}")
    return Flag(module, code, severity, seq_id, start, end, str(value), message)


def write(path: Path, flags: list[Flag]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLUMNS)
        w.writerows(astuple(f) for f in flags)


def read(path: Path) -> list[Flag]:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    flags = []
    for r in rows:
        f = make(r["module"], r["code"], r["seq_id"], _int(r["start"]), _int(r["end"]),
                 r["value"], r["message"])
        if f.severity != r["severity"]:
            raise ValueError(f"{path}: severity of {r['code']} must be {f.severity}")
        flags.append(f)
    return flags


def _int(s: str):
    return int(s) if s.lstrip("-").isdigit() else s
