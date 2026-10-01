"""`asmqc test`: built-in smoke test (SPEC §11.1); also the image's %test.

Generates the seeded synthetic data set, runs asmqc on it and compares the
results with the data set's truth. BUSCO is not run (no real genes); the
installed lineage is asserted instead. Long reads are not part of the data
set, so M9 does not run.
"""

import csv
import gzip
import importlib.util
import json
import os
import sys
import tempfile
from pathlib import Path

from asmqc import cli, summary
from asmqc.m06_busco import check_lineage
from asmqc.params import PARAMS
from asmqc.runner import asmqc_home, refs_dir

MODULES = "1,2,4,5,7,8,11"


def _rows(path: Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def check_lineage_installed(refs: Path) -> list[str]:
    lineage = PARAMS["m06"]["lineage"]
    problems = check_lineage(refs / "busco_downloads" / "lineages" / lineage)
    verified = refs / "verified.json"
    table = {line.split("\t")[0]: line.split("\t")[2]
             for line in (asmqc_home() / "refs" / "refs.tsv").read_text().splitlines()[1:]}
    got = json.loads(verified.read_text()) if verified.exists() else {}
    for name, md5 in table.items():
        if got.get(name) != md5:
            problems.append(f"reference {name}: verified md5 {got.get(name)}, expected {md5}")
    return problems


def compare(res: Path, exp: dict) -> list[str]:
    """Problems found comparing a result directory with expected.json."""
    _, (row,) = summary.read_summary(res / "qc_summary.tsv")
    problems = []

    def want(cond, msg):
        if not cond:
            problems.append(msg)

    for key, value in exp["m01"].items():
        if key == "flags_ena_blocking":
            continue
        col = key if key == "ena_rules" else f"m01_{key}"
        ok = (abs(float(row[col]) - value) < 0.006 if isinstance(value, float)
              else row[col] == str(value))
        want(ok, f"{col}: {row[col]}, expected {value}")

    arms = {(r["chromosome"], r["arm"]): r["status"]
            for r in _rows(res / "m04_telomeres" / "telomeres.tsv")}
    for chrom, by_arm in exp["m04"]["arms"].items():
        for arm, a in by_arm.items():
            got = arms.get((chrom, arm))
            want(got == a["status"], f"M4 {chrom} {arm}: {got}, expected {a['status']}")
    want(row["m04_interstitial_arrays"] == str(len(exp["m04"]["interstitial"])),
         f"m04_interstitial_arrays: {row['m04_interstitial_arrays']}")

    org = sorted(r["seq_id"] for r in _rows(res / "m05_organelle_rdna" /
                                            "organelle_scaffolds.tsv")
                 if r["organelle"] != "none")
    want(org == sorted(exp["m05"]["plastid_scaffolds"]), f"M5 organelle scaffolds: {org}")
    for fam in ("45s", "5s"):
        n = exp["m05"][f"rdna{fam}"][0]["copies"]
        want(row[f"m05_rdna{fam}_copies"] == str(n),
             f"m05_rdna{fam}_copies: {row[f'm05_rdna{fam}_copies']}, expected {n}")

    classes = {r["scaffold"]: r["class"] for r in _rows(res / "m07_redundancy" /
                                                         "redundancy.tsv")}
    want(classes == exp["m07"]["classes"], f"M7 classes: {classes}")

    bed = _rows_bed(res / "m11_homopolymer" / "errors.bed.gz")
    hp = sum(any(b[0] == h["seq"] and int(b[1]) + 1 == h["run_start"] and b[3] == "hp"
                 and b[6] == f"+{h['change']}{h['base']}" for b in bed)
             for h in exp["m11"]["hp"])
    st = sum(any(b[0] == t["seq"] and int(b[1]) + 1 == t["run_start"] and b[3] == "str2"
                 and b[6] == f"{'+' if t['change_copies'] > 0 else '-'}"
                             f"{2 * abs(t['change_copies'])}" for b in bed)
             for t in exp["m11"]["str2"])
    want(hp >= 28, f"M11: {hp} of 30 planted HP errors recovered (>= 28 required)")
    want(st >= 4, f"M11: {st} of 5 planted STR2 errors recovered (>= 4 required)")
    return problems


def _rows_bed(path: Path) -> list[list[str]]:
    with gzip.open(path, "rt") as fh:
        return [line.rstrip("\n").split("\t") for line in fh]


def run(directory: Path | None, threads: int) -> int:
    refs = refs_dir()
    problems = check_lineage_installed(refs)
    spec = importlib.util.spec_from_file_location(
        "make_testdata", asmqc_home() / "tests" / "make_testdata.py")
    make_testdata = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(make_testdata)

    with tempfile.TemporaryDirectory(prefix="asmqc_test_", dir=directory) as tmp:
        tmp = Path(tmp)
        data, out = tmp / "data", tmp / "out"
        print("asmqc test: generating the synthetic data set", file=sys.stderr)
        make_testdata.main([str(data), "--refs", str(refs)])
        rc = cli.main(["run", "--assembly", str(data / "asm.fa"), "--agp", str(data / "asm.agp"),
                       "--illumina", f"{data}/reads_R1.fq.gz,{data}/reads_R2.fq.gz",
                       "--reads-used-in-assembly", "no", "--label", "smoke",
                       "--outdir", str(out), "--threads", str(threads), "--mem-gb", "8",
                       "--modules", MODULES])
        if rc != 0:
            problems.append(f"asmqc run exited {rc}; logs in {out / 'smoke' / 'logs'}")
            log = out / "smoke" / "logs" / "snakemake.log"
            if log.exists():
                print(log.read_text()[-3000:], file=sys.stderr)
        else:
            problems += compare(out / "smoke",
                                json.loads((data / "expected.json").read_text()))
    for p in problems:
        print(f"asmqc test: FAIL: {p}", file=sys.stderr)
    print(f"asmqc test: {'FAILED' if problems else 'PASSED'}", file=sys.stderr)
    return 1 if problems else 0


def default_threads() -> int:
    return max(1, min(8, len(os.sched_getaffinity(0))))
