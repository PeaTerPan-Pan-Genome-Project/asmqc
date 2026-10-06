import json

import pytest

from asmqc import flags as fl
from asmqc import schema, summary

IDENTITY = {"label": "S1", "asmqc_version": "1.0.0", "assembly_md5": "x" * 32,
            "run_date": "2026-10-01"}


def plan(run=("m01",)) -> dict:
    return {m: {"status": "run" if m in run else "skipped_by_user", "reason": None}
            for m in schema.MODULES}


def write_module(work, m, values=None, status="ok", flags=()):
    (work / m).mkdir(parents=True)
    if values is None:
        values = dict.fromkeys(schema.module_columns(m)) if status == "ok" else {}
    (work / m / "summary.json").write_text(json.dumps({"status": status, "values": values}))
    if flags:
        fl.write(work / m / "flags.tsv", list(flags))


def test_statuses_and_na(tmp_path):
    write_module(tmp_path, "m01", flags=[fl.make("m01", "iupac_present")])
    r = summary.collect(tmp_path, plan(run=("m01", "m04")))
    assert r["m01"].status == "ok" and r["m04"].status == "failed"
    row = summary.summary_row(IDENTITY, r)
    assert row["m04_status"] == "failed" and row["m04_capped_arms"] == "NA"
    assert row["m02_status"] == "skipped_by_user"
    assert row["ena_rules"] == "PASS" and row["n_flags_info"] == "1"


def test_ena_fail_iff_blocking_flag(tmp_path):
    write_module(tmp_path, "m01", flags=[fl.make("m01", "seq_lt_20bp", "s")])
    row = summary.summary_row(IDENTITY, summary.collect(tmp_path, plan()))
    assert row["ena_rules"] == "FAIL" and row["n_flags_ena_blocking"] == "1"


def test_ena_na_without_m01(tmp_path):
    row = summary.summary_row(IDENTITY, summary.collect(tmp_path, plan(run=())))
    assert row["ena_rules"] == "NA"


def test_module_skips_itself(tmp_path):
    write_module(tmp_path, "m07", status="skipped_no_input")
    r = summary.collect(tmp_path, plan(run=("m07",)))
    assert r["m07"].status == "skipped_no_input"


def test_missing_column_is_an_error(tmp_path):
    values = dict.fromkeys(schema.module_columns("m01")[1:])
    write_module(tmp_path, "m01", values=values)
    with pytest.raises(ValueError, match="m01_n_seq"):
        summary.collect(tmp_path, plan())


def test_write_read_round_trip(tmp_path):
    write_module(tmp_path, "m01")
    row = summary.summary_row(IDENTITY, summary.collect(tmp_path, plan()))
    summary.write_summary(tmp_path / "qc.tsv", row)
    header, rows = summary.read_summary(tmp_path / "qc.tsv")
    assert header == schema.HEADER and rows == [row]


def test_wall_seconds_concurrent_rules(tmp_path):
    """Two map rules of 100 s each that ran at the same time: 100 s elapsed,
    200 s of rule time; one later rule of another stage."""
    import os

    from asmqc import summary as sm

    b = tmp_path / "benchmarks"
    b.mkdir()
    for name, secs, end in [("map.hifi.tsv", 100.0, 1000.0), ("map.illumina.tsv", 100.0, 1003.0),
                            ("m09.craq.tsv", 50.0, 1100.0)]:
        (b / name).write_text(f"s\th:m:s\n{secs}\t0:01:40\n")
        os.utime(b / name, (end, end))
    wall, rules = sm.wall_seconds(tmp_path)
    assert wall == {"map": 103.0, "m09": 50.0}
    assert rules == {"map": 200.0, "m09": 50.0}
