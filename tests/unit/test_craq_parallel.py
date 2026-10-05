"""Thread/memory budget and the CRAQ bash shim (parallel LR and SR passes)."""

import os
import subprocess
import time
from pathlib import Path

import pytest

from asmqc import plan as planner
from asmqc import runner

SHIM = Path(__file__).parents[2] / "workflow" / "bin" / "craq_shim"


@pytest.mark.parametrize("cores", [2, 3, 8, 16, 17, 64, 384])
def test_budget_splits_cores_and_memory(cores):
    b = planner.budget(cores, 256 * 1024, craq=True)
    assert b["craq_threads"] >= 1 and b["side_threads"] >= 1
    assert b["craq_threads"] + b["side_threads"] == cores
    assert b["craq_threads"] <= 8
    assert b["craq_mem_mb"] + b["side_mem_mb"] == 256 * 1024


def test_budget_without_craq():
    assert planner.budget(64, 1000, craq=False) == {
        "craq_threads": 0, "side_threads": 64, "craq_mem_mb": 0, "side_mem_mb": 1000}


def fake_passes(d: Path, lr_rc: int = 0) -> None:
    (d / "runLR.sh").write_text(f"sleep 2; date +%s.%N > lr.end; exit {lr_rc}\n")
    (d / "runSR.sh").write_text("date +%s.%N > sr.start; sleep 2; date +%s.%N > sr.end\n")


def driver(d: Path, parallel: int) -> subprocess.CompletedProcess:
    """CRAQ's driver order: system('bash runLR.sh'), system('bash runSR.sh'), then AQI."""
    env = os.environ | {"PATH": f"{SHIM}:{os.environ['PATH']}",
                        "ASMQC_CRAQ_PARALLEL": str(parallel),
                        "ASMQC_CRAQ_LR_STATUS": str(d / "lr.status"), "ASMQC_CRAQ_POLL": "0.2"}
    script = (f"bash {d}/runLR.sh -t 1; bash {d}/runSR.sh -t 1; rc=$?;"
              " test -e lr.end && echo aqi-sees-lr; exit $rc")
    return subprocess.run(["/bin/sh", "-c", script], cwd=d, env=env, capture_output=True, check=False,
                          text=True)


def test_shim_runs_passes_concurrently(tmp_path):
    fake_passes(tmp_path)
    t0 = time.monotonic()
    r = driver(tmp_path, parallel=1)
    elapsed = time.monotonic() - t0
    assert r.returncode == 0 and "aqi-sees-lr" in r.stdout  # AQI waits for both passes
    assert float((tmp_path / "sr.start").read_text()) < float((tmp_path / "lr.end").read_text())
    assert elapsed < 3.8  # serial would take >= 4 s


def test_shim_serial_without_flag(tmp_path):
    fake_passes(tmp_path)
    r = driver(tmp_path, parallel=0)
    assert r.returncode == 0 and "aqi-sees-lr" in r.stdout
    assert float((tmp_path / "sr.start").read_text()) > float((tmp_path / "lr.end").read_text())


def test_shim_reports_long_read_failure(tmp_path):
    fake_passes(tmp_path, lr_rc=3)
    r = driver(tmp_path, parallel=1)
    assert r.returncode == 3 and "long-read pass exit 3" in r.stderr


def test_copy_benchmarks(tmp_path):
    (tmp_path / "work" / "benchmarks").mkdir(parents=True)
    (tmp_path / "work" / "benchmarks" / "m09.craq.tsv").write_text("s\tmax_rss\n1\t2\n")
    runner.copy_benchmarks(tmp_path / "work", tmp_path / "logs" / "benchmarks")
    assert (tmp_path / "logs" / "benchmarks" / "m09.craq.tsv").exists()
    runner.copy_benchmarks(tmp_path / "none", tmp_path / "x")  # no benchmarks: nothing
    assert not (tmp_path / "x").exists()
